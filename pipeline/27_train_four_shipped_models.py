#!/usr/bin/env python3
"""Train and evaluate the four model variants that VarTrace ships.

Usage: 27_train_four_shipped_models.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import json
import datetime
import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.impute import SimpleImputer
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             precision_score, recall_score, f1_score,
                             accuracy_score)

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"])
y = (df["label"] == "RNA_only").astype(int).values
variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
              + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
lines = df["group"].values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]

CORE = BASELINE + READLEVEL + SEQCTX
VARIANTS = [
    ("data-only, no cross-cell",         CORE,                          "single VCF"),
    ("data-only, with cross-cell",       CORE + CROSSCELL,              "many VCFs"),
    ("database-informed, no cross-cell", CORE + DATABASE,               "single annotated"),
    ("database-informed, with cross-cell", CORE + CROSSCELL + DATABASE, "many annotated"),
]

for c in CORE + CROSSCELL + DATABASE:
    df[c] = pd.to_numeric(df[c], errors="coerce")

PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
LINES = ["molm13", "hg001", "hg002"]
gkf = GroupKFold(n_splits=5)
out = []


def log(s=""):
    print(s)
    out.append(s)


def model():
    return GradientBoostingClassifier(**PARAMS)


def metrics(yt, p, thr=0.5):
    pred = (p >= thr).astype(int)
    return dict(roc_auc=roc_auc_score(yt, p),
                pr_auc=average_precision_score(yt, p),
                precision=precision_score(yt, pred, zero_division=0),
                recall=recall_score(yt, pred, zero_division=0),
                f1=f1_score(yt, pred, zero_division=0),
                accuracy=accuracy_score(yt, pred),
                prevalence=float(np.mean(yt)))


def oof(feats):
    """Out-of-fold predictions over the variant-grouped folds."""
    X = df[feats].values
    p = np.zeros(len(df))
    for itr, ite in gkf.split(X, y, groups=variant_id):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        p[ite] = model().fit(Xtr, y[itr]).predict_proba(Xte)[:, 1]
    return p


def logo(feats, held):
    X = df[feats].values
    te = lines == held
    imp = SimpleImputer(strategy="median")
    Xtr, Xte = imp.fit_transform(X[~te]), imp.transform(X[te])
    return y[te], model().fit(Xtr, y[~te]).predict_proba(Xte)[:, 1]


HDR = ("%-38s %7s %7s %9s %7s %7s %8s" %
       ("", "ROC-AUC", "PR-AUC", "precision", "recall", "F1", "accuracy"))


def row(lbl, m):
    return ("%-38s %7.3f %7.3f %9.3f %7.3f %7.3f %8.3f" %
            (lbl, m["roc_auc"], m["pr_auc"], m["precision"],
             m["recall"], m["f1"], m["accuracy"]))


log("The four model variants")
log("rows=%d  unique variants=%d  prevalence=%.3f"
    % (len(df), len(set(variant_id)), y.mean()))
log()
log("1. cross-validation split BY variant (pooled out-of-fold)")
log(HDR)
scores = {}
for name, feats, _ in VARIANTS:
    scores[name] = oof(feats)
    log(row("%s (%d)" % (name, len(feats)), metrics(y, scores[name])))
log()
log("2. leave-one-cell-line-out (ROC-AUC / PR-AUC)")
log("%-38s %16s %16s %16s" % ("", "molm13", "hg001", "hg002"))
for name, feats, _ in VARIANTS:
    cells = []
    for held in LINES:
        yt, p = logo(feats, held)
        cells.append("%.3f / %.3f" % (roc_auc_score(yt, p),
                                      average_precision_score(yt, p)))
    log("%-38s %16s %16s %16s" % ("%s (%d)" % (name, len(feats)), *cells))
log()

# Imputation

FULL = CORE + CROSSCELL + DATABASE
DEDICATED = CORE + DATABASE
Xfull = df[FULL].values
p_imp = np.zeros(len(df))
cc_idx = [FULL.index(c) for c in CROSSCELL]
for itr, ite in gkf.split(Xfull, y, groups=variant_id):
    imp = SimpleImputer(strategy="median")
    Xtr = imp.fit_transform(Xfull[itr])
    m = model().fit(Xtr, y[itr])
    # the test rows arrive with the cross-cell columns undefined
    Xte_blank = Xfull[ite].copy()
    Xte_blank[:, cc_idx] = np.nan
    Xte = imp.transform(Xte_blank)          # filled with the TRAINING median
    p_imp[ite] = m.predict_proba(Xte)[:, 1]

m_imp = metrics(y, p_imp)
m_ded = metrics(y, scores["database-informed, no cross-cell"])
m_full = metrics(y, scores["database-informed, with cross-cell"])

log(HDR)
log(row("full 35, cross-cell imputed", m_imp))
log(row("dedicated 31", m_ded))
log(row("full 35, cross-cell present", m_full))
log()
d = m_ded["roc_auc"] - m_imp["roc_auc"]
log("dedicated minus imputed: %+.4f ROC-AUC, %+.4f PR-AUC"
    % (d, m_ded["pr_auc"] - m_imp["pr_auc"]))
log("cost of losing the cross-cell features altogether: %+.4f ROC-AUC"
    % (m_ded["roc_auc"] - m_full["roc_auc"]))
log()
log("verdict: %s" % ("dedicated models are worth shipping - imputing loses more "
                     "than a thousandth" if d > 0.001 else
                     "imputing is as good as a dedicated model on this data"))
log()

# Save models
log("4. saved models")
today = datetime.date.today().isoformat()
for name, feats, usage in VARIANTS:
    X = df[feats].values
    imp = SimpleImputer(strategy="median")
    Xi = imp.fit_transform(X)
    cal = CalibratedClassifierCV(model(), method="isotonic", cv=3).fit(Xi, y)
    slug = name.replace(", ", "__").replace(" ", "-").replace("-", "_")
    pkl = os.path.join(outdir, "%s.pkl" % slug)
    joblib.dump({"model": cal, "imputer": imp, "features": feats,
                 "positive_class": "RNA_only"}, pkl)
    meta = {
        "name": name,
        "use_when": usage,
        "algorithm": "GradientBoostingClassifier",
        "hyperparameters": {k: v for k, v in PARAMS.items() if k != "random_state"},
        "calibration": "isotonic, cv=3",
        "features": feats,
        "n_features": len(feats),
        "trained_on": {"rows": int(len(df)),
                       "unique_variants": int(len(set(variant_id))),
                       "cell_lines": sorted(set(lines.tolist())),
                       "prevalence_rna_only": round(float(y.mean()), 4)},
        "metrics_cv_by_variant": {k: round(v, 4)
                                  for k, v in metrics(y, scores[name]).items()},
        "date": today,
    }
    with open(os.path.join(outdir, "%s.json" % slug), "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    log("  %-40s %3d features -> %s" % (name, len(feats), os.path.basename(pkl)))

path = os.path.join(outdir, "model_variants.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
