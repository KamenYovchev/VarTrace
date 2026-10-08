#!/usr/bin/env python3
"""Final evaluation for the SNV-origin classifier (thesis results chapter).

Usage: 24_evaluate_final_model.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import (roc_auc_score, roc_curve, accuracy_score,
                             precision_score, recall_score, f1_score)

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"])
y = (df["label"] == "RNA_only").astype(int).values
variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
              + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
SIGNATURES = [c for c in df.columns if c.startswith("sig_")]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
ALL = BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES + DATABASE
for c in ALL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

STEPS = [
    ("baseline", BASELINE),
    ("+read-level", BASELINE + READLEVEL),
    ("+seq-context", BASELINE + READLEVEL + SEQCTX),
    ("+cross-cell", BASELINE + READLEVEL + SEQCTX + CROSSCELL),
    ("+signatures", BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES),
    ("+database", ALL),
]
DATA_ONLY = BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES

# The ablation (26_audit_leakage_and_ablation.py, section 3) showed the 86 signature
# features add nothing measurable: 0.969 without them vs 0.968 with, against a fold sd
# of 0.007, and the same to within 0.001 under leave-one-cell-line-out.
FINAL = BASELINE + READLEVEL + SEQCTX + CROSSCELL + DATABASE
FINAL_DATA_ONLY = BASELINE + READLEVEL + SEQCTX + CROSSCELL


def impute(Xtr, Xte):
    imp = SimpleImputer(strategy="median")
    return imp.fit_transform(Xtr), imp.transform(Xte)


# 1.
gkf = GroupKFold(n_splits=5)
lines = ["5-fold cross-validation, folds split BY VARIANT (GroupKFold)",
         "model: GradientBoosting;  metric: ROC-AUC (mean +- sd over folds)", ""]
cv_means, cv_sds = [], []
for name, feats in STEPS:
    X = df[feats].values
    aucs = []
    for itr, ite in gkf.split(X, y, groups=variant_id):
        Xtr, Xte = impute(X[itr], X[ite])
        m = GradientBoostingClassifier(random_state=42).fit(Xtr, y[itr])
        aucs.append(roc_auc_score(y[ite], m.predict_proba(Xte)[:, 1]))
    mu, sd = float(np.mean(aucs)), float(np.std(aucs))
    cv_means.append(mu)
    cv_sds.append(sd)
    lines.append("%-14s %.3f +- %.3f" % (name, mu, sd))
with open(os.path.join(outdir, "cv_results.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))

# 2.
plt.figure(figsize=(7, 4.2))
xs = np.arange(len(STEPS))
plt.errorbar(xs, cv_means, yerr=cv_sds, marker="o", capsize=4)
plt.xticks(xs, [n for n, _ in STEPS], rotation=20)
plt.ylabel("ROC-AUC (5-fold CV by variant)")
plt.title("Model improvement by feature group")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "progression.png"), dpi=200)
plt.close()

# 3.
gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=42)
itr, ite = next(gss.split(df[FINAL].values, y, groups=variant_id))
plt.figure(figsize=(5.2, 5))
for label_name, feats in [("data-only", FINAL_DATA_ONLY), ("with databases", FINAL)]:
    X = df[feats].values
    Xtr, Xte = impute(X[itr], X[ite])
    m = GradientBoostingClassifier(random_state=42).fit(Xtr, y[itr])
    p = m.predict_proba(Xte)[:, 1]
    fpr, tpr, _ = roc_curve(y[ite], p)
    plt.plot(fpr, tpr, label="%s (AUC=%.3f)" % (label_name, roc_auc_score(y[ite], p)))
plt.plot([0, 1], [0, 1], "k--", lw=0.8)
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.title("ROC, final models (variant split)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "roc.png"), dpi=200)
plt.close()

# 4+6.
X = df[FINAL].values
oof = np.zeros(len(y))
for itr_f, ite_f in gkf.split(X, y, groups=variant_id):
    Xtr, Xte = impute(X[itr_f], X[ite_f])
    cal = CalibratedClassifierCV(GradientBoostingClassifier(random_state=42),
                                 method="isotonic", cv=3)
    cal.fit(Xtr, y[itr_f])
    oof[ite_f] = cal.predict_proba(Xte)[:, 1]

frac_pos, mean_pred = calibration_curve(y, oof, n_bins=10)
plt.figure(figsize=(5.2, 5))
plt.plot(mean_pred, frac_pos, marker="o", label="Variant Origin Score")
plt.plot([0, 1], [0, 1], "k--", lw=0.8, label="perfectly calibrated")
plt.xlabel("Predicted P(RNA-only)")
plt.ylabel("Observed fraction RNA-only")
plt.title("Calibration of the Variant Origin Score")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, "calibration.png"), dpi=200)
plt.close()

out = df[["cell", "group", "chr", "pos", "ref", "alt", "label"]].copy()
out["origin_score_rna_only"] = np.round(oof, 4)
out.to_csv(os.path.join(outdir, "origin_scores.csv"), index=False)
oof_auc = roc_auc_score(y, oof)
print("\nVariant Origin Score: out-of-fold AUC=%.3f  (scores in origin_scores.csv)" % oof_auc)

# 5.
Xtr, Xte = impute(X[itr], X[ite])
try:
    import shap
    m = GradientBoostingClassifier(random_state=42).fit(Xtr, y[itr])
    expl = shap.TreeExplainer(m)
    sv = expl.shap_values(Xte)
    plt.figure()
    shap.summary_plot(sv, features=Xte, feature_names=FINAL, show=False, max_display=15)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "shap_summary.png"), dpi=200, bbox_inches="tight")
    plt.close()
    print("SHAP summary written (shap_summary.png)")
except ImportError:
    rf = RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                random_state=42, n_jobs=-1).fit(Xtr, y[itr])
    imp = sorted(zip(FINAL, rf.feature_importances_), key=lambda t: -t[1])[:15]
    plt.figure(figsize=(6.5, 5))
    plt.barh([f for f, _ in imp][::-1], [v for _, v in imp][::-1])
    plt.xlabel("Feature importance (RandomForest)")
    plt.title("Top 15 features")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "importance.png"), dpi=200)
    plt.close()
    print("shap not installed - wrote importance.png instead")

# Final held-out metrics of calibrated model on the variant split
cal = CalibratedClassifierCV(GradientBoostingClassifier(random_state=42),
                             method="isotonic", cv=3).fit(Xtr, y[itr])
p = cal.predict_proba(Xte)[:, 1]
pred = (p >= 0.5).astype(int)
print("\nfinal calibrated model (variant-split test): "
      "acc=%.3f prec=%.3f rec=%.3f f1=%.3f auc=%.3f" % (
          accuracy_score(y[ite], pred), precision_score(y[ite], pred, zero_division=0),
          recall_score(y[ite], pred, zero_division=0),
          f1_score(y[ite], pred, zero_division=0), roc_auc_score(y[ite], p)))
print("\nall outputs in: %s" % outdir)
