#!/usr/bin/env python3
"""How should VarTrace handle database features a user may or may not have?

Usage: 33_measure_database_tiers.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, average_precision_score

try:
    from xgboost import XGBClassifier
    HAVE_XGB = True
except ImportError:
    HAVE_XGB = False

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"]).reset_index(drop=True)
y = (df["label"] == "RNA_only").astype(int).values
vid = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
       + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
line = df["group"].values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]

DB_TIERS = {
    "none": [],
    "free_api": ["in_dbsnp", "gnomad_af", "in_gnomad"],
    "free_full": ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal"],
    "all": ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"],
}
TIER_ORDER = ["none", "free_api", "free_full", "all"]
CORE = BASELINE + READLEVEL + SEQCTX
FULL = CORE + CROSSCELL + DB_TIERS["all"]
for c in FULL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

GB = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
XGB = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42,
           eval_metric="logloss", n_jobs=8)
out = []


def log(s=""):
    print(s)
    out.append(s)


def oof_dedicated(feats):
    """A model trained on exactly these features, nothing else."""
    X = df[feats].values
    p = np.zeros(len(df))
    for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=vid):
        im = SimpleImputer(strategy="median")
        Xtr, Xte = im.fit_transform(X[itr]), im.transform(X[ite])
        p[ite] = GradientBoostingClassifier(**GB).fit(
            Xtr, y[itr]).predict_proba(Xte)[:, 1]
    return p


log("Missing database features")
log()

# 1.
redi = pd.to_numeric(df["in_rediportal"], errors="coerce").fillna(0).values > 0
log("1. How often the REDIportal feature is set")
log()
log("  %-12s %8s %12s %8s" % ("line", "rows", "in REDIportal", "share"))
for ln in sorted(set(line)):
    m = line == ln
    log("  %-12s %8d %12d %7.2f%%"
        % (ln, m.sum(), (redi & m).sum(), 100.0 * (redi & m).mean()))
log("  %-12s %8d %12d %7.2f%%" % ("total", len(df), redi.sum(), 100.0 * redi.mean()))
log()
log("by class:")
log("  genomic     %5d of %5d  (%.2f%%)"
    % ((redi & (y == 0)).sum(), (y == 0).sum(), 100.0 * redi[y == 0].mean()))
log("  RNA only    %5d of %5d  (%.2f%%)"
    % ((redi & (y == 1)).sum(), (y == 1).sum(), 100.0 * redi[y == 1].mean()))
log()

# 2.
log("2. Eight dedicated models, one per situation")
log()
log("  %-14s %-8s %6s %9s %9s" % ("tier", "cells", "feat", "ROC-AUC", "PR-AUC"))
dedicated = {}
for tier in TIER_ORDER:
    for cells, extra in [("yes", CROSSCELL), ("no", [])]:
        feats = CORE + extra + DB_TIERS[tier]
        p = oof_dedicated(feats)
        dedicated[(tier, cells)] = (p, feats)
        log("  %-14s %-8s %6d %9.4f %9.4f"
            % (tier, cells, len(feats), roc_auc_score(y, p),
               average_precision_score(y, p)))
log()
log("cost of each step down, with cross-cell features:")
for a, b in zip(TIER_ORDER[1:], TIER_ORDER[:-1]):
    ra = roc_auc_score(y, dedicated[(a, "yes")][0])
    rb = roc_auc_score(y, dedicated[(b, "yes")][0])
    log("  %-10s -> %-10s %+.4f" % (a, b, rb - ra))
log()

# 3.
log("3. One shared model with missing values")
log()
if not HAVE_XGB:
    log("xgboost not available, section skipped.")
else:
    Xall = df[FULL].values.astype(float)
    # one shared OOF loop; inside each fold predict once per tier
    preds = {k: np.zeros(len(df)) for k in dedicated}
    for itr, ite in GroupKFold(n_splits=5).split(Xall, y, groups=vid):
        mdl = XGBClassifier(**XGB).fit(Xall[itr], y[itr])
        for (tier, cells) in dedicated:
            keep = set(CORE + (CROSSCELL if cells == "yes" else []) + DB_TIERS[tier])
            Xte = Xall[ite].copy()
            for j, name in enumerate(FULL):
                if name not in keep:
                    Xte[:, j] = np.nan
            preds[(tier, cells)][ite] = mdl.predict_proba(Xte)[:, 1]
    log("  %-14s %-8s %10s %10s %10s"
        % ("tier", "cells", "dedicated", "shared model", "difference"))
    worst = 0.0
    for tier in TIER_ORDER:
        for cells in ("yes", "no"):
            a = roc_auc_score(y, dedicated[(tier, cells)][0])
            b = roc_auc_score(y, preds[(tier, cells)])
            worst = min(worst, b - a)
            log("  %-14s %-8s %10.4f %10.4f %+10.4f" % (tier, cells, a, b, b - a))
    log()
    log("  worst gap of the shared model: %+.4f ROC-AUC" % worst)
    log()
    if worst > -0.005:
        log("One shared model is enough.")
    else:
        log("The shared model falls behind, keep the dedicated models.")

path = os.path.join(outdir, "database_tiers.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
