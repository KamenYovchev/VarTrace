#!/usr/bin/env python3
"""One model trained with feature groups withheld at random.

Usage: 37_train_missing_aware_model.py <dataset_sig.csv> <outdir>
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
except ImportError:
    sys.exit("xgboost is missing from the environment")

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"]).reset_index(drop=True)
y = (df["label"] == "RNA_only").astype(int).values
vid = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
       + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DBSNP = ["in_dbsnp"]
GNOMAD = ["gnomad_af", "in_gnomad"]
REDI = ["in_rediportal"]
COSMIC = ["in_cosmic"]

CORE = BASELINE + READLEVEL + SEQCTX
ALLDB = DBSNP + GNOMAD + REDI + COSMIC
FULL = CORE + CROSSCELL + ALLDB
for c in FULL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# The optional blocks, each independently present or absent at prediction time
BLOCKS = {"cells": CROSSCELL, "dbsnp": DBSNP, "gnomad": GNOMAD,
          "redi": REDI, "cosmic": COSMIC}

# Situations to report on.
SITUATIONS = [
    ("nothing, one file", []),
    ("nothing, many files", ["cells"]),
    ("dbSNP+gnomAD, one file", ["dbsnp", "gnomad"]),
    ("dbSNP+gnomAD, many", ["cells", "dbsnp", "gnomad"]),
    ("everything", ["cells", "dbsnp", "gnomad", "redi", "cosmic"]),
    ("gnomAD+REDIportal without dbSNP", ["cells", "gnomad", "redi"]),
    ("dbSNP only", ["cells", "dbsnp"]),
    ("gnomAD only", ["cells", "gnomad"]),
    ("everything but gnomAD", ["cells", "dbsnp", "redi", "cosmic"]),
]

GB = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
XGB = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42,
           eval_metric="logloss", n_jobs=8)
COPIES = 4          # how many masked replicas of the training set
rng = np.random.default_rng(42)
out = []


def log(s=""):
    print(s)
    out.append(s)


def cols_for(blocks):
    keep = set(CORE)
    for b in blocks:
        keep |= set(BLOCKS[b])
    return keep


def mask_to(X, blocks):
    """Return a copy of X with every column outside blocks set to missing."""
    keep = cols_for(blocks)
    Xm = X.copy()
    for j, name in enumerate(FULL):
        if name not in keep:
            Xm[:, j] = np.nan
    return Xm


log("One model trained with missing groups")
log()

Xall = df[FULL].values.astype(float)

# Dedicated models, reference
log("1. Dedicated models, for comparison")
log()
log("One model per situation, trained on exactly its features.")
log()
dedicated = {}
log("  %-32s %6s %9s %9s" % ("situation", "feat", "ROC-AUC", "PR-AUC"))
for name, blocks in SITUATIONS:
    feats = [c for c in FULL if c in cols_for(blocks)]
    X = df[feats].values
    p = np.zeros(len(df))
    for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=vid):
        im = SimpleImputer(strategy="median")
        Xtr, Xte = im.fit_transform(X[itr]), im.transform(X[ite])
        p[ite] = GradientBoostingClassifier(**GB).fit(
            Xtr, y[itr]).predict_proba(Xte)[:, 1]
    dedicated[name] = p
    log("  %-32s %6d %9.4f %9.4f"
        % (name, len(feats), roc_auc_score(y, p), average_precision_score(y, p)))
log()

# The missing-aware model
log("2. The single model trained with gaps")
log()
log("Training set copied %d times, first copy intact, groups withheld in the rest." % COPIES)

shared = {name: np.zeros(len(df)) for name, _ in SITUATIONS}
for itr, ite in GroupKFold(n_splits=5).split(Xall, y, groups=vid):
    parts_X, parts_y = [Xall[itr]], [y[itr]]
    for _ in range(COPIES - 1):
        Xc = Xall[itr].copy()
        for b, cols in BLOCKS.items():
            drop = rng.random(len(itr)) < 0.5
            for c in cols:
                j = FULL.index(c)
                Xc[drop, j] = np.nan
        parts_X.append(Xc)
        parts_y.append(y[itr])
    Xtr = np.vstack(parts_X)
    ytr = np.concatenate(parts_y)
    mdl = XGBClassifier(**XGB).fit(Xtr, ytr)
    for name, blocks in SITUATIONS:
        shared[name][ite] = mdl.predict_proba(mask_to(Xall[ite], blocks))[:, 1]

log("  %-32s %10s %10s %10s" % ("situation", "dedicated", "shared", "difference"))
diffs = []
for name, _ in SITUATIONS:
    a = roc_auc_score(y, dedicated[name])
    b = roc_auc_score(y, shared[name])
    diffs.append(b - a)
    log("  %-32s %10.4f %10.4f %+10.4f" % (name, a, b, b - a))
log()
worst = min(diffs)
log("  worst gap: %+.4f ROC-AUC" % worst)
log("  mean difference: %+.4f" % float(np.mean(diffs)))
log()

log("Result")
if worst > -0.005:
    log("One model is enough.")
elif worst > -0.02:
    log("The single model falls slightly behind.")
else:
    log("The single model falls clearly behind, keep the dedicated models.")

path = os.path.join(outdir, "missing_aware.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
