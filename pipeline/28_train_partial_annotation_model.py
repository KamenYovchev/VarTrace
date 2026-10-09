#!/usr/bin/env python3
"""Train the fifth VarTrace model, the one for partial annotation.

Usage: 28_train_partial_annotation_model.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import json
import datetime
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
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
DB_BLOCKS = {
    "dbsnp": ["in_dbsnp"],
    "gnomad": ["gnomad_af", "in_gnomad"],
    "redi": ["in_rediportal"],
    "cosmic": ["in_cosmic"],
}
CORE = BASELINE + READLEVEL + SEQCTX
OPTIONAL = dict(DB_BLOCKS)
OPTIONAL["cells"] = CROSSCELL
FULL = CORE + CROSSCELL + [c for b in DB_BLOCKS.values() for c in b]
for c in FULL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

XGB = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42,
           eval_metric="logloss", n_jobs=8)
COPIES = 4
rng = np.random.default_rng(42)
X = df[FULL].values.astype(float)
out = []


def log(s=""):
    print(s)
    out.append(s)


def augment(idx, gen):
    """Replicate the rows at idx, withholding optional blocks at random."""
    parts_X, parts_y = [X[idx]], [y[idx]]
    for _ in range(COPIES - 1):
        Xc = X[idx].copy()
        for cols in OPTIONAL.values():
            drop = gen.random(len(idx)) < 0.5
            for c in cols:
                Xc[drop, FULL.index(c)] = np.nan
        parts_X.append(Xc)
        parts_y.append(y[idx])
    return np.vstack(parts_X), np.concatenate(parts_y)


log("Fifth model, for partial annotation")
log()
log("rows %d, unique variants %d, RNA_only share %.3f"
    % (len(df), len(set(vid)), y.mean()))
log("features %d, optional groups %d" % (len(FULL), len(OPTIONAL)))
log()

# Honest metrics, out of fold
log("Out-of-fold metrics")
log()

SITUATIONS = [
    ("all four databases, many files", ["cells", "dbsnp", "gnomad", "redi", "cosmic"]),
    ("dbSNP+gnomAD, many files", ["cells", "dbsnp", "gnomad"]),
    ("dbSNP+gnomAD, one file", ["dbsnp", "gnomad"]),
    ("gnomAD only, many files", ["cells", "gnomad"]),
    ("dbSNP only, many files", ["cells", "dbsnp"]),
    ("gnomAD+REDIportal without dbSNP", ["cells", "gnomad", "redi"]),
]


def keep_cols(blocks):
    keep = set(CORE)
    for b in blocks:
        keep |= set(OPTIONAL[b])
    return keep


scores = {n: np.zeros(len(df)) for n, _ in SITUATIONS}
for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=vid):
    Xtr, ytr = augment(itr, rng)
    mdl = XGBClassifier(**XGB).fit(Xtr, ytr)
    for name, blocks in SITUATIONS:
        keep = keep_cols(blocks)
        Xte = X[ite].copy()
        for j, col in enumerate(FULL):
            if col not in keep:
                Xte[:, j] = np.nan
        scores[name][ite] = mdl.predict_proba(Xte)[:, 1]

log("  %-34s %9s %9s" % ("situation", "ROC-AUC", "PR-AUC"))
metrics = {}
for name, _ in SITUATIONS:
    r = roc_auc_score(y, scores[name])
    a = average_precision_score(y, scores[name])
    metrics[name] = {"roc_auc": round(float(r), 4), "pr_auc": round(float(a), 4)}
    log("  %-34s %9.4f %9.4f" % (name, r, a))
log()

# Final model, all rows
log("Final model on all rows")
log()
Xtr, ytr = augment(np.arange(len(df)), np.random.default_rng(42))
final = XGBClassifier(**XGB).fit(Xtr, ytr)
log("trained on %d rows after copying (%d original, %d copies)"
    % (len(Xtr), len(df), COPIES))
log()

name = "partial_annotation"
bundle = {"model": final, "features": FULL, "imputer": None,
          "handles_missing": True}
pkl = os.path.join(outdir, name + ".pkl")
joblib.dump(bundle, pkl)

desc = {
    "name": name,
    "purpose": ("Partial annotation: any subset of the database features, "
                "with or without the cross-cell features."),
    "algorithm": "XGBClassifier",
    "hyperparameters": {k: v for k, v in XGB.items() if k != "n_jobs"},
    "features": FULL,
    "optional_blocks": {k: v for k, v in OPTIONAL.items()},
    "missing_values": ("Passed as NaN and never filled; the model was trained with "
                       "feature groups withheld at random."),
    "training": {
        "rows": int(len(df)),
        "unique_variants": int(len(set(vid))),
        "positive_rate": round(float(y.mean()), 4),
        "augmentation": ("%d copies of the training set; the first intact, in the "
                         "others each optional group withheld with "
                         "probability 0.5" % COPIES),
        "cell_lines": sorted(set(df["group"].astype(str))),
        "cells": int(df["cell"].nunique()),
    },
    "metrics_out_of_fold": metrics,
    "not_for": ("No database features at all; there it trails by about "
                "0.017, so the 26- and 30-feature models are used."),
    "date": datetime.date.today().isoformat(),
}
js = os.path.join(outdir, name + ".json")
with open(js, "w") as f:
    json.dump(desc, f, ensure_ascii=False, indent=2)

log("saved:")
log("  %s" % pkl)
log("  %s" % js)

path = os.path.join(outdir, "train_partial_model.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
