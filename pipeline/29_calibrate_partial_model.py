#!/usr/bin/env python3
"""Calibrate the partial-annotation model, and measure whether it helped.

Usage: 29_calibrate_partial_model.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import json
import datetime
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             brier_score_loss, precision_score, recall_score)

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
INNER_FOLDS = 3          # same number of calibration folds as the thesis models
THRESHOLDS = [0.5, 0.867]
X = df[FULL].values.astype(float)
out = []

# Situations: name for the model description, name for the report, groups present
SITUATIONS = [
    ("all four databases, many files", "all four databases, many files",
     ["cells", "dbsnp", "gnomad", "redi", "cosmic"]),
    ("dbSNP+gnomAD, many files", "dbSNP+gnomAD, many files",
     ["cells", "dbsnp", "gnomad"]),
    ("dbSNP+gnomAD, one file", "dbSNP+gnomAD, one file",
     ["dbsnp", "gnomad"]),
    ("gnomAD only, many files", "gnomAD only, many files",
     ["cells", "gnomad"]),
    ("dbSNP only, many files", "dbSNP only, many files",
     ["cells", "dbsnp"]),
    ("gnomAD+REDIportal without dbSNP", "gnomAD+REDIportal without dbSNP",
     ["cells", "gnomad", "redi"]),
]
KEY = "dbSNP+gnomAD, one file"   # the case produced by Ensembl self-annotation


def log(s=""):
    print(s)
    out.append(s)


def augment(idx, gen):
    """Four copies of the rows at idx: the first intact, in every later copy each
optional block withheld with probability one half.
    """
    parts_X, parts_y, parts_g = [X[idx]], [y[idx]], [vid[idx]]
    for _ in range(COPIES - 1):
        Xc = X[idx].copy()
        for cols in OPTIONAL.values():
            drop = gen.random(len(idx)) < 0.5
            for c in cols:
                Xc[drop, FULL.index(c)] = np.nan
        parts_X.append(Xc)
        parts_y.append(y[idx])
        parts_g.append(vid[idx])
    return np.vstack(parts_X), np.concatenate(parts_y), np.concatenate(parts_g)


def calibrated(Xa, ya, ga):
    splits = list(GroupKFold(n_splits=INNER_FOLDS).split(Xa, ya, groups=ga))
    return CalibratedClassifierCV(XGBClassifier(**XGB), method="isotonic",
                                  cv=splits).fit(Xa, ya)


def mask_to(Xm, blocks):
    keep = set(CORE)
    for b in blocks:
        keep |= set(OPTIONAL[b])
    Xm = Xm.copy()
    for j, name in enumerate(FULL):
        if name not in keep:
            Xm[:, j] = np.nan
    return Xm


def ece(yt, p, bins=10):
    """Expected calibration error: average gap between predicted and observed, weighted
by how many variants fall in each score band.
    """
    edges = np.linspace(0, 1, bins + 1)
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, bins - 1)
    e = 0.0
    for b in range(bins):
        m = idx == b
        if m.any():
            e += m.mean() * abs(p[m].mean() - yt[m].mean())
    return e


log("Calibration of the fifth model")
log()
log("rows %d, unique variants %d, RNA_only share %.3f"
    % (len(df), len(set(vid)), y.mean()))
log()

# Honest comparison, out of fold
raw = {n[0]: np.zeros(len(df)) for n in SITUATIONS}
cal = {n[0]: np.zeros(len(df)) for n in SITUATIONS}
rng = np.random.default_rng(42)
for fold, (itr, ite) in enumerate(GroupKFold(n_splits=5).split(X, y, groups=vid), 1):
    Xa, ya, ga = augment(itr, rng)
    m_raw = XGBClassifier(**XGB).fit(Xa, ya)
    m_cal = calibrated(Xa, ya, ga)
    for en, _, blocks in SITUATIONS:
        Xt = mask_to(X[ite], blocks)
        raw[en][ite] = m_raw.predict_proba(Xt)[:, 1]
        cal[en][ite] = m_cal.predict_proba(Xt)[:, 1]
    print("  fold %d of 5 done" % fold)

log("1. Raw against calibrated scores, out of fold")
log()
log("  %-32s %15s %17s %15s" % ("situation", "ROC-AUC", "Brier", "ECE"))
log("  %-32s %15s %17s %15s" % ("", "raw, calibrated", "raw, calibrated", "raw, calibrated"))
metrics = {}
for en, bg, _ in SITUATIONS:
    r, c = raw[en], cal[en]
    m = {
        "roc_auc": round(float(roc_auc_score(y, c)), 4),
        "pr_auc": round(float(average_precision_score(y, c)), 4),
        "brier": round(float(brier_score_loss(y, c)), 4),
        "ece": round(float(ece(y, c)), 4),
    }
    for t in THRESHOLDS:
        pred = (c >= t).astype(int)
        m["precision_at_%s" % t] = round(float(precision_score(y, pred, zero_division=0)), 4)
        m["recall_at_%s" % t] = round(float(recall_score(y, pred, zero_division=0)), 4)
    metrics[en] = m
    log("  %-32s %6.4f -> %6.4f %7.4f -> %6.4f %6.4f -> %6.4f"
        % (bg, roc_auc_score(y, r), m["roc_auc"], brier_score_loss(y, r),
           m["brier"], ece(y, r), m["ece"]))
log()

log("2. Operating points after calibration")
log()
log("For the four main models threshold 0.867 gives precision 0.990.")
log()
log("  %-32s %20s %20s" % ("situation", "threshold 0.5", "threshold 0.867"))
log("  %-32s %20s %20s" % ("", "precision / recall", "precision / recall"))
for en, bg, _ in SITUATIONS:
    m = metrics[en]
    log("  %-32s %9.3f / %-8.3f %9.3f / %-8.3f"
        % (bg, m["precision_at_0.5"], m["recall_at_0.5"],
           m["precision_at_0.867"], m["recall_at_0.867"]))
log()

log("3. Reliability by bin, %s" % [s[1] for s in SITUATIONS if s[0] == KEY][0].upper())
log()
edges = np.linspace(0, 1, 11)
log("  %-11s %22s %22s" % ("bin", "raw: n  pred  obs", "calibrated: n  pred  obs"))
for b in range(10):
    lo, hi = edges[b], edges[b + 1]
    row = "  %4.1f - %3.1f  " % (lo, hi)
    for scores in (raw[KEY], cal[KEY]):
        sel = (scores >= lo) & ((scores < hi) if b < 9 else (scores <= hi))
        if sel.sum():
            row += " %6d %5.2f %6.2f     " % (sel.sum(), scores[sel].mean(), y[sel].mean())
        else:
            row += " %6d %5s %6s     " % (0, "-", "-")
    log(row)
log()

# Final model, all rows
log("4. Final calibrated model on all rows")
log()
Xa, ya, ga = augment(np.arange(len(df)), np.random.default_rng(42))
final = calibrated(Xa, ya, ga)
log("trained on %d rows after copying (%d original, %d copies), %d calibrated submodels"
    % (len(Xa), len(df), COPIES, len(final.calibrated_classifiers_)))
log()

name = "partial_annotation"
bundle = {"model": final, "features": FULL, "imputer": None,
          "handles_missing": True, "positive_class": "RNA_only"}
pkl = os.path.join(outdir, name + ".pkl")
joblib.dump(bundle, pkl)

desc = {
    "name": "partial annotation",
    "use_when": ("some but not all database features present, including after "
                 "self-annotation through Ensembl; one or many files"),
    "algorithm": "XGBClassifier",
    "hyperparameters": {k: v for k, v in XGB.items() if k != "n_jobs"},
    "calibration": "isotonic, cv=%d, folds grouped by variant" % INNER_FOLDS,
    "features": FULL,
    "n_features": len(FULL),
    "optional_blocks": {k: v for k, v in OPTIONAL.items()},
    "missing_values": ("Passed as NaN and never filled. The model was trained "
                       "with feature groups withheld at random, so it has "
                       "learned how to handle their absence."),
    "trained_on": {
        "rows": int(len(df)),
        "unique_variants": int(len(set(vid))),
        "cell_lines": sorted(set(df["group"].astype(str))),
        "cells": int(df["cell"].nunique()),
        "prevalence_rna_only": round(float(y.mean()), 4),
        "augmentation": ("%d copies of the training set; the first intact, in "
                         "each later copy every optional block withheld with "
                         "probability 0.5" % COPIES),
    },
    "metrics_cv_by_variant": metrics,
    "not_for": ("No database features at all. There this model trails the "
                "dedicated 26- and 30-feature models by about 0.017, so those "
                "are used instead."),
    "date": datetime.date.today().isoformat(),
}
js = os.path.join(outdir, name + ".json")
with open(js, "w") as f:
    json.dump(desc, f, ensure_ascii=False, indent=2)

log("saved:")
log("  %s" % pkl)
log("  %s" % js)

path = os.path.join(outdir, "calibrate_partial_model.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
