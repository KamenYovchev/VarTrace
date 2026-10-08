#!/usr/bin/env python3
"""Five closing analyses, all from the same out-of-fold predictions.

Usage: 30_analyse_errors_and_thresholds.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             precision_score, recall_score, f1_score,
                             accuracy_score, precision_recall_curve)

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"]).reset_index(drop=True)
y = (df["label"] == "RNA_only").astype(int).values
variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
              + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
FEATS = BASELINE + READLEVEL + SEQCTX + CROSSCELL + DATABASE
for c in FEATS:
    df[c] = pd.to_numeric(df[c], errors="coerce")

PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
out = []


def log(s=""):
    print(s)
    out.append(s)


# Out-of-fold predictions, the basis for everything below
X = df[FEATS].values
p = np.zeros(len(df))
for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=variant_id):
    imp = SimpleImputer(strategy="median")
    Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
    p[ite] = GradientBoostingClassifier(**PARAMS).fit(Xtr, y[itr]).predict_proba(Xte)[:, 1]

pred = (p >= 0.5).astype(int)
wrong = pred != y
sub = df["ref"].astype(str) + ">" + df["alt"].astype(str)
in_any_db = ((df["in_dbsnp"].fillna(0) > 0) | (df["in_gnomad"].fillna(0) > 0)
             | (df["in_cosmic"].fillna(0) > 0) | (df["in_rediportal"].fillna(0) > 0))

log("Baseline numbers")
log("rows %d   unique variants %d   RNA_only share %.3f"
    % (len(df), len(set(variant_id)), y.mean()))
log("ROC-AUC %.3f   PR-AUC %.3f   errors at threshold 0.5: %d (%.1f%%)"
    % (roc_auc_score(y, p), average_precision_score(y, p),
       wrong.sum(), 100.0 * wrong.mean()))
log()

# 1.
log("1. Where the errors are")
log()
log("FP: called RNA but in the genome.  FN: called genomic but RNA.")
log()


def breakdown(title, keys, order=None):
    log(title)
    log("  %-26s %7s %8s %7s %7s" % ("", "n", "errors", "FP", "FN"))
    grp = pd.DataFrame({"k": keys, "wrong": wrong, "y": y, "pred": pred})
    it = order if order is not None else sorted(grp["k"].unique(), key=str)
    for k in it:
        m = grp["k"] == k
        if m.sum() < 20:
            continue
        fp = ((grp["pred"] == 1) & (grp["y"] == 0) & m).sum()
        fn = ((grp["pred"] == 0) & (grp["y"] == 1) & m).sum()
        log("  %-26s %7d %7.1f%% %7d %7d"
            % (str(k), m.sum(), 100.0 * grp.loc[m, "wrong"].mean(), fp, fn))
    log()


def qbin(col, n=4):
    v = pd.to_numeric(df[col], errors="coerce")
    try:
        return pd.qcut(v, n, duplicates="drop").astype(str)
    except ValueError:
        return v.astype(str)


breakdown("by cell line", df["group"].values)
breakdown("by depth (quartiles)", qbin("dp").values)
breakdown("by allele fraction (quartiles)", qbin("vaf").values)
breakdown("by number of cells with the variant (quartiles)", qbin("n_cells_variant").values)
breakdown("by database presence",
          np.where(in_any_db, "in at least one database", "in no database"))
breakdown("by substitution type", sub.values)

# 2.
log("2. Operating points")
log()
prec, rec, thr = precision_recall_curve(y, p)
log("  mode                 threshold  precision   recall")
for target in [0.99, 0.95, 0.90]:
    ok = np.where(prec[:-1] >= target)[0]
    if len(ok):
        i = ok[np.argmax(rec[:-1][ok])]
        log("  precision >= %.2f        %.3f      %.3f    %.3f"
            % (target, thr[i], prec[i], rec[i]))
for target in [0.99, 0.95, 0.90]:
    ok = np.where(rec[:-1] >= target)[0]
    if len(ok):
        i = ok[np.argmax(prec[:-1][ok])]
        log("  recall    >= %.2f        %.3f      %.3f    %.3f"
            % (target, thr[i], prec[i], rec[i]))
log()
log("  at threshold 0.50        0.500      %.3f    %.3f"
    % (precision_score(y, pred, zero_division=0),
       recall_score(y, pred, zero_division=0)))
log()

# 3.
log("3. Model against common filtering rules")
log()

vaf = pd.to_numeric(df["vaf"], errors="coerce").fillna(0).values
dp = pd.to_numeric(df["dp"], errors="coerce").fillna(0).values
redi = df["in_rediportal"].fillna(0).values > 0

RULES = [
    ("R1 outside all databases", (~in_any_db).values.astype(int)),
    ("R2 R1 or in REDIportal", ((~in_any_db).values | redi).astype(int)),
    ("R3 R2 or VAF < 0.25", ((~in_any_db).values | redi | (vaf < 0.25)).astype(int)),
    ("R4 R3 and depth >= 10",
     (((~in_any_db).values | redi | (vaf < 0.25)) & (dp >= 10)).astype(int)),
]
log("  %-26s %9s %11s %8s %7s" % ("", "accuracy", "precision", "recall", "F1"))
for name, r in RULES + [("model (threshold 0.5)", pred)]:
    log("  %-26s %9.3f %11.3f %8.3f %7.3f"
        % (name, accuracy_score(y, r), precision_score(y, r, zero_division=0),
           recall_score(y, r, zero_division=0), f1_score(y, r, zero_division=0)))
log()

# 4.
log("4. Substitution spectrum")
log()

groups = [("true RNA", y == 1), ("true genome", y == 0),
          ("predicted RNA", pred == 1), ("predicted genome", pred == 0)]
types = sorted(sub.unique())
log("  %-24s %6s " % ("", "n") + " ".join("%5s" % t for t in types))
for name, m in groups:
    counts = sub[m].value_counts()
    tot = m.sum()
    log("  %-24s %6d " % (name, tot)
        + " ".join("%5.1f" % (100.0 * counts.get(t, 0) / tot) for t in types))
log()
for name, m in groups:
    ag = ((sub == "A>G") & m).sum()
    tc = ((sub == "T>C") & m).sum()
    log("  %-24s A>G + T>C = %.1f%%" % (name, 100.0 * (ag + tc) / m.sum()))
log()

# 5.
log("5. Protocol versions inside HG001")
log()
hg1 = df["group"].values == "hg001"
v21 = hg1 & df["cell"].astype(str).str.startswith("ResolveOMEv2.1").values
v10 = hg1 & ~df["cell"].astype(str).str.startswith("ResolveOMEv2.1").values
log("  %-14s %7s %7s %6s %8s %8s" % ("", "cells", "n", "share", "ROC-AUC", "PR-AUC"))
for name, m in [("v2.1", v21), ("v1.0", v10), ("both together", hg1)]:
    if m.sum() > 20 and 0 < y[m].mean() < 1:
        log("  %-14s %7d %7d %6.3f %8.3f %8.3f"
            % (name, df.loc[m, "cell"].nunique(), m.sum(), y[m].mean(),
               roc_auc_score(y[m], p[m]), average_precision_score(y[m], p[m])))
log()
log("  variants per cell: v2.1 %.0f, v1.0 %.0f"
    % (v21.sum() / max(df.loc[v21, "cell"].nunique(), 1),
       v10.sum() / max(df.loc[v10, "cell"].nunique(), 1)))
log()

path = os.path.join(outdir, "analysis_extras.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
