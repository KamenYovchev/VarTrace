#!/usr/bin/env python3
"""Does the model still work on variants that are in no database?

Usage: 35_check_variants_outside_databases.py <dataset_sig.csv> <outdir>
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

for c in CORE + CROSSCELL + DATABASE:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# A variant no lookup can resolve: absent from every database
novel = ((df["in_dbsnp"].fillna(0) == 0) & (df["in_gnomad"].fillna(0) == 0)
         & (df["in_cosmic"].fillna(0) == 0) & (df["in_rediportal"].fillna(0) == 0)).values
known = ~novel

PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
gkf = GroupKFold(n_splits=5)
out = []


def log(s=""):
    print(s)
    out.append(s)


def oof(feats):
    X = df[feats].values
    p = np.zeros(len(df))
    for itr, ite in gkf.split(X, y, groups=variant_id):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        m = GradientBoostingClassifier(**PARAMS).fit(Xtr, y[itr])
        p[ite] = m.predict_proba(Xte)[:, 1]
    return p


def metrics(yt, p, thr=0.5):
    pred = (p >= thr).astype(int)
    return dict(n=len(yt), prevalence=float(np.mean(yt)),
                roc_auc=roc_auc_score(yt, p),
                pr_auc=average_precision_score(yt, p),
                precision=precision_score(yt, pred, zero_division=0),
                recall=recall_score(yt, pred, zero_division=0),
                f1=f1_score(yt, pred, zero_division=0),
                accuracy=accuracy_score(yt, pred))


log("Does the model work on variants that are in no database?")
log()
log("  all variants : %6d   RNA_only share %.3f" % (len(df), y.mean()))
log("  known        : %6d   RNA_only share %.3f" % (known.sum(), y[known].mean()))
log("  novel        : %6d   RNA_only share %.3f" % (novel.sum(), y[novel].mean()))
log()
log("  novel variants are %.1f%% of all rows" % (100.0 * novel.sum() / len(df)))
log()

SETS = [
    ("data-only (30)", CORE + CROSSCELL),
    ("database-informed (35)", CORE + CROSSCELL + DATABASE),
]

HDR = ("%-26s %-10s %6s %6s %7s %7s %9s %7s %7s" %
       ("model", "subset", "n", "share", "ROC-AUC", "PR-AUC",
        "precision", "recall", "F1"))


def row(name, subset_name, m):
    return ("%-26s %-10s %6d %6.3f %7.3f %7.3f %9.3f %7.3f %7.3f" %
            (name, subset_name, m["n"], m["prevalence"], m["roc_auc"],
             m["pr_auc"], m["precision"], m["recall"], m["f1"]))


log(HDR)
preds = {}
for name, feats in SETS:
    p = oof(feats)
    preds[name] = p
    for subset_name, mask in [("all", np.ones(len(df), bool)),
                              ("known", known), ("novel", novel)]:
        if mask.sum() > 0 and 0 < y[mask].mean() < 1:
            log(row(name, subset_name, metrics(y[mask], p[mask])))
    log()

# Per cell line, novel only
log("Novel variants by cell line (database-informed, out of fold)")
p = preds["database-informed (35)"]
log("%-10s %7s %6s %8s %8s" % ("line", "n", "share", "ROC-AUC", "PR-AUC"))
for ln in ["molm13", "hg001", "hg002"]:
    m = novel & (lines == ln)
    if m.sum() > 10 and 0 < y[m].mean() < 1:
        log("%-10s %7d %6.3f %8.3f %8.3f"
            % (ln, m.sum(), y[m].mean(),
               roc_auc_score(y[m], p[m]), average_precision_score(y[m], p[m])))
log()

# The novel subset cannot settle the question: its RNA_only prevalence is essentially 1,
# so every metric on it is degenerate.
log("Model against database lookup")
log()

lookup_pred = novel.astype(int)          # novel, then RNA_only, known, then DNA_origin
model_pred = (p >= 0.5).astype(int)


def report(name, pred):
    return ("%-22s accuracy %.3f  precision %.3f  recall %.3f  F1 %.3f"
            % (name, accuracy_score(y, pred),
               precision_score(y, pred, zero_division=0),
               recall_score(y, pred, zero_division=0),
               f1_score(y, pred, zero_division=0)))


log(report("database lookup", lookup_pred))
log(report("model (threshold 0.5)", model_pred))
log()

# Where the lookup necessarily fails: in a database, yet RNA-origin
missed = known & (y == 1)
caught = missed & (model_pred == 1)
log("Variants in a database that are RNA-only:")
log("  total %d, the lookup gets all of them wrong" % missed.sum())
log("  the model gets %d of them (%.1f%%)"
    % (caught.sum(), 100.0 * caught.sum() / max(missed.sum(), 1)))
log()

# And the reverse error: not in any database, yet genomic
missed2 = novel & (y == 0)
caught2 = missed2 & (model_pred == 0)
log("Variants in no database that are genomic:")
log("  total %d, the lookup gets all of them wrong" % missed2.sum())
log("  the model gets %d of them" % caught2.sum())
log()

log("Result")
log("Novel variants are almost all one class (RNA_only share %.3f),")
log("so metrics there say little."
    % y[novel].mean())
log()
log("Known variants: RNA_only share %.3f."
    % y[known].mean())
log("On them the model gives ROC-AUC %.3f against a baseline of %.3f."
    % (roc_auc_score(y[known], p[known]), y[known].mean()))

path = os.path.join(outdir, "novel_variants_check.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
