#!/usr/bin/env python3
"""Is real ADAR editing present in the RNA_only class?

Usage: 31_split_rna_class_by_rediportal.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, average_precision_score

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"]).reset_index(drop=True)
y = (df["label"] == "RNA_only").astype(int).values
vid = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
       + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
sub = (df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
redi = pd.to_numeric(df["in_rediportal"], errors="coerce").fillna(0).values > 0

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
FEATS = BASELINE + READLEVEL + SEQCTX + CROSSCELL + DATABASE
for c in FEATS:
    df[c] = pd.to_numeric(df[c], errors="coerce")

out = []


def log(s=""):
    print(s)
    out.append(s)


TYPES = sorted(set(sub))


def spectrum(mask, name):
    n = mask.sum()
    if n == 0:
        return
    cnt = pd.Series(sub[mask]).value_counts()
    log("  %-34s %6d " % (name, n)
        + " ".join("%5.1f" % (100.0 * cnt.get(t, 0) / n) for t in TYPES))


log("Is there real editing in the RNA class")
log()

rna = y == 1
dna = y == 0
log("RNA-only variants: %d" % rna.sum())
log("  in REDIportal: %d (%.1f%%)"
    % ((rna & redi).sum(), 100.0 * (rna & redi).sum() / rna.sum()))
log("  not in REDIportal: %d (%.1f%%)"
    % ((rna & ~redi).sum(), 100.0 * (rna & ~redi).sum() / rna.sum()))
log()
log("genomic variants: %d, in REDIportal: %d (%.1f%%)"
    % (dna.sum(), (dna & redi).sum(), 100.0 * (dna & redi).sum() / dna.sum()))
log()

log("Substitution spectrum")
log()
log("  %-34s %6s " % ("", "n") + " ".join("%5s" % t for t in TYPES))
spectrum(rna & redi, "RNA only, in REDIportal")
spectrum(rna & ~redi, "RNA only, not in REDIportal")
spectrum(dna & redi, "genomic, in REDIportal")
spectrum(dna & ~redi, "genomic, not in REDIportal")
log()

log("ADAR signature, A>G plus T>C:")
for name, m in [("RNA only, in REDIportal", rna & redi),
                ("RNA only, not in REDIportal", rna & ~redi),
                ("genomic, in REDIportal", dna & redi),
                ("genomic, not in REDIportal", dna & ~redi)]:
    if m.sum():
        ag = ((sub == "A>G") & m).sum()
        tc = ((sub == "T>C") & m).sum()
        log("  %-34s %5.1f%%   (A>G %.1f%%, T>C %.1f%%)"
            % (name, 100.0 * (ag + tc) / m.sum(),
               100.0 * ag / m.sum(), 100.0 * tc / m.sum()))
log()

# How well does the model do on each half
log("Model on the two halves")
log()
X = df[FEATS].values
p = np.zeros(len(df))
for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=vid):
    imp = SimpleImputer(strategy="median")
    Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
    p[ite] = GradientBoostingClassifier(
        learning_rate=0.05, max_depth=3, n_estimators=300,
        random_state=42).fit(Xtr, y[itr]).predict_proba(Xte)[:, 1]

log("  %-24s %7s %7s %9s %9s" % ("subset", "n", "share", "ROC-AUC", "PR-AUC"))
for name, m in [("in REDIportal", redi), ("not in REDIportal", ~redi)]:
    if m.sum() > 20 and 0 < y[m].mean() < 1:
        log("  %-24s %7d %7.3f %9.3f %9.3f"
            % (name, m.sum(), y[m].mean(),
               roc_auc_score(y[m], p[m]), average_precision_score(y[m], p[m])))
log()
log("mean model score:")
for name, m in [("RNA only, in REDIportal", rna & redi),
                ("RNA only, not in REDIportal", rna & ~redi),
                ("genomic, in REDIportal", dna & redi),
                ("genomic, not in REDIportal", dna & ~redi)]:
    if m.sum():
        log("  %-34s %.3f" % (name, p[m].mean()))
log()

log("Result")
in_redi = 100.0 * (((sub == "A>G") | (sub == "T>C")) & rna & redi).sum() / max((rna & redi).sum(), 1)
out_redi = 100.0 * (((sub == "A>G") | (sub == "T>C")) & rna & ~redi).sum() / max((rna & ~redi).sum(), 1)
log("A>G + T>C, RNA variants in REDIportal:     %.1f%%" % in_redi)
log("A>G + T>C, RNA variants not in REDIportal: %.1f%%" % out_redi)
log("difference: %+.1f points" % (in_redi - out_redi))
log()
if in_redi - out_redi > 15:
    log("Catalogued editing sites carry the editing signature, the rest of the class does not.")
else:
    log("Even catalogued editing sites show no clear signature here.")

path = os.path.join(outdir, "rediportal_split.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
