#!/usr/bin/env python3
"""How much does the in_cosmic feature carry, and what does a user lose without it?

Usage: 32_measure_cosmic_contribution.py <dataset_sig.csv> <outdir>
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
line = df["group"].values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
FULL = BASELINE + READLEVEL + SEQCTX + CROSSCELL + DATABASE
NOCOSMIC = [f for f in FULL if f != "in_cosmic"]
NOCELL = BASELINE + READLEVEL + SEQCTX + DATABASE
NOCELL_NOCOSMIC = [f for f in NOCELL if f != "in_cosmic"]
for c in FULL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
CANCER = "molm13"
out = []


def log(s=""):
    print(s)
    out.append(s)


def oof(feats):
    X = df[feats].values
    p = np.zeros(len(df))
    for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=vid):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        p[ite] = GradientBoostingClassifier(**PARAMS).fit(
            Xtr, y[itr]).predict_proba(Xte)[:, 1]
    return p


log("What the COSMIC feature adds")
log()

# 1.
cos = pd.to_numeric(df["in_cosmic"], errors="coerce").fillna(0).values > 0
dbsnp = pd.to_numeric(df["in_dbsnp"], errors="coerce").fillna(0).values > 0

log("1. How often it is set")
log()
log("  %-12s %8s %10s %10s" % ("line", "rows", "in COSMIC", "share"))
for ln in sorted(set(line)):
    m = line == ln
    log("  %-12s %8d %10d %9.1f%%"
        % (ln, m.sum(), (cos & m).sum(), 100.0 * (cos & m).mean()))
log("  %-12s %8d %10d %9.1f%%" % ("total", len(df), cos.sum(), 100.0 * cos.mean()))
log()
log("by class:")
log("  genomic     %6d of %6d  (%.1f%%)"
    % ((cos & (y == 0)).sum(), (y == 0).sum(), 100.0 * cos[y == 0].mean()))
log("  RNA only    %6d of %6d  (%.1f%%)"
    % ((cos & (y == 1)).sum(), (y == 1).sum(), 100.0 * cos[y == 1].mean()))
log()
log("overlap with dbSNP:")
if cos.sum():
    log("  COSMIC variants also in dbSNP: %d of %d (%.1f%%)"
        % ((cos & dbsnp).sum(), cos.sum(), 100.0 * dbsnp[cos].mean()))
    log("  in COSMIC only:                %d" % (cos & ~dbsnp).sum())
log()

# 2.
log("2. Feature importance")
log()
imp_all = SimpleImputer(strategy="median")
Xf = imp_all.fit_transform(df[FULL].values)
gb = GradientBoostingClassifier(**PARAMS).fit(Xf, y)
imps = pd.Series(gb.feature_importances_, index=FULL).sort_values(ascending=False)
log("top ten features:")
for name, v in imps.head(10).items():
    log("  %-22s %.4f" % (name, v))
log()
rank = list(imps.index).index("in_cosmic") + 1
log("  in_cosmic is ranked %d of %d, importance %.4f"
    % (rank, len(FULL), imps["in_cosmic"]))
log()
log("other database features:")
for f in DATABASE:
    log("  %-22s %.4f   (rank %d)"
        % (f, imps[f], list(imps.index).index(f) + 1))
log()

# 3.
log("3. With and without, cross-validation by variant")
log()
p_full = oof(FULL)
p_noco = oof(NOCOSMIC)
p_ncell = oof(NOCELL)
p_ncell_noco = oof(NOCELL_NOCOSMIC)

log("  %-34s %9s %9s" % ("", "ROC-AUC", "PR-AUC"))
for name, p in [("35 features, with COSMIC", p_full),
                ("34 features, without COSMIC", p_noco),
                ("31 features, with COSMIC, no cross-cell", p_ncell),
                ("30 features, without COSMIC, no cross-cell", p_ncell_noco)]:
    log("  %-34s %9.4f %9.4f"
        % (name, roc_auc_score(y, p), average_precision_score(y, p)))
log()
log("  cost of a missing COSMIC column:")
log("    full input:       %+.4f ROC-AUC"
    % (roc_auc_score(y, p_noco) - roc_auc_score(y, p_full)))
log("    single file:      %+.4f ROC-AUC"
    % (roc_auc_score(y, p_ncell_noco) - roc_auc_score(y, p_ncell)))
log()

# 4.
log("4. By cell line")
log()
log("  %-12s %8s %9s %9s %10s" % ("line", "rows", "with COSMIC", "without", "difference"))
for ln in sorted(set(line)):
    m = line == ln
    if m.sum() < 50 or not (0 < y[m].mean() < 1):
        continue
    a = roc_auc_score(y[m], p_full[m])
    b = roc_auc_score(y[m], p_noco[m])
    log("  %-12s %8d %9.4f %9.4f %+10.4f" % (ln, m.sum(), a, b, b - a))
norm = line != CANCER
log("  %-12s %8d %9.4f %9.4f %+10.4f"
    % ("normal lines", norm.sum(), roc_auc_score(y[norm], p_full[norm]),
       roc_auc_score(y[norm], p_noco[norm]),
       roc_auc_score(y[norm], p_noco[norm]) - roc_auc_score(y[norm], p_full[norm])))
log()

# 5.
log("5. Tumour line only")
log()
sub = df[line == CANCER].reset_index(drop=True)
ys = (sub["label"] == "RNA_only").astype(int).values
vs = (sub["chr"].astype(str) + ":" + sub["pos"].astype(str) + ":"
      + sub["ref"].astype(str) + ">" + sub["alt"].astype(str)).values
log("  rows %d, RNA_only share %.3f, in COSMIC %d"
    % (len(sub), ys.mean(), (pd.to_numeric(sub["in_cosmic"], errors="coerce")
                             .fillna(0).values > 0).sum()))
log()
log("  %-26s %9s %9s" % ("", "ROC-AUC", "PR-AUC"))
for name, feats in [("with COSMIC", FULL), ("without COSMIC", NOCOSMIC)]:
    X = sub[feats].values
    p = np.zeros(len(sub))
    for itr, ite in GroupKFold(n_splits=5).split(X, ys, groups=vs):
        im = SimpleImputer(strategy="median")
        Xtr, Xte = im.fit_transform(X[itr]), im.transform(X[ite])
        p[ite] = GradientBoostingClassifier(**PARAMS).fit(
            Xtr, ys[itr]).predict_proba(Xte)[:, 1]
    log("  %-26s %9.4f %9.4f"
        % (name, roc_auc_score(ys, p), average_precision_score(ys, p)))
log()


path = os.path.join(outdir, "cosmic_contribution.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
