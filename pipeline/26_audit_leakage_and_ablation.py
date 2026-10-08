#!/usr/bin/env python3
"""Three validation checks that close out the analysis.

Usage: 26_audit_leakage_and_ablation.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import math
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.model_selection import GroupKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (roc_auc_score, average_precision_score, precision_score,
                             recall_score, f1_score, accuracy_score)

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"])
y = (df["label"] == "RNA_only").astype(int).values
variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
              + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
groups_line = df["group"].values
# The exact aggregation key used by 21_build_feature_table.py for the cross-cell
# features
agg_key = (df["group"].astype(str) + "|" + pd.Series(variant_id, index=df.index)).values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
SIGNATURES = [c for c in df.columns if c.startswith("sig_")]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
ALL = BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES + DATABASE
DATA_ONLY = BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES
for c in ALL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

LINES = ["molm13", "hg001", "hg002"]
gkf = GroupKFold(n_splits=5)
out = []


def log(s=""):
    print(s)
    out.append(s)


def model():
    return GradientBoostingClassifier(random_state=42)


# 1.
log("1. leakage audit OF the cross-cell features")
log("rows=%d  unique variants=%d  unique aggregation keys=%d  cell lines=%d"
    % (len(df), len(set(variant_id)), len(set(agg_key)), len(set(groups_line))))
log()

# 1a.
k2v = defaultdict(set)
for k, v in zip(agg_key, variant_id):
    k2v[k].add(v)
bad_v = {k: vs for k, vs in k2v.items() if len(vs) > 1}
log("1a. aggregation key -> variant_id must be many-to-one")
log("    aggregation keys mapping to >1 variant: %d  %s"
    % (len(bad_v), "OK" if not bad_v else "*** FAIL ***"))

# 1b.
k2l = defaultdict(set)
for k, g in zip(agg_key, groups_line):
    k2l[k].add(g)
bad_l = {k: gs for k, gs in k2l.items() if len(gs) > 1}
log("1b. aggregation key -> cell line must be many-to-one")
log("    aggregation keys spanning >1 cell line: %d  %s"
    % (len(bad_l), "OK" if not bad_l else "*** FAIL ***"))
log()

# 1c.
log("1c. per fold, aggregation keys shared between train and test (must be 0)")
tot_shared = 0
for i, (itr, ite) in enumerate(gkf.split(df[ALL].values, y, groups=variant_id), 1):
    shared = set(agg_key[itr]) & set(agg_key[ite])
    v_shared = set(variant_id[itr]) & set(variant_id[ite])
    tot_shared += len(shared)
    log("    fold %d: n_train=%-6d n_test=%-6d shared agg keys=%d  shared variants=%d"
        % (i, len(itr), len(ite), len(shared), len(v_shared)))
log("    total shared across folds: %d  %s"
    % (tot_shared, "OK" if tot_shared == 0 else "*** FAIL ***"))
log()

# 1d.
log("1d. leave-one-cell-line-out, aggregation keys shared with the rest (must be 0)")
for held in LINES:
    te = groups_line == held
    shared = set(agg_key[te]) & set(agg_key[~te])
    log("    held out %-7s n_test=%-6d shared agg keys=%d" % (held, te.sum(), len(shared)))
log()

# 1e.
log("1e. independent recomputation of the cross-cell features")
cells_per_line = df.groupby("group")["cell"].nunique().to_dict()
vafs = defaultdict(list)
for k, v in zip(agg_key, pd.to_numeric(df["vaf"], errors="coerce").values):
    vafs[k].append(v)
rec_n, rec_frac, rec_mean, rec_sd = [], [], [], []
for k, g in zip(agg_key, groups_line):
    vs_all = vafs[k]
    vs = [v for v in vs_all if not (v is None or (isinstance(v, float) and math.isnan(v)))]
    rec_n.append(len(vs_all))
    rec_frac.append(round(len(vs_all) / cells_per_line[g], 4) if cells_per_line[g] else np.nan)
    if vs:
        m = sum(vs) / len(vs)
        s = math.sqrt(sum((v - m) ** 2 for v in vs) / len(vs)) if len(vs) > 1 else 0.0
        rec_mean.append(round(m, 4))
        rec_sd.append(round(s, 4))
    else:
        rec_mean.append(np.nan)
        rec_sd.append(np.nan)
for name, rec in [("n_cells_variant", rec_n), ("frac_cells", rec_frac),
                  ("vaf_mean_cells", rec_mean), ("vaf_sd_cells", rec_sd)]:
    stored = pd.to_numeric(df[name], errors="coerce").values
    r = np.asarray(rec, dtype=float)
    both_nan = np.isnan(stored) & np.isnan(r)
    close = np.isclose(stored, r, rtol=0, atol=1e-3, equal_nan=False) | both_nan
    log("    %-16s identical in %d/%d rows  %s"
        % (name, close.sum(), len(df), "OK" if close.all() else "*** MISMATCH ***"))
log()


# 2.
def metrics(ytrue, proba, thr=0.5):
    pred = (proba >= thr).astype(int)
    return dict(roc_auc=roc_auc_score(ytrue, proba),
                pr_auc=average_precision_score(ytrue, proba),
                precision=precision_score(ytrue, pred, zero_division=0),
                recall=recall_score(ytrue, pred, zero_division=0),
                f1=f1_score(ytrue, pred, zero_division=0),
                accuracy=accuracy_score(ytrue, pred),
                prevalence=float(np.mean(ytrue)))


def oof_predict(feats):
    """Out-of-fold probabilities over the 5 variant-grouped folds."""
    X = df[feats].values
    proba = np.zeros(len(df))
    for itr, ite in gkf.split(X, y, groups=variant_id):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        proba[ite] = model().fit(Xtr, y[itr]).predict_proba(Xte)[:, 1]
    return proba


def logo_predict(feats, held):
    X = df[feats].values
    te = groups_line == held
    imp = SimpleImputer(strategy="median")
    Xtr, Xte = imp.fit_transform(X[~te]), imp.transform(X[te])
    return y[te], model().fit(Xtr, y[~te]).predict_proba(Xte)[:, 1]


HDR = ("%-28s %7s %7s %9s %7s %7s %8s %10s"
       % ("", "ROC-AUC", "PR-AUC", "precision", "recall", "F1", "accuracy", "prevalence"))


def row(label, m):
    return ("%-28s %7.3f %7.3f %9.3f %7.3f %7.3f %8.3f %10.3f"
            % (label, m["roc_auc"], m["pr_auc"], m["precision"],
               m["recall"], m["f1"], m["accuracy"], m["prevalence"]))


log("2. full metrics (GradientBoosting, threshold 0.5, positive class = RNA_only)")
log()
log("2a. 5-fold cross-validation split BY variant (pooled out-of-fold predictions)")
log(HDR)
oof = {}
for feats, fname in [(DATA_ONLY, "data-only"), (ALL, "database-informed")]:
    oof[fname] = oof_predict(feats)
    log(row(fname, metrics(y, oof[fname])))
log()
log("2b. leave-one-cell-line-out")
log(HDR)
for feats, fname in [(DATA_ONLY, "data-only"), (ALL, "database-informed")]:
    for held in LINES:
        yt, pr = logo_predict(feats, held)
        log(row("%s / %s" % (fname, held), metrics(yt, pr)))
log()


# 3.
log("3. ablation BY feature group")
log()

GROUPS = [("baseline", BASELINE), ("read-level", READLEVEL), ("seq-context", SEQCTX),
          ("cross-cell", CROSSCELL), ("signatures", SIGNATURES), ("database", DATABASE)]


def cv_auc_sd(feats):
    X = df[feats].values
    aucs = []
    for itr, ite in gkf.split(X, y, groups=variant_id):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        m = model().fit(Xtr, y[itr])
        aucs.append(roc_auc_score(y[ite], m.predict_proba(Xte)[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


log("3a. cumulative (CV by variant, ROC-AUC mean +- sd over folds)")
cum = []
for name, g in GROUPS:
    cum += g
    mu, sd = cv_auc_sd(cum)
    log("    + %-12s (%3d feats)  %.3f +- %.3f" % (name, len(cum), mu, sd))
log()

log("3b. leave-one-group-out from the full set (CV by variant)")
mu_all, sd_all = cv_auc_sd(ALL)
log("    all           (%3d feats)  %.3f +- %.3f" % (len(ALL), mu_all, sd_all))
for name, g in GROUPS:
    feats = [c for c in ALL if c not in set(g)]
    mu, sd = cv_auc_sd(feats)
    log("    all - %-9s (%3d feats)  %.3f +- %.3f   delta=%+.3f"
        % (name, len(feats), mu, sd, mu - mu_all))
log()

log("3c. do the mutational signatures still earn their place?")
NO_SIG = [c for c in ALL if not c.startswith("sig_")]
mu_ns, sd_ns = cv_auc_sd(NO_SIG)
log("    CV by variant:  all %.3f +- %.3f   vs  all-without-signatures %.3f +- %.3f"
    % (mu_all, sd_all, mu_ns, sd_ns))
log("                    delta = %+.4f, pooled sd = %.4f -> %s"
    % (mu_ns - mu_all, (sd_all + sd_ns) / 2,
       "WITHIN NOISE, signatures can be dropped"
       if abs(mu_ns - mu_all) <= (sd_all + sd_ns) / 2 else
       "OUTSIDE NOISE, signatures contribute"))
log()
log("    leave-one-cell-line-out (ROC-AUC):")
log("    %-10s %10s %10s %8s" % ("held out", "all", "no-signat.", "delta"))
for held in LINES:
    _, p_all = logo_predict(ALL, held)
    yt, p_ns = logo_predict(NO_SIG, held)
    a_all = roc_auc_score(yt, p_all)
    a_ns = roc_auc_score(yt, p_ns)
    log("    %-10s %10.3f %10.3f %+8.3f" % (held, a_all, a_ns, a_ns - a_all))
log()
log("    the same, for the data-only variant:")
NO_SIG_DO = [c for c in DATA_ONLY if not c.startswith("sig_")]
log("    %-10s %10s %10s %8s" % ("held out", "data-only", "no-signat.", "delta"))
for held in LINES:
    _, p_do = logo_predict(DATA_ONLY, held)
    yt, p_ns = logo_predict(NO_SIG_DO, held)
    a_do = roc_auc_score(yt, p_do)
    a_ns = roc_auc_score(yt, p_ns)
    log("    %-10s %10.3f %10.3f %+8.3f" % (held, a_do, a_ns, a_ns - a_do))
log()

path = os.path.join(outdir, "validation_checks.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
