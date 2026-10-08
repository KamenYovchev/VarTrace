#!/usr/bin/env python3
"""XGBoost, held-out cell lines and a neural net for the SNV-origin classifier.

Usage: 25_compare_algorithms_and_lines.py <dataset_sig.csv> <outdir>
"""
import sys
import os
import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, ParameterGrid
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import roc_auc_score

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"])
y = (df["label"] == "RNA_only").astype(int).values
variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
              + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
groups_line = df["group"].values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
SIGNATURES = [c for c in df.columns if c.startswith("sig_")]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
# The 86 signature features were dropped after the ablation in
# 26_audit_leakage_and_ablation.py showed they add nothing measurable (0.969 without vs
# 0.968 with, fold sd 0.007; identical to within 0.001 under leave-one-cell-line-out).
ALL = BASELINE + READLEVEL + SEQCTX + CROSSCELL + DATABASE
DATA_ONLY = BASELINE + READLEVEL + SEQCTX + CROSSCELL
for c in ALL:
    df[c] = pd.to_numeric(df[c], errors="coerce")

gkf = GroupKFold(n_splits=5)
lines_out = []


def log(s):
    print(s)
    lines_out.append(s)


def cv_auc(model_fn, feats):
    X = df[feats].values
    aucs = []
    for itr, ite in gkf.split(X, y, groups=variant_id):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        m = model_fn().fit(Xtr, y[itr])
        aucs.append(roc_auc_score(y[ite], m.predict_proba(Xte)[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


# 1.
GRID = ParameterGrid({"n_estimators": [300, 600],
                      "max_depth": [3, 5],
                      "learning_rate": [0.05, 0.1]})

log("1. algorithm comparison (5-fold CV by variant, same grid for both)")


def search(make, label):
    best, best_mu, best_sd = None, -1.0, 0.0
    for params in GRID:
        mu, sd = cv_auc(lambda p=params: make(**p), ALL)
        log("    %-18s %-52s %.3f +- %.3f" % (label, str(params), mu, sd))
        if mu > best_mu:
            best, best_mu, best_sd = params, mu, sd
    return best, best_mu, best_sd


gb_best, gb_mu, gb_sd = search(
    lambda **kw: GradientBoostingClassifier(random_state=42, **kw),
    "GradientBoosting")
log("GradientBoosting best %-40s %.3f +- %.3f" % (str(gb_best), gb_mu, gb_sd))

try:
    from xgboost import XGBClassifier

    def xgb(**kw):
        return XGBClassifier(eval_metric="logloss", random_state=42,
                             n_jobs=-1, **kw)

    xgb_best, xgb_mu, xgb_sd = search(xgb, "XGBoost")
    log("XGBoost best          %-40s %.3f +- %.3f" % (str(xgb_best), xgb_mu, xgb_sd))
    have_xgb = True
except ImportError:
    log("xgboost not installed - GradientBoosting kept")
    xgb_best, xgb_mu, xgb_sd, have_xgb = None, -1.0, 0.0, False

pooled_sd = (gb_sd + xgb_sd) / 2 if have_xgb else gb_sd
log("\ndifference = %+.4f, pooled sd = %.4f" % (xgb_mu - gb_mu, pooled_sd))

if have_xgb and (xgb_mu - gb_mu) > pooled_sd:
    final_model_fn = lambda: xgb(**xgb_best)
    final_name = "XGBoost %s" % xgb_best
    best_mu = xgb_mu
    log("XGBoost exceeds by more than one pooled sd -> chosen")
else:
    final_model_fn = lambda: GradientBoostingClassifier(random_state=42, **gb_best)
    final_name = "GradientBoosting %s" % gb_best
    best_mu = gb_mu
    log("difference within noise -> the simpler model is kept")
log("chosen final algorithm: %s" % final_name)
mu_gb = gb_mu

# 2.
log("\n2. generalization: leave-one-cell-line-out")
for feats, fname in [(DATA_ONLY, "data-only"), (ALL, "all features")]:
    X = df[feats].values
    for held in ["molm13", "hg001", "hg002"]:
        tr = groups_line != held
        te = ~tr
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[tr]), imp.transform(X[te])
        m = final_model_fn().fit(Xtr, y[tr])
        auc = roc_auc_score(y[te], m.predict_proba(Xte)[:, 1])
        log("train on others, test on %-7s (%s, n=%d)  AUC=%.3f"
            % (held, fname, te.sum(), auc))

# 3.
log("\n3. neural net (mlp) vs boosting, same features (5-fold CV)")
mu, sd = cv_auc(lambda: make_pipeline(
    StandardScaler(),
    MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=500,
                  early_stopping=True, random_state=42)), ALL)
log("mlp (64,32)                     %.3f +- %.3f" % (mu, sd))
log("(boosting above: %.3f) -> neural net %s at this data size"
    % (best_mu if best_mu > mu_gb else mu_gb,
       "does NOT beat boosting" if mu <= max(best_mu, mu_gb) else "beats boosting"))

# 4.
X = df[ALL].values
imp = SimpleImputer(strategy="median")
Xi = imp.fit_transform(X)
cal = CalibratedClassifierCV(final_model_fn(), method="isotonic", cv=3).fit(Xi, y)
joblib.dump({"model": cal, "imputer": imp, "features": ALL,
             "algorithm": final_name,
             "positive_class": "RNA_only"},
            os.path.join(outdir, "origin_model.pkl"))
log("\nfinal calibrated model saved: %s (algorithm: %s, %d features)"
    % (os.path.join(outdir, "origin_model.pkl"), final_name, len(ALL)))

with open(os.path.join(outdir, "improvements_results.txt"), "w") as f:
    f.write("\n".join(lines_out) + "\n")
