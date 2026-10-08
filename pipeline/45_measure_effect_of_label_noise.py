#!/usr/bin/env python3
"""Does label noise change the numbers, and does it change the calibration?

Usage: 45_measure_effect_of_label_noise.py <dataset_sig.csv> <corrected.csv>
"""
import sys

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (average_precision_score, brier_score_loss, f1_score,
                             recall_score, roc_auc_score)
from sklearn.model_selection import GroupKFold

orig_csv, corr_csv = sys.argv[1], sys.argv[2]
df = pd.read_csv(orig_csv).dropna(subset=["label"]).reset_index(drop=True)
corr = pd.read_csv(corr_csv, usecols=["cell", "chr", "pos", "ref", "alt", "label"])
corr = corr.rename(columns={"label": "label_corrected"})
df = df.merge(corr, on=["cell", "chr", "pos", "ref", "alt"], how="left")
assert df.label_corrected.notna().all()

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
CORE = BASELINE + READLEVEL + SEQCTX
FEATS = [("data-only (30)", CORE + CROSSCELL), ("database-informed (35)", CORE + CROSSCELL + DATABASE)]
for c in CORE + CROSSCELL + DATABASE:
    df[c] = pd.to_numeric(df[c], errors="coerce")

PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)
LINES = ["molm13", "hg001", "hg002"]
y_orig = (df.label == "RNA_only").astype(int).values
y_fix = (df.label_corrected == "RNA_only").astype(int).values
wrong = y_orig != y_fix
vid = (df.chr.astype(str) + ":" + df.pos.astype(str) + ":" + df.ref + ">" + df.alt).values
lines = df.group.values
everyone, kept = np.ones(len(df), bool), ~wrong


def fit_predict(X, y, tr, te, calibrated=False):
    imp = SimpleImputer(strategy="median")
    gb = GradientBoostingClassifier(**PARAMS)
    est = CalibratedClassifierCV(gb, method="isotonic", cv=3) if calibrated else gb
    est.fit(imp.fit_transform(X[tr]), y[tr])
    return est.predict_proba(imp.transform(X[te]))[:, 1]


def oof(X, y, mask, calibrated=False):
    idx = np.where(mask)[0]
    p = np.full(len(df), np.nan)
    for a, b in GroupKFold(n_splits=5).split(X[idx], y[idx], groups=vid[idx]):
        p[idx[b]] = fit_predict(X, y, idx[a], idx[b], calibrated)
    return p


def logo(X, y, mask, held):
    tr = np.where(mask & (lines != held))[0]
    te = np.where(mask & (lines == held))[0]
    p = np.full(len(df), np.nan)
    p[te] = fit_predict(X, y, tr, te)
    return p


def met(y, p, mask):
    y, p = y[mask], p[mask]
    pred = (p >= 0.5).astype(int)
    return "%.3f  %.3f  %.3f  %.3f" % (roc_auc_score(y, p), average_precision_score(y, p),
                                       recall_score(y, pred), f1_score(y, pred))


def ece(y, p, bins=10):
    e, edges = 0.0, np.linspace(0, 1, bins + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (p >= lo) & ((p < hi) if hi < 1 else (p <= hi))
        if s.any():
            e += s.mean() * abs(p[s].mean() - y[s].mean())
    return e


for l in LINES:
    s = lines == l
    print("%-7s rows %5d, RNA_only %5d, corrected %4d" % (l, s.sum(), y_orig[s].sum(), (wrong & s).sum()))
print()

for name, feats in FEATS:
    X = df[feats].values
    p_a, p_c, p_d = oof(X, y_orig, everyone), oof(X, y_fix, everyone), oof(X, y_orig, kept)
    print("%s, cross-validation by variant: ROC-AUC  PR-AUC  recall  F1" % name)
    print("A ", met(y_orig, p_a, everyone))
    print("B ", met(y_fix, p_a, everyone))
    print("C ", met(y_fix, p_c, everyone))
    print("D ", met(y_orig, p_d, kept))
    print("   mean P(RNA_only) under A: wrong rows %.3f, true RNA_only %.3f, DNA_origin %.3f"
          % (p_a[wrong].mean(), p_a[(y_orig == 1) & ~wrong].mean(), p_a[y_orig == 0].mean()))
    print("leave one line out, ROC-AUC A / B / C / D")
    for held in LINES:
        te = lines == held
        pa, pc, pd_ = logo(X, y_orig, everyone, held), logo(X, y_fix, everyone, held), logo(X, y_orig, kept, held)
        print("   %-7s %.3f / %.3f / %.3f / %.3f" % (
            held, roc_auc_score(y_orig[te], pa[te]), roc_auc_score(y_fix[te], pa[te]),
            roc_auc_score(y_fix[te], pc[te]), roc_auc_score(y_orig[te & kept], pd_[te & kept])))
    q_a, q_c = oof(X, y_orig, everyone, True), oof(X, y_fix, everyone, True)
    print("calibration (isotonic, cv=3): Brier / ECE")
    print("   trained original, judged original:   %.4f / %.4f" % (brier_score_loss(y_orig, q_a), ece(y_orig, q_a)))
    print("   trained original, judged corrected:  %.4f / %.4f" % (brier_score_loss(y_fix, q_a), ece(y_fix, q_a)))
    print("   trained corrected, judged corrected: %.4f / %.4f" % (brier_score_loss(y_fix, q_c), ece(y_fix, q_c)))
    print()
