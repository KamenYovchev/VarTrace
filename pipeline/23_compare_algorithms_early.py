#!/usr/bin/env python3
"""Early comparison of algorithms for the SNV-origin classifier.

Usage: 23_compare_algorithms_early.py <dataset.csv>
"""
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.metrics import (accuracy_score, precision_score, recall_score,
                             f1_score, roc_auc_score)

df = pd.read_csv(sys.argv[1])
df = df.dropna(subset=["label"])
y = (df["label"] == "RNA_only").astype(int).values
variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
              + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + \
         (["is_transition"] if "is_transition" in df.columns else [])
CROSSCELL = [c for c in ["n_cells_variant", "frac_cells", "vaf_mean_cells",
                         "vaf_sd_cells"] if c in df.columns]
DATABASE = [c for c in ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal",
                        "in_cosmic"] if c in df.columns]
SIGNATURES = [c for c in df.columns if c.startswith("sig_")]
for c in BASELINE + READLEVEL + SEQCTX + CROSSCELL + DATABASE + SIGNATURES:
    df[c] = pd.to_numeric(df[c], errors="coerce")

idx = np.arange(len(df))
# Row split (continuity with earlier results)
itr_r, ite_r = train_test_split(idx, test_size=0.25, random_state=42, stratify=y)
# Variant split (honest: same variant never in both sides)
gss = GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=42)
itr_v, ite_v = next(gss.split(idx, y, groups=variant_id))


def models():
    return {
        "LogReg": make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                LogisticRegression(max_iter=1000, class_weight="balanced")),
        "RF": make_pipeline(SimpleImputer(strategy="median"),
                            RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                                   random_state=42, n_jobs=-1)),
        "GB": make_pipeline(SimpleImputer(strategy="median"),
                            GradientBoostingClassifier(random_state=42)),
    }


def auc_by_model(feats, itr, ite):
    X = df[feats].values
    out = {}
    for name, m in models().items():
        m.fit(X[itr], y[itr])
        out[name] = roc_auc_score(y[ite], m.predict_proba(X[ite])[:, 1])
    return out


print("features: baseline=%d read-level=%d seq-context=%d cross-cell=%d database=%d signatures=%d"
      % (len(BASELINE), len(READLEVEL), len(SEQCTX), len(CROSSCELL), len(DATABASE),
         len(SIGNATURES)))
print("rows=%d  unique variants=%d  RNA_only fraction=%.3f"
      % (len(df), len(set(variant_id)), y.mean()))

experiments = [
    ("baseline (6)", BASELINE),
    ("+ read-level", BASELINE + READLEVEL),
    ("+ seq-context", BASELINE + READLEVEL + SEQCTX),
    ("+ cross-cell", BASELINE + READLEVEL + SEQCTX + CROSSCELL),
    ("+ signatures", BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES),
    ("+ database", BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES + DATABASE),
]
for split_name, itr, ite in [("ROW split", itr_r, ite_r),
                             ("VARIANT split (honest)", itr_v, ite_v)]:
    print("\nROC-AUC, %s  (train=%d test=%d)" % (split_name, len(itr), len(ite)))
    print("%-22s %7s %7s %7s" % ("experiment", "LogReg", "RF", "GB"))
    for name, feats in experiments:
        r = auc_by_model(feats, itr, ite)
        print("%-22s %7.3f %7.3f %7.3f" % (name, r["LogReg"], r["RF"], r["GB"]))

# Full metrics + importance: all features, honest variant split
feats = BASELINE + READLEVEL + SEQCTX + CROSSCELL + SIGNATURES + DATABASE
X = df[feats].values
rf = make_pipeline(SimpleImputer(strategy="median"),
                   RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                          random_state=42, n_jobs=-1))
rf.fit(X[itr_v], y[itr_v])
pred, proba = rf.predict(X[ite_v]), rf.predict_proba(X[ite_v])[:, 1]
print("\nfull metrics: RF, all features, variant split")
print("acc=%.3f  prec=%.3f  rec=%.3f  f1=%.3f  rocauc=%.3f" % (
    accuracy_score(y[ite_v], pred), precision_score(y[ite_v], pred, zero_division=0),
    recall_score(y[ite_v], pred, zero_division=0), f1_score(y[ite_v], pred, zero_division=0),
    roc_auc_score(y[ite_v], proba)))
imp = sorted(zip(feats, rf.named_steps["randomforestclassifier"].feature_importances_),
             key=lambda t: -t[1])
print("\ntop 15 feature importance (RF):")
for f, v in imp[:15]:
    print("  %-16s %.3f" % (f, v))
