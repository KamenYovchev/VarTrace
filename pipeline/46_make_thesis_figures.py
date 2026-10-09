#!/usr/bin/env python3
"""The figures of the results chapter.

ROC, calibration and SHAP are computed here. The feature-group progression is
plotted from the numbers step 26 printed.

Usage: 46_make_thesis_figures.py <dataset_sig.csv> <outdir>
"""
# Plot labels are in Bulgarian because the figures go into the thesis.
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GroupKFold

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"]).reset_index(drop=True)
y = (df.label == "RNA_only").astype(int).values
vid = (df.chr.astype(str) + ":" + df.pos.astype(str) + ":" + df.ref + ">" + df.alt).values

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = [c for c in ["fs", "sor", "mq", "mqranksum", "readposranksum",
                         "baseqranksum", "qd"] if c in df.columns]
SEQCTX = [c for c in df.columns if c.startswith("sub_")] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
DATA_ONLY = BASELINE + READLEVEL + SEQCTX + CROSSCELL
WITH_DB = DATA_ONLY + DATABASE
for c in WITH_DB:
    df[c] = pd.to_numeric(df[c], errors="coerce")
PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)


def oof(feats, calibrated=False):
    X = df[feats].values
    p = np.zeros(len(df))
    for tr, te in GroupKFold(n_splits=5).split(X, y, groups=vid):
        imp = SimpleImputer(strategy="median")
        gb = GradientBoostingClassifier(**PARAMS)
        est = CalibratedClassifierCV(gb, method="isotonic", cv=3) if calibrated else gb
        est.fit(imp.fit_transform(X[tr]), y[tr])
        p[te] = est.predict_proba(imp.transform(X[te]))[:, 1]
    return p


# Feature-group progression, mean and sd over the five folds (step 26)
steps = ["базови", "+ прочити", "+ контекст", "+ между клетките", "+ сигнатури", "+ бази данни"]
means = [0.804, 0.831, 0.838, 0.910, 0.915, 0.968]
sds = [0.012, 0.009, 0.010, 0.012, 0.013, 0.007]
plt.figure(figsize=(7, 4.2))
xs = np.arange(len(steps))
plt.errorbar(xs, means, yerr=sds, marker="o", capsize=4)
plt.xticks(xs, steps, rotation=20)
plt.ylabel("ROC-AUC, кръстосана проверка по вариант")
plt.title("Принос на групите признаци")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(outdir, "progression.png"), dpi=200)
plt.close()

p30, p35 = oof(DATA_ONLY), oof(WITH_DB)
plt.figure(figsize=(5.2, 5))
for name, p in [("само от данните", p30), ("с бази данни", p35)]:
    fpr, tpr, _ = roc_curve(y, p)
    plt.plot(fpr, tpr, label="%s (ROC-AUC %.3f)" % (name, roc_auc_score(y, p)))
plt.plot([0, 1], [0, 1], "k--", lw=0.8)
plt.xlabel("дял на неверно положителните")
plt.ylabel("дял на верно положителните")
plt.title("Работна характеристика на двата модела")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "roc.png"), dpi=200)
plt.close()
print("ROC-AUC 30: %.4f  35: %.4f" % (roc_auc_score(y, p30), roc_auc_score(y, p35)))

q35 = oof(WITH_DB, calibrated=True)
frac, mean = calibration_curve(y, q35, n_bins=10)
plt.figure(figsize=(5.2, 5))
plt.plot(mean, frac, marker="o", label="модел с бази данни")
plt.plot([0, 1], [0, 1], "k--", lw=0.8, label="съвършено калибриране")
plt.xlabel("предсказана вероятност за възникнал в РНК")
plt.ylabel("действителен дял на възникналите в РНК")
plt.title("Калибриране на оценката")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "calibration.png"), dpi=200)
plt.close()
print("calibrated, ROC-AUC %.4f" % roc_auc_score(y, q35))
for name, p in [("uncalibrated", p35), ("calibrated", q35)]:
    hit = p >= 0.867
    print("%s at 0.867: precision %.4f  recall %.4f"
          % (name, y[hit].mean(), hit[y == 1].mean()))
for a, b in zip(mean, frac):
    print("  predicted %.3f  observed %.3f" % (a, b))

X = SimpleImputer(strategy="median").fit_transform(df[WITH_DB].values)
m = GradientBoostingClassifier(**PARAMS).fit(X, y)
sv = shap.TreeExplainer(m).shap_values(X)
plt.figure()
shap.summary_plot(sv, features=X, feature_names=WITH_DB, show=False, max_display=15)
plt.xlabel("стойност по SHAP, влияние върху решението")
cb = plt.gcf().axes[-1]
cb.set_ylabel("стойност на признака")
cb.set_yticklabels(["ниска", "висока"])
plt.tight_layout()
plt.savefig(os.path.join(outdir, "shap_summary.png"), dpi=200, bbox_inches="tight")
plt.close()
order = np.argsort(-np.abs(sv).mean(axis=0))
print("SHAP, top features:", [WITH_DB[i] for i in order[:6]])

# Importance table: the built-in importance of the calibrated model, averaged
# over its three calibration parts. Same model as the one saved in step 27.
imp = SimpleImputer(strategy="median")
cal = CalibratedClassifierCV(GradientBoostingClassifier(**PARAMS), method="isotonic",
                             cv=3).fit(imp.fit_transform(df[WITH_DB].values), y)
fi = np.mean([c.estimator.feature_importances_ for c in cal.calibrated_classifiers_], axis=0)
order = np.argsort(-fi)
print("importance:")
for i in order[:6]:
    print("  %-18s %.3f" % (WITH_DB[i], fi[i]))
print("  other %d together %.3f" % (len(order) - 6, fi[order[6:]].sum()))
