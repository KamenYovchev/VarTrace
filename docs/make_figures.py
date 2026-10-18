#!/usr/bin/env python3
"""Precision-recall and calibration figures for the README.

Same features, settings and variant-grouped folds as
pipeline/27_train_four_shipped_models.py, on the table with corrected labels
(pipeline step 44), which is what the shipped models are trained on. Inside
each fold the model is calibrated the way the shipped models are, so the
calibration curve shows how the shipped models behave on variants they did
not see.

Usage:
    make_figures.py scores <dataset_corrected.csv> <scores.csv>
    make_figures.py plot <scores.csv> <outdir>
"""
import sys
import os
import numpy as np
import pandas as pd

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
READLEVEL = ["fs", "sor", "mq", "mqranksum", "readposranksum", "baseqranksum", "qd"]
SUBS = ["A_C", "A_G", "A_T", "C_A", "C_G", "C_T", "G_A", "G_C", "G_T", "T_A", "T_C", "T_G"]
SEQCTX = ["sub_" + s for s in SUBS] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
CORE = BASELINE + READLEVEL + SEQCTX
MODELS = {"with_databases": CORE + CROSSCELL + DATABASE,
          "data_only": CORE + CROSSCELL}
PARAMS = dict(learning_rate=0.05, max_depth=3, n_estimators=300, random_state=42)


def scores(csv_in, csv_out):
    from sklearn.model_selection import GroupKFold
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.impute import SimpleImputer

    df = pd.read_csv(csv_in).dropna(subset=["label"])
    y = (df["label"] == "RNA_only").astype(int).values
    variant_id = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
                  + df["ref"].astype(str) + ">" + df["alt"].astype(str)).values
    out = pd.DataFrame({"label": y})
    for name, feats in MODELS.items():
        X = df[feats].apply(pd.to_numeric, errors="coerce").values
        raw, cal = np.zeros(len(df)), np.zeros(len(df))
        for itr, ite in GroupKFold(n_splits=5).split(X, y, groups=variant_id):
            imp = SimpleImputer(strategy="median")
            Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
            raw[ite] = GradientBoostingClassifier(**PARAMS).fit(
                Xtr, y[itr]).predict_proba(Xte)[:, 1]
            cal[ite] = CalibratedClassifierCV(
                GradientBoostingClassifier(**PARAMS), method="isotonic",
                cv=3).fit(Xtr, y[itr]).predict_proba(Xte)[:, 1]
        out[name + "_raw"], out[name + "_cal"] = raw, cal
    out.to_csv(csv_out, index=False)


def ece(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, bins - 1)
    total = 0.0
    for b in range(bins):
        m = idx == b
        if m.any():
            total += m.mean() * abs(p[m].mean() - y[m].mean())
    return total


def plot(csv_in, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import (precision_recall_curve, average_precision_score,
                                 roc_auc_score, brier_score_loss)

    s = pd.read_csv(csv_in)
    y = s["label"].values
    os.makedirs(outdir, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    labels = {"with_databases": "with databases, 35 features",
              "data_only": "RNA data only, 30 features"}
    styles = {"with_databases": "-", "data_only": "--"}

    for name in MODELS:
        p, c = s[name + "_raw"].values, s[name + "_cal"].values
        print("%-15s ROC-AUC %.3f  PR-AUC %.3f  | calibrated: ECE %.3f  Brier %.3f"
              % (name, roc_auc_score(y, p), average_precision_score(y, p),
                 ece(y, c), brier_score_loss(y, c)))

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for name in MODELS:
        p = s[name + "_raw"].values
        prec, rec, _ = precision_recall_curve(y, p)
        ax.plot(rec, prec, styles[name], color="black", lw=1.6,
                label="%s, PR-AUC %.3f" % (labels[name], average_precision_score(y, p)))
    ax.axhline(y.mean(), color="grey", lw=1, ls=":",
               label="share of RNA-level variants, %.3f" % y.mean())
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "precision_recall.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    ax.plot([0, 1], [0, 1], color="grey", lw=1, ls=":", label="perfect calibration")
    edges = np.linspace(0, 1, 11)
    for name in MODELS:
        c = s[name + "_cal"].values
        idx = np.clip(np.digitize(c, edges[1:-1]), 0, 9)
        xs = [c[idx == b].mean() for b in range(10) if (idx == b).sum() >= 20]
        ys = [y[idx == b].mean() for b in range(10) if (idx == b).sum() >= 20]
        ax.plot(xs, ys, styles[name], marker="o", ms=4, color="black", lw=1.6,
                label="%s, ECE %.3f" % (labels[name], ece(y, c)))
    ax.set_xlabel("predicted probability of being absent from the genome")
    ax.set_ylabel("observed share absent from the genome")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "calibration.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    if sys.argv[1] == "scores":
        scores(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "plot":
        plot(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
