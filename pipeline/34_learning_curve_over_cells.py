#!/usr/bin/env python3
"""Would more cells have helped?

Usage: 34_learning_curve_over_cells.py <dataset_sig.csv> <outdir>
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
from sklearn.metrics import roc_auc_score, average_precision_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csv_in, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

df = pd.read_csv(csv_in).dropna(subset=["label"]).reset_index(drop=True)
df["y"] = (df["label"] == "RNA_only").astype(int)
df["vid"] = (df["chr"].astype(str) + ":" + df["pos"].astype(str) + ":"
             + df["ref"].astype(str) + ">" + df["alt"].astype(str))

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
SIZES = [5, 10, 15, 20, 25, 30, 35, 40, 45]
REPEATS = 5
rng = np.random.default_rng(42)
out = []


def log(s=""):
    print(s)
    out.append(s)


cells_by_line = {ln: sorted(g["cell"].unique())
                 for ln, g in df.groupby("group")}
total_cells = sum(len(v) for v in cells_by_line.values())


def draw(n):
    """Draw n cells, proportionally across the three lines."""
    picked = []
    for ln, cells in cells_by_line.items():
        k = int(round(n * len(cells) / total_cells))
        k = max(1, min(k, len(cells)))
        picked += list(rng.choice(cells, size=k, replace=False))
    return picked


def recompute_crosscell(sub):
    """Exactly 21_build_feature_table.py's logic, but over this subsample only."""
    n_cells_in_line = sub.groupby("group")["cell"].nunique().to_dict()
    key = sub["group"].astype(str) + "|" + sub["vid"]
    vafs = defaultdict(list)
    for k, v in zip(key, sub["vaf"].values):
        vafs[k].append(v)
    n_c, frac, mean, sd = [], [], [], []
    for k, ln in zip(key, sub["group"].values):
        vs_all = vafs[k]
        vs = [v for v in vs_all if not (v is None or (isinstance(v, float) and math.isnan(v)))]
        n_c.append(len(vs_all))
        frac.append(round(len(vs_all) / n_cells_in_line[ln], 4))
        if vs:
            m = sum(vs) / len(vs)
            s = math.sqrt(sum((v - m) ** 2 for v in vs) / len(vs)) if len(vs) > 1 else 0.0
            mean.append(round(m, 4))
            sd.append(round(s, 4))
        else:
            mean.append(np.nan)
            sd.append(np.nan)
    sub = sub.copy()
    sub["n_cells_variant"] = n_c
    sub["frac_cells"] = frac
    sub["vaf_mean_cells"] = mean
    sub["vaf_sd_cells"] = sd
    return sub


def evaluate(sub):
    yv = sub["y"].values
    if len(sub) < 200 or yv.mean() in (0.0, 1.0):
        return None
    X = sub[FEATS].values
    p = np.zeros(len(sub))
    n_splits = min(5, len(set(sub["vid"])) // 20)
    if n_splits < 2:
        return None
    for itr, ite in GroupKFold(n_splits=n_splits).split(X, yv, groups=sub["vid"].values):
        imp = SimpleImputer(strategy="median")
        Xtr, Xte = imp.fit_transform(X[itr]), imp.transform(X[ite])
        p[ite] = GradientBoostingClassifier(**PARAMS).fit(
            Xtr, yv[itr]).predict_proba(Xte)[:, 1]
    return roc_auc_score(yv, p), average_precision_score(yv, p), len(sub)


log("Learning curve over cells")
log()
log("%7s %8s %9s %17s %17s" % ("cells", "draws", "rows", "ROC-AUC", "PR-AUC"))

res = {}
for n in SIZES:
    reps = 1 if n >= total_cells else REPEATS
    aucs, prs, sizes = [], [], []
    for _ in range(reps):
        cells = draw(n) if n < total_cells else sorted(df["cell"].unique())
        sub = recompute_crosscell(df[df["cell"].isin(cells)])
        r = evaluate(sub)
        if r:
            aucs.append(r[0])
            prs.append(r[1])
            sizes.append(r[2])
    if aucs:
        res[n] = (np.mean(aucs), np.std(aucs), np.mean(prs), np.std(prs))
        log("%7d %8d %9.0f   %.3f +- %.3f   %.3f +- %.3f"
            % (n, len(aucs), np.mean(sizes), np.mean(aucs), np.std(aucs),
               np.mean(prs), np.std(prs)))

log()
ks = sorted(res)
if len(ks) >= 3:
    last_gain = res[ks[-1]][0] - res[ks[-3]][0]
    log("gain from the last %d cells: %+.4f ROC-AUC" % (ks[-1] - ks[-3], last_gain))
    log()
    if abs(last_gain) < 0.005:
        log("The curve has reached a plateau.")
    else:
        log("The curve is still rising.")

# Figure
if len(ks) >= 3:
    plt.figure(figsize=(6, 4.2))
    for idx, lbl, col in [(0, "ROC-AUC", None), (2, "PR-AUC", None)]:
        mu = [res[k][idx] for k in ks]
        sd = [res[k][idx + 1] for k in ks]
        plt.errorbar(ks, mu, yerr=sd, marker="o", capsize=3, label=lbl)
    # labels in Bulgarian, the figure goes into the thesis
    plt.xlabel("брой клетки в обучението")
    plt.ylabel("стойност (кръстосана проверка по вариант)")
    plt.title("Крива на учене")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "learning_curve.png"), dpi=200)
    plt.close()

path = os.path.join(outdir, "learning_curve.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s" % path)
