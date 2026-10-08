#!/usr/bin/env python3
"""Score the variants the paired DNA could not resolve.

Usage: 36_score_indeterminate_variants.py <labeled_dir> <model.pkl> <outdir>
"""
import sys
import os
import glob
import csv
import math
import joblib
import numpy as np
import pandas as pd
from collections import defaultdict

labeled_dir, model_pkl, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)

READLEVEL = [("fs", "FS"), ("sor", "SOR"), ("mq", "MQ"), ("mqranksum", "MQRankSum"),
             ("readposranksum", "ReadPosRankSum"), ("baseqranksum", "BaseQRankSum"),
             ("qd", "QD")]
SUBS = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
        "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
SETNAME = {"1": "DNA_origin", "2": "RNA_only"}


def group_of(base):
    if "MOLM13" in base:
        return "molm13"
    if base.startswith("ResolveOMEv2.1"):
        return "hg001"
    if "NA12878" in base:
        return "hg001"
    if base.startswith("ResolveOMEv2.0-384"):
        return "hg002"
    return "other"


def to_int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return np.nan


def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


rows = []
group_cells = defaultdict(set)
var_vafs = defaultdict(list)

for f in sorted(glob.glob(os.path.join(labeled_dir, "*.labeled.tsv"))):
    base = os.path.basename(f).replace(".labeled.tsv", "")
    grp = group_of(base)
    group_cells[grp].add(base)
    with open(f) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        i_chr, i_pos = idx["Otherinfo1"], idx["Otherinfo2"]
        i_ref, i_alt = idx["Otherinfo4"], idx["Otherinfo5"]
        i_qual, i_info = idx["Otherinfo6"], idx["Otherinfo8"]
        i_fmt, i_smp = idx["Otherinfo9"], idx["Otherinfo10"]
        i_set = idx["set"]
        i_dbsnp = idx.get("avsnp150")
        i_afs = [i for i, n in enumerate(header) if n == "AF"]
        i_cosmic = idx.get("cosmic104")
        i_redi = idx.get("rediportal")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            setv = c[i_set]
            ref, alt = c[i_ref], c[i_alt]
            d = dict(zip(c[i_fmt].split(":"), c[i_smp].split(":")))
            ref_ct = alt_ct = np.nan
            ad = d.get("AD", ".")
            if "," in ad:
                pp = ad.split(",")
                ref_ct, alt_ct = to_int(pp[0]), to_int(pp[1])
            dp = to_int(d.get("DP", "."))
            gq = to_int(d.get("GQ", "."))
            vaf = d.get("VAF")
            if vaf in (None, ".", ""):
                if not (np.isnan(ref_ct) or np.isnan(alt_ct)) and (ref_ct + alt_ct) > 0:
                    vaf = round(alt_ct / (ref_ct + alt_ct), 4)
                else:
                    vaf = np.nan
            else:
                try:
                    vaf = float(vaf)
                except ValueError:
                    vaf = np.nan
            try:
                qual = float(c[i_qual])
            except ValueError:
                qual = np.nan
            info = parse_info(c[i_info])
            rl = []
            for _, key in READLEVEL:
                v = info.get(key, "")
                try:
                    rl.append(float(v))
                except ValueError:
                    rl.append(np.nan)
            key = ref + ">" + alt
            sc = [1 if key == s else 0 for s in SUBS]
            sc.append(1 if (ref, alt) in TRANSITIONS else 0)
            in_dbsnp = 1 if (i_dbsnp is not None and c[i_dbsnp] not in (".", "")) else 0
            g_af = 0.0
            for ia in i_afs:
                v = c[ia]
                if v not in (".", "", "NA"):
                    try:
                        g_af = max(g_af, float(v))
                    except ValueError:
                        pass
            db = [in_dbsnp, round(g_af, 6), 1 if g_af > 0 else 0,
                  1 if (i_redi is not None and c[i_redi] not in (".", "")) else 0,
                  1 if (i_cosmic is not None and c[i_cosmic] not in (".", "")) else 0]
            vkey = (grp, c[i_chr], c[i_pos], ref, alt)
            var_vafs[vkey].append(vaf if not (isinstance(vaf, float) and np.isnan(vaf)) else None)
            rows.append(dict(cell=base, group=grp, chr=c[i_chr], pos=c[i_pos],
                             ref=ref, alt=alt,
                             label=SETNAME.get(setv, "indeterminate"),
                             set_raw=setv,
                             base=[ref_ct, alt_ct, vaf, dp, qual, gq],
                             rl=rl, sc=sc, db=db, vkey=vkey))

# Cross-cell over all rows, exactly as 21_build_feature_table.py aggregates
recs = []
for r in rows:
    vs_all = var_vafs[r["vkey"]]
    vs = [v for v in vs_all if v is not None]
    n_cells = len(vs_all)
    n_group = len(group_cells[r["group"]])
    frac = round(n_cells / n_group, 4) if n_group else np.nan
    if vs:
        m = sum(vs) / len(vs)
        s = math.sqrt(sum((v - m) ** 2 for v in vs) / len(vs)) if len(vs) > 1 else 0.0
        mean, sd = round(m, 4), round(s, 4)
    else:
        mean = sd = np.nan
    recs.append(r["base"] + r["rl"] + r["sc"] + r["db"] + [n_cells, frac, mean, sd])

cols = (["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
        + [n for n, _ in READLEVEL]
        + ["sub_" + s.replace(">", "_") for s in SUBS] + ["is_transition"]
        + ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
        + ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"])
X = pd.DataFrame(recs, columns=cols)
meta = pd.DataFrame([{k: r[k] for k in ("cell", "group", "chr", "pos", "ref", "alt", "label")}
                     for r in rows])

bundle = joblib.load(model_pkl)
feats, model, imp = bundle["features"], bundle["model"], bundle["imputer"]
missing = [f for f in feats if f not in X.columns]
if missing:
    sys.exit("missing features: %s" % missing)
p = model.predict_proba(imp.transform(X[feats].values))[:, 1]

meta["origin_score_rna_only"] = np.round(p, 4)
ind = meta["label"] == "indeterminate"

out = []


def log(s=""):
    print(s)
    out.append(s)


log("Variants the paired DNA could not decide")
log()
log("rows in the labelled files: %d" % len(meta))
log("  labelled (used for training):  %d" % (~ind).sum())
log("  indeterminate:                 %d" % ind.sum())
log()
log("by cell line:")
for g, n in meta.loc[ind, "group"].value_counts().sort_index().items():
    log("  %-10s %5d" % (g, n))
log()

pi = meta.loc[ind, "origin_score_rna_only"].values
hi = (pi > 0.9).sum()
lo = (pi < 0.1).sum()
mid = len(pi) - hi - lo
log("score distribution:")
log("  confident RNA    (>0.9)  %5d  (%.1f%%)" % (hi, 100.0 * hi / len(pi)))
log("  confident genome (<0.1)  %5d  (%.1f%%)" % (lo, 100.0 * lo / len(pi)))
log("  in between               %5d  (%.1f%%)" % (mid, 100.0 * mid / len(pi)))
log()
log("  confident in total:      %5d  (%.1f%%)" % (hi + lo, 100.0 * (hi + lo) / len(pi)))
log()

log("same for the labelled variants:")
pk = meta.loc[~ind, "origin_score_rna_only"].values
log("  confident RNA    %5d  (%.1f%%)" % ((pk > 0.9).sum(), 100.0 * (pk > 0.9).mean()))
log("  confident genome %5d  (%.1f%%)" % ((pk < 0.1).sum(), 100.0 * (pk < 0.1).mean()))
log()

path_csv = os.path.join(outdir, "indeterminate_scores.csv")
meta.loc[ind].to_csv(path_csv, index=False)
path = os.path.join(outdir, "indeterminate.txt")
with open(path, "w") as f:
    f.write("\n".join(out) + "\n")
print("\nwritten to: %s\n            %s" % (path, path_csv))
