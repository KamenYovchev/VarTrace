#!/usr/bin/env python3
"""Build the training table: features and the DNA_origin or RNA_only label.

Usage: 21_build_feature_table.py <out.csv> <labeled_dir>
"""
import sys
import os
import glob
import csv
import math
from collections import defaultdict

out_csv, labeled_dir = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(os.path.join(labeled_dir, "*.labeled.tsv")))

READLEVEL = [("fs", "FS"), ("sor", "SOR"), ("mq", "MQ"), ("mqranksum", "MQRankSum"),
             ("readposranksum", "ReadPosRankSum"), ("baseqranksum", "BaseQRankSum"),
             ("qd", "QD")]

SUBS = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def group_of(base):
    if "MOLM13" in base:
        return "molm13"
    if base.startswith("ResolveOMEv2.1"):
        return "hg001"
    # the 10 HG001 cells from ResolveOME v1.0, same donor, so same group for the
    # cross-cell features
    if "NA12878" in base:
        return "hg001"
    if base.startswith("ResolveOMEv2.0-384"):
        return "hg002"
    return "other"


def to_int(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return ""


def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d


def seqcontext(ref, alt):
    key = ref + ">" + alt
    onehot = [1 if key == s else 0 for s in SUBS]
    transition = 1 if (ref, alt) in TRANSITIONS else 0
    return onehot + [transition]


raw = []
group_cells = defaultdict(set)
var_vafs = defaultdict(list)

for f in files:
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
        # database annotation columns (ANNOVAR): avsnp150, the two gnomAD "AF" columns
        # (genome + exome), cosmic104, rediportal
        i_dbsnp = idx.get("avsnp150")
        i_afs = [i for i, n in enumerate(header) if n == "AF"]
        i_cosmic = idx.get("cosmic104")
        i_redi = idx.get("rediportal")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            setv = c[i_set]
            if setv not in ("1", "2"):
                continue
            ref, alt = c[i_ref], c[i_alt]
            d = dict(zip(c[i_fmt].split(":"), c[i_smp].split(":")))
            ref_ct = alt_ct = ""
            ad = d.get("AD", ".")
            if "," in ad:
                p = ad.split(",")
                ref_ct, alt_ct = to_int(p[0]), to_int(p[1])
            dp = to_int(d.get("DP", "."))
            gq = to_int(d.get("GQ", "."))
            vaf = d.get("VAF")
            if vaf in (None, ".", ""):
                if ref_ct != "" and alt_ct != "" and (ref_ct + alt_ct) > 0:
                    vaf = round(alt_ct / (ref_ct + alt_ct), 4)
                else:
                    vaf = ""
            else:
                try:
                    vaf = float(vaf)
                except ValueError:
                    vaf = ""
            try:
                qual = float(c[i_qual])
            except ValueError:
                qual = ""
            info = parse_info(c[i_info])
            rl = [info.get(key, "") for _, key in READLEVEL]
            sc = seqcontext(ref, alt)
            # database features
            in_dbsnp = 1 if (i_dbsnp is not None and c[i_dbsnp] not in (".", "")) else 0
            g_af = 0.0
            for ia in i_afs:
                v = c[ia]
                if v not in (".", "", "NA"):
                    try:
                        g_af = max(g_af, float(v))
                    except ValueError:
                        pass
            in_gnomad = 1 if g_af > 0 else 0
            in_cosmic = 1 if (i_cosmic is not None and c[i_cosmic] not in (".", "")) else 0
            in_redi = 1 if (i_redi is not None and c[i_redi] not in (".", "")) else 0
            db = [in_dbsnp, round(g_af, 6), in_gnomad, in_redi, in_cosmic]
            label = "DNA_origin" if setv == "1" else "RNA_only"
            vkey = (grp, c[i_chr], c[i_pos], ref, alt)
            var_vafs[vkey].append(float(vaf) if vaf != "" else None)
            raw.append((base, grp, c[i_chr], c[i_pos], ref, alt,
                        [ref_ct, alt_ct, vaf, dp, qual, gq] + rl + sc + db, label))

rows = []
for base, grp, chrom, pos, ref, alt, feats, label in raw:
    vkey = (grp, chrom, pos, ref, alt)
    vafs = [v for v in var_vafs[vkey] if v is not None]
    n_cells = len(var_vafs[vkey])
    n_group = len(group_cells[grp])
    frac = round(n_cells / n_group, 4) if n_group else ""
    if vafs:
        mean = sum(vafs) / len(vafs)
        sd = math.sqrt(sum((v - mean) ** 2 for v in vafs) / len(vafs)) if len(vafs) > 1 else 0.0
        mean, sd = round(mean, 4), round(sd, 4)
    else:
        mean = sd = ""
    rows.append([base, grp, chrom, pos, ref, alt] + feats
                + [n_cells, frac, mean, sd] + [label])

sub_cols = ["sub_" + s.replace(">", "_") for s in SUBS] + ["is_transition"]
db_cols = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]
cols = (["cell", "group", "chr", "pos", "ref", "alt",
         "ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
        + [name for name, _ in READLEVEL] + sub_cols + db_cols
        + ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
        + ["label"])
with open(out_csv, "w", newline="") as o:
    w = csv.writer(o)
    w.writerow(cols)
    w.writerows(rows)

n = len(rows)
n1 = sum(1 for r in rows if r[-1] == "DNA_origin")
sys.stderr.write("total rows: %d\n" % n)
sys.stderr.write("DNA_origin (Set1): %d\n" % n1)
sys.stderr.write("RNA_only   (Set2): %d\n" % (n - n1))
sys.stderr.write("unique variants: %d\n" % len(var_vafs))
sys.stderr.write("database features: %s\n" % ", ".join(db_cols))
sys.stderr.write("written to: %s\n" % out_csv)
