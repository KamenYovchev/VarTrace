#!/usr/bin/env python3
"""How many HG001 and HG002 labels does Genome in a Bottle contradict?

Usage: 40_check_labels_giab.py <dataset_sig.csv> <giab ref dir> <outdir>
"""
import bisect
import os
import shlex
import subprocess
import sys
from collections import defaultdict

import pandas as pd

csv_in, ref_dir, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)
BCFTOOLS = shlex.split(os.environ.get("BCFTOOLS", "bcftools"))
LINES = {
    "hg001": ("HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
              "HG001_GRCh38_1_22_v4.2.1_benchmark.bed"),
    "hg002": ("HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
              "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"),
}


def read_bed(path):
    iv = defaultdict(list)
    for line in open(path):
        c, s, e = line.split()[:3]
        iv[c].append((int(s), int(e)))
    starts = {c: [s for s, _ in sorted(v)] for c, v in iv.items()}
    ends = {c: [e for _, e in sorted(v)] for c, v in iv.items()}
    return starts, ends


def in_bed(starts, ends, c, pos):
    # BED starts at 0, VCF positions at 1
    if c not in starts:
        return False
    i = bisect.bisect_right(starts[c], pos - 1) - 1
    return i >= 0 and pos <= ends[c][i]


def read_truth(path):
    q = subprocess.run(BCFTOOLS + ["query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n", path],
                       capture_output=True, text=True, check=True).stdout
    snv, anypos = {}, set()
    for line in q.splitlines():
        c, p, r, a, gt = line.split("\t")
        p = int(p)
        anypos.add((c, p))
        for alt in a.split(","):
            if len(r) == 1 and len(alt) == 1:
                snv[(c, p, r, alt)] = gt
    return snv, anypos


df = pd.read_csv(csv_in, usecols=["cell", "group", "chr", "pos", "ref", "alt", "label"])
rows = []
for line, (vcf, bed) in LINES.items():
    starts, ends = read_bed(os.path.join(ref_dir, bed))
    snv, anypos = read_truth(os.path.join(ref_dir, vcf))
    sub = df[df.group == line].copy()
    keys = list(zip(sub.chr, sub.pos, sub.ref, sub.alt))
    sub["in_region"] = [in_bed(starts, ends, c, p) for c, p, _, _ in keys]
    sub["giab_match"] = [k in snv for k in keys]
    sub["giab_gt"] = [snv.get(k, "") for k in keys]
    sub["giab_other_at_pos"] = [(c, p) in anypos and not m
                                for (c, p, _, _), m in zip(keys, sub.giab_match)]
    sub.to_csv(os.path.join(outdir, "%s_rows.tsv" % line), sep="\t", index=False)

    for label in ["RNA_only", "DNA_origin"]:
        s = sub[sub.label == label]
        r = s[s.in_region]
        u = r.drop_duplicates(["chr", "pos", "ref", "alt"])
        rows.append({
            "line": line, "label": label,
            "rows": len(s), "rows_in_region": len(r),
            "rows_match_giab": int(r.giab_match.sum()),
            "share_match": round(r.giab_match.mean(), 4) if len(r) else None,
            "variants_in_region": len(u),
            "variants_match_giab": int(u.giab_match.sum()),
            "share_variants_match": round(u.giab_match.mean(), 4) if len(u) else None,
        })

res = pd.DataFrame(rows)
res.to_csv(os.path.join(outdir, "summary.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(res.to_string(index=False))
