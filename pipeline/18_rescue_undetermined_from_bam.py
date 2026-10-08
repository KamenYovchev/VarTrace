#!/usr/bin/env python3
"""Rescue no_dna RNA positions using DNA BAM depth (streamed over http).

Usage: 18_rescue_undetermined_from_bam.py <in.labeled.tsv> <out.labeled.tsv> <bam_url|none> <samtools>
"""
import sys
import os
import subprocess

in_tsv, out_tsv, bam_url, samtools = sys.argv[1:5]
DP_MIN = 5

with open(in_tsv) as fh:
    lines = fh.read().rstrip("\n").split("\n")
header = lines[0].split("\t")

def idx(n):
    return header.index(n)

i_chr, i_pos = idx("Otherinfo1"), idx("Otherinfo2")
i_ref, i_alt = idx("Otherinfo4"), idx("Otherinfo5")
i_redi = idx("rediportal")
n = len(header)
i_set, i_label, i_status, i_dp = n - 8, n - 7, n - 6, n - 5   # appended extra cols

rows = [l.split("\t") for l in lines[1:]]

targets = []
for ri, r in enumerate(rows):
    if r[i_status] == "no_dna":
        chrom = r[i_chr] if r[i_chr].startswith("chr") else "chr" + r[i_chr]
        targets.append((ri, chrom, int(r[i_pos])))

depth = {}
if targets and bam_url != "NONE":
    bedfile = out_tsv + ".nodna.bed"
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(bedfile, "w") as b:
        for _, chrom, pos in sorted(targets, key=lambda t: (t[1], t[2])):
            b.write("%s\t%d\t%d\n" % (chrom, pos - 1, pos))
    res = subprocess.run([samtools, "bedcov", bedfile, bam_url],
                         capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write("SAMTOOLS_ERR %s: %s\n" % (out_tsv, res.stderr[:300]))
    for line in res.stdout.rstrip("\n").split("\n"):
        if not line:
            continue
        c = line.split("\t")
        depth[(c[0], int(c[2]))] = int(c[-1])
    os.remove(bedfile)

for ri, chrom, pos in targets:
    d = depth.get((chrom, pos))
    if d is None:
        continue
    r = rows[ri]
    r[i_dp] = str(d)
    if d >= DP_MIN:
        ref, alt = r[i_ref], r[i_alt]
        ag = (ref == "A" and alt == "G") or (ref == "T" and alt == "C")
        redi = r[i_redi] not in (".", "")
        r[i_set] = "2"
        r[i_status] = "dna_ref_covered_bam"
        r[i_label] = "RNA_editing" if (redi or ag) else "RNA_novel"

counts = {}
for r in rows:
    counts[r[i_label]] = counts.get(r[i_label], 0) + 1

os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
with open(out_tsv, "w") as out:
    out.write("\t".join(header) + "\n")
    for r in rows:
        out.write("\t".join(r) + "\n")

order = ["germline", "somatic", "novel_DNA", "RNA_editing", "RNA_novel", "indeterminate"]
sys.stderr.write("%s\t%s\t%d\n" % (
    os.path.basename(out_tsv).replace(".labeled.tsv", ""),
    "\t".join(str(counts.get(k, 0)) for k in order), sum(counts.values())))
