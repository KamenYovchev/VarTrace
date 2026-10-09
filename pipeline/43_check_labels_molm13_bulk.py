#!/usr/bin/env python3
"""How many MOLM13 labels does the bulk DNA from step 42 contradict?

Usage: 43_check_labels_molm13_bulk.py <dataset_sig.csv> <bulk.vcf.gz> <bulk.bam> <reference.fa> <outdir>
"""
import os
import shlex
import subprocess
import sys

import pandas as pd

csv_in, vcf, bam, ref, outdir = sys.argv[1:6]
os.makedirs(outdir, exist_ok=True)
BCFTOOLS = shlex.split(os.environ.get("BCFTOOLS", "bcftools"))
MIN_DP = 10

df = pd.read_csv(csv_in, usecols=["cell", "group", "chr", "pos", "ref", "alt", "label"])
df = df[df.group == "molm13"].copy()
positions = os.path.join(outdir, "molm13_positions.tsv")
df[["chr", "pos"]].drop_duplicates().sort_values(["chr", "pos"]).to_csv(
    positions, sep="\t", header=False, index=False)


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


calls_vcf = os.path.join(outdir, "bulk_calls.vcf.gz")
run(BCFTOOLS + ["view", "-R", positions, vcf, "-Oz", "-o", calls_vcf])
calls = run(BCFTOOLS + ["query", "-f",
             "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL[\t%GT\t%GQ\t%AD\t%DP\t%VAF]\n", calls_vcf])
pile = subprocess.run(
    f"{shlex.join(BCFTOOLS)} mpileup -f {ref} -R {positions} -a AD,DP -q 20 -Q 20 --max-depth 2000 {bam} 2>/dev/null"
    f" | {shlex.join(BCFTOOLS)} query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD\\t%DP]\\n'",
    shell=True, capture_output=True, text=True, check=True).stdout

called = {}
for line in calls.splitlines():
    c, p, r, alts, flt, qual, gt, gq, ad, dp, vaf = line.split("\t")
    alts = alts.split(",")
    for i in set(gt.replace("|", "/").split("/")) - {"0", "."}:
        called[(c, int(p), r, alts[int(i) - 1])] = (flt, gt, vaf, gq)

depth, altreads = {}, {}
for line in pile.splitlines():
    c, p, r, alts, ad, dp = line.split("\t")
    ads = [int(x) for x in ad.split(",")]
    depth[(c, int(p))] = sum(ads)
    for a, n in zip(alts.split(","), ads[1:]):
        if a != "<*>":
            altreads[(c, int(p), r, a)] = n

k = list(zip(df.chr, df.pos, df.ref, df.alt))
df["bulk_dp"] = [depth.get((c, p), 0) for c, p, _, _ in k]
df["bulk_alt_reads"] = [altreads.get(x, 0) for x in k]
df["bulk_call"] = [called.get(x, ("", "", "", ""))[0] for x in k]
df["bulk_gt"] = [called.get(x, ("", "", "", ""))[1] for x in k]
df["bulk_vaf"] = [called.get(x, ("", "", "", ""))[2] for x in k]
df["bulk_gq"] = pd.to_numeric([called.get(x, ("", "", "", ""))[3] for x in k], errors="coerce")
df["covered"] = df.bulk_dp >= MIN_DP
df["bulk_pass"] = df.bulk_call == "PASS"
df["wrong"] = (df.label == "RNA_only") & df.covered & df.bulk_pass
df["alt_seen_not_called"] = df.covered & ~df.bulk_pass & (df.bulk_alt_reads >= 2)
df.to_csv(os.path.join(outdir, "molm13_rows.tsv"), sep="\t", index=False)

out = []
for label in ["RNA_only", "DNA_origin"]:
    s = df[df.label == label]
    c = s[s.covered]
    u = c.drop_duplicates(["chr", "pos", "ref", "alt"])
    out.append({
        "label": label, "rows": len(s), "rows_covered": len(c),
        "rows_bulk_pass": int(c.bulk_pass.sum()),
        "share_bulk_pass": round(c.bulk_pass.mean(), 4) if len(c) else None,
        "variants_covered": len(u), "variants_bulk_pass": int(u.bulk_pass.sum()),
        "share_variants": round(u.bulk_pass.mean(), 4) if len(u) else None,
        "rows_alt_seen_not_called": int(c.alt_seen_not_called.sum()),
    })
res = pd.DataFrame(out)
res.to_csv(os.path.join(outdir, "summary.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(res.to_string(index=False))

w = df[df.wrong]
u = w.drop_duplicates(["chr", "pos", "ref", "alt"])
print("\nwrong RNA_only by bulk genotype:", w.bulk_gt.value_counts().to_dict())
a = df[(df.label == "RNA_only") & df.covered]
r = a.groupby("cell").wrong.mean()
print("per cell, %d cells: min %.3f, median %.3f, max %.3f" % (len(r), r.min(), r.median(), r.max()))
print("borderline among %d wrong variants: depth 10-14: %d, GQ below 20: %d"
      % (len(u), (u.bulk_dp < 15).sum(), (u.bulk_gq < 20).sum()))
clear = u[(u.bulk_dp >= 15) & ~(u.bulk_gq < 20)]
print("clear cases: %d, median depth %.0f" % (len(clear), clear.bulk_dp.median()))
print("RNA_only not covered (depth < %d): %d" % (MIN_DP, ((df.label == "RNA_only") & ~df.covered).sum()))
