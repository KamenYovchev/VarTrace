#!/usr/bin/env python3
"""Label each annotated RNA variant against the DNA VCF of the same cell.

Usage: 17_label_variants_by_origin.py <rna_multianno.txt> <dna_vcf(.gz)> <group> <out.tsv>
"""
import sys
import gzip

rna_file, dna_vcf, group, out_file = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
DP_MIN = 5

with open(rna_file) as fh:
    header = fh.readline().rstrip("\n").split("\t")

def idx(name):
    return header.index(name)

af_cols = [i for i, h in enumerate(header) if h == "AF"]   # gnomAD genome + exome
i_chr, i_pos = idx("Otherinfo1"), idx("Otherinfo2")
i_ref, i_alt = idx("Otherinfo4"), idx("Otherinfo5")
i_cosmic, i_redi, i_dbsnp = idx("cosmic104"), idx("rediportal"), idx("avsnp150")

# Read RNA rows + collect positions
rna_rows, positions = [], set()
with open(rna_file) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        chrom = f[i_chr] if f[i_chr].startswith("chr") else "chr" + f[i_chr]
        rna_rows.append((f, chrom, f[i_pos], f[i_ref], f[i_alt]))
        positions.add((chrom, f[i_pos]))

# DNA lookup at RNA positions only (memory-light)
op = gzip.open(dna_vcf, "rt") if dna_vcf.endswith(".gz") else open(dna_vcf)
dna = {}
with op as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        chrom, pos = f[0], f[1]
        if (chrom, pos) not in positions:
            continue
        d = dict(zip(f[8].split(":"), f[9].split(":")))
        try:
            dp_i = int(d.get("DP", "-1"))
        except ValueError:
            dp_i = -1
        dna[(chrom, pos)] = {"alts": f[4].split(","), "dp": dp_i, "gt": d.get("GT", "."),
                             "ad": d.get("AD", "."), "vaf": d.get("VAF", "."),
                             "gq": d.get("GQ", "."), "filt": f[6]}

def in_gnomad(f):
    for c in af_cols:
        v = f[c]
        if v not in (".", "", "NA"):
            try:
                if float(v) > 0:
                    return True
            except ValueError:
                pass
    return False

extra = ["set", "label", "DNA_status", "DNA_DP", "DNA_GT", "DNA_AD", "DNA_VAF", "DNA_FILTER"]
counts = {}
with open(out_file, "w") as out:
    out.write("\t".join(header + extra) + "\n")
    for f, chrom, pos, ref, alt in rna_rows:
        rec = dna.get((chrom, pos))
        g_in = in_gnomad(f)
        cosmic = f[i_cosmic] not in (".", "")
        redi = f[i_redi] not in (".", "")
        dbsnp = f[i_dbsnp] not in (".", "")
        ag = (ref == "A" and alt == "G") or (ref == "T" and alt == "C")
        setv, label, status = "indeterminate", "indeterminate", "no_dna"
        dp = gt = ad = vaf = filt = "."
        if rec:
            dp, gt, ad, vaf, filt = str(rec["dp"]), rec["gt"], rec["ad"], rec["vaf"], rec["filt"]
            matched = (alt in rec["alts"]) and (rec["gt"] not in ("0/0", "0|0", "./.", "."))
            if matched:
                setv, status = "1", "dna_variant"
                label = "germline" if (g_in or dbsnp) else ("somatic" if cosmic else "novel_DNA")
            elif rec["dp"] >= DP_MIN:
                setv, status = "2", "dna_ref_covered"
                label = "RNA_editing" if (redi or ag) else "RNA_novel"
            else:
                setv, status, label = "indeterminate", "dna_lowcov", "indeterminate"
        counts[label] = counts.get(label, 0) + 1
        out.write("\t".join(f + [setv, label, status, dp, gt, ad, vaf, filt]) + "\n")

order = ["germline", "somatic", "novel_DNA", "RNA_editing", "RNA_novel", "indeterminate"]
total = sum(counts.values())
sys.stderr.write("%s\t%s\t%s\t%d\n" % (
    group, out_file.split("/")[-1].replace(".labeled.tsv", ""),
    "\t".join(str(counts.get(k, 0)) for k in order), total))
