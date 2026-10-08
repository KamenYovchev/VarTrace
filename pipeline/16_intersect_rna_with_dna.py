#!/usr/bin/env python3
"""Left-join RNA annotated variants with the matching DNA VCF by chr+pos+ref+alt.

Usage: 16_intersect_rna_with_dna.py <rna_multianno.txt> <dna_vcf.gz> <output.tsv>
"""
import sys, gzip

rna_file, dna_vcf, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

# 1) RNA header, then locate the original-VCF columns (Otherinfo1=chrom ...
with open(rna_file) as fh:
    header = fh.readline().rstrip("\n").split("\t")
try:
    i_chr = header.index("Otherinfo1"); i_pos = header.index("Otherinfo2")
    i_ref = header.index("Otherinfo4"); i_alt = header.index("Otherinfo5")
except ValueError:                       # fallback: first 5 cols (Chr is chr-stripped)
    i_chr, i_pos, i_ref, i_alt = 0, 1, 3, 4

# 2) read RNA rows + collect variant keys (chr-prefixed to match the DNA VCF)
rna_rows, rna_keys = [], set()
with open(rna_file) as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        chrom = f[i_chr] if f[i_chr].startswith("chr") else "chr" + f[i_chr]
        key = (chrom, f[i_pos], f[i_ref], f[i_alt])
        rna_rows.append((f, key)); rna_keys.add(key)

# 3) scan DNA VCF once; keep only records whose chr:pos:ref:alt matches an RNA variant
op = gzip.open(dna_vcf, "rt") if dna_vcf.endswith(".gz") else open(dna_vcf)
dna = {}
with op as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        chrom, pos, ref, alts, filt, fmt, samp = f[0], f[1], f[3], f[4], f[6], f[8], f[9]
        for alt in alts.split(","):
            key = (chrom, pos, ref, alt)
            if key in rna_keys:
                d = dict(zip(fmt.split(":"), samp.split(":")))
                dna[key] = {"GT": d.get("GT", "."), "DP": d.get("DP", "."),
                            "AD": d.get("AD", "."), "GQ": d.get("GQ", "."),
                            "VAF": d.get("VAF", "."), "FILTER": filt}

# 4) write RNA rows + appended DNA columns
extra = ["DNA_matched", "DNA_FILTER", "DNA_GT", "DNA_DP", "DNA_AD", "DNA_GQ", "DNA_VAF"]
n_match = 0
with open(out_file, "w") as out:
    out.write("\t".join(header + extra) + "\n")
    for f, key in rna_rows:
        if key in dna:
            m = dna[key]; n_match += 1
            e = ["yes", m["FILTER"], m["GT"], m["DP"], m["AD"], m["GQ"], m["VAF"]]
        else:
            e = ["no", ".", ".", ".", ".", ".", "."]
        out.write("\t".join(f + e) + "\n")

sys.stderr.write("%s : %d RNA variants, %d matched in DNA\n"
                 % (out_file, len(rna_rows), n_match))
