#!/usr/bin/env python3
"""Add variant-pattern (mutational-signature) features to the model dataset.

Usage: 22_add_signature_features.py <in.csv> <ref.fa> <decomposed.txt> <out.csv> [samtools]
"""
import sys
import csv
import subprocess
import tempfile
import os

in_csv, ref_fa, decomp_txt, out_csv = sys.argv[1:5]
samtools = sys.argv[5] if len(sys.argv) > 5 else "samtools"

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def revcomp(s):
    return "".join(COMP.get(b, "N") for b in reversed(s))


# Read dataset
with open(in_csv) as fh:
    r = csv.reader(fh)
    header = next(r)
    rows = list(r)
i_chr, i_pos = header.index("chr"), header.index("pos")
i_ref, i_alt = header.index("ref"), header.index("alt")
i_grp = header.index("group")

# Unique loci, then regions file
loci = []
seen = set()
for row in rows:
    key = (row[i_chr], row[i_pos])
    if key not in seen:
        seen.add(key)
        loci.append(key)

with tempfile.NamedTemporaryFile("w", suffix=".regions", delete=False) as tf:
    regions_path = tf.name
    for chrom, pos in loci:
        p = int(pos)
        tf.write("%s:%d-%d\n" % (chrom, max(1, p - 1), p + 1))

# Fetch contexts with samtools faidx
res = subprocess.run([samtools, "faidx", ref_fa, "-r", regions_path],
                     capture_output=True, text=True)
os.unlink(regions_path)
if res.returncode != 0:
    sys.stderr.write("samtools faidx failed: %s\n" % res.stderr[:500])
    sys.exit(1)

contexts = {}
name, seq = None, []
for line in res.stdout.splitlines() + [">END"]:
    if line.startswith(">"):
        if name is not None:
            contexts[name] = "".join(seq).upper()
        name, seq = line[1:].strip(), []
    else:
        seq.append(line.strip())

ctx_by_locus = {}
for chrom, pos in loci:
    p = int(pos)
    key = "%s:%d-%d" % (chrom, max(1, p - 1), p + 1)
    ctx_by_locus[(chrom, pos)] = contexts.get(key, "")

# SBS96 mutation type per row
def sbs96(ctx, ref, alt):
    if len(ctx) != 3 or "N" in ctx or ctx[1] != ref:
        return None
    if ref in ("A", "G"):                       # pyrimidine-center convention
        ctx = revcomp(ctx)
        ref, alt = COMP[ref], COMP[alt]
    return "%s[%s>%s]%s" % (ctx[0], ref, alt, ctx[2])


# Read decomposition: per (sample/group, mutation type), then sig probabilities
with open(decomp_txt) as fh:
    dr = csv.reader(fh, delimiter="\t")
    dheader = next(dr)
    signatures = dheader[2:]                    # col0 sample, col1 mutation type
    probs = {}
    for line in dr:
        probs[(line[0], line[1])] = [float(x) for x in line[2:]]

sig_cols = ["sig_" + s for s in signatures]

# Write output
n_matched = n_ctx_fail = 0
with open(out_csv, "w", newline="") as o:
    w = csv.writer(o)
    lab = header.index("label")
    w.writerow(header[:lab] + sig_cols + [header[lab]])
    for row in rows:
        ctx = ctx_by_locus.get((row[i_chr], row[i_pos]), "")
        mt = sbs96(ctx, row[i_ref], row[i_alt])
        pv = probs.get((row[i_grp], mt)) if mt else None
        if pv is None:
            pv = [0.0] * len(signatures)
            n_ctx_fail += 1
        else:
            n_matched += 1
        w.writerow(row[:lab] + ["%.4f" % x for x in pv] + [row[lab]])

sys.stderr.write("signatures found in decomposition: %s\n" % ", ".join(signatures))
sys.stderr.write("rows with signature probabilities: %d\n" % n_matched)
sys.stderr.write("rows without (context/type missing): %d\n" % n_ctx_fail)
sys.stderr.write("written to: %s\n" % out_csv)
