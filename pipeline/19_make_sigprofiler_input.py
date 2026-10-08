#!/usr/bin/env python3
"""Build SigProfiler input from the QC-filtered RNA SNV VCFs."""
import gzip
import glob
import os

ROOT = os.environ.get("ROOT", "/data")
IN = ROOT + "/results/rna_qc_filtered"
OUT = ROOT + "/results/step8_signatures"
PC = os.path.join(OUT, "per_cell_input")
PG = os.path.join(OUT, "per_group_input")
os.makedirs(PC, exist_ok=True)
os.makedirs(PG, exist_ok=True)

HEADER = "Project\tSample\tID\tGenome\tmut_type\tchrom\tpos_start\tpos_end\tref\talt\tType\n"
ROW = "sceSNV\t%s\t.\tGRCh38\tSNP\t%s\t%s\t%s\t%s\t%s\tRNA\n"


def group_of(base):
    if "MOLM13" in base:
        return "molm13"
    if base.startswith("ResolveOMEv2.1"):
        return "hg001"
    # 10 additional HG001 cells (ResolveOME v1.0, added 2026-08-15): same donor NA12878,
    # so same group as the v2.1 cells.
    if "NA12878" in base:
        return "hg001"
    if base.startswith("ResolveOMEv2.0-384"):
        return "hg002"
    return "other"


cell_rows, group_rows = [], []
for f in sorted(glob.glob(os.path.join(IN, "*.qc.vcf.gz"))):
    base = os.path.basename(f).replace(".qc.vcf.gz", "")
    grp = group_of(base)
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = c[0], c[1], c[3], c[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            cell_rows.append(ROW % (base, chrom, pos, pos, ref, alt))
            group_rows.append(ROW % (grp, chrom, pos, pos, ref, alt))

with open(os.path.join(PC, "per_cell.txt"), "w") as o:
    o.write(HEADER)
    o.writelines(cell_rows)
with open(os.path.join(PG, "per_group.txt"), "w") as o:
    o.write(HEADER)
    o.writelines(group_rows)

print("per_cell rows: %d, per_group rows: %d" % (len(cell_rows), len(group_rows)))
