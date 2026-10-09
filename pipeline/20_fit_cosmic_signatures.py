#!/usr/bin/env python3
"""SBS96 spectra per cell and per group, and COSMIC fitting per group."""
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerAssignment import Analyzer as Analyze
import os

ROOT = os.environ.get("ROOT", "/data")
OUT = ROOT + "/results/signatures"

# Per-cell SBS96 (45 cells)
matGen.SigProfilerMatrixGeneratorFunc(
    "percell", "GRCh38", os.path.join(OUT, "per_cell_input") + "/",
    plot=False, exome=False, bed_file=None, chrom_based=False,
    tsb_stat=False, seqInfo=True, cushion=100)

# Per-group SBS96 (3 samples: molm13, hg001, hg002)
matGen.SigProfilerMatrixGeneratorFunc(
    "pergroup", "GRCh38", os.path.join(OUT, "per_group_input") + "/",
    plot=False, exome=False, bed_file=None, chrom_based=False,
    tsb_stat=False, seqInfo=True, cushion=100)

# COSMIC v3.4 fit at the group level (+ per-variant probabilities)
group_matrix = os.path.join(OUT, "per_group_input", "output", "SBS", "pergroup.SBS96.all")
Analyze.cosmic_fit(
    samples=group_matrix,
    output=os.path.join(OUT, "group_cosmic_fit"),
    input_type="matrix", context_type="96",
    cosmic_version=3.4, exome=False, genome_build="GRCh38",
    signature_database=None, exclude_signature_subgroups=None,
    export_probabilities=True, make_plots=False,
    sample_reconstruction_plots=False, verbose=False)

print("done: per-cell SBS96, per-group SBS96, group COSMIC fit.")
