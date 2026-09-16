"""Reading the two kinds of input.

A plain VCF, or a table annotated with ANNOVAR. In the annotated table the VCF
fields sit in the Otherinfo columns, which is how ANNOVAR carries them over.

Only biallelic single base substitutions are kept, because that is what the
models were trained on. Everything else is counted and skipped.
"""

import os

from vartrace.variant import Variant

BASES = {"A", "C", "G", "T"}

# ANNOVAR column names, as used when the models were trained
DBSNP_COL = "avsnp150"
COSMIC_COL = "cosmic104"
REDI_COL = "rediportal"
GNOMAD_COL = "AF"


def is_snv(ref, alt):
    return ref in BASES and alt in BASES


def cell_name(path):
    name = os.path.basename(path)
    for suffix in (".labeled.tsv", ".tsv", ".vcf.gz", ".vcf", ".txt"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_vcf(path):
    variants, skipped = [], 0
    cell = cell_name(path)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                skipped += 1
                continue
            chrom, pos, _, ref, alt, qual, _, info, fmt, sample = c[:10]
            if not is_snv(ref, alt):
                skipped += 1
                continue
            variants.append(Variant(cell=cell, chrom=chrom, pos=pos, ref=ref,
                                    alt=alt, qual=qual, info=info,
                                    format=fmt, sample=sample))
    return variants, skipped


def read_annovar(path):
    variants, skipped = [], 0
    cell = cell_name(path)
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        gnomad_cols = [i for i, name in enumerate(header) if name == GNOMAD_COL]

        databases = set()
        if DBSNP_COL in idx:
            databases.add("dbsnp")
        if gnomad_cols:
            databases.add("gnomad")
        if REDI_COL in idx:
            databases.add("redi")
        if COSMIC_COL in idx:
            databases.add("cosmic")

        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < len(header):
                skipped += 1
                continue
            ref, alt = c[idx["Otherinfo4"]], c[idx["Otherinfo5"]]
            if not is_snv(ref, alt):
                skipped += 1
                continue
            variants.append(Variant(
                cell=cell,
                chrom=c[idx["Otherinfo1"]],
                pos=c[idx["Otherinfo2"]],
                ref=ref,
                alt=alt,
                qual=c[idx["Otherinfo6"]],
                info=c[idx["Otherinfo8"]],
                format=c[idx["Otherinfo9"]],
                sample=c[idx["Otherinfo10"]],
                databases=set(databases),
                dbsnp=c[idx[DBSNP_COL]] if "dbsnp" in databases else "",
                gnomad_af=[c[i] for i in gnomad_cols],
                rediportal=c[idx[REDI_COL]] if "redi" in databases else "",
                cosmic=c[idx[COSMIC_COL]] if "cosmic" in databases else "",
            ))
    return variants, skipped


def looks_annotated(path):
    with open(path) as fh:
        first = fh.readline()
    return "Otherinfo1" in first


def read(path):
    if looks_annotated(path):
        return read_annovar(path)
    return read_vcf(path)


def read_all(paths):
    """Returns every variant from every file, and how many cells that was."""
    variants, skipped = [], 0
    for path in paths:
        got, lost = read(path)
        variants.extend(got)
        skipped += lost
    return variants, len(paths), skipped
