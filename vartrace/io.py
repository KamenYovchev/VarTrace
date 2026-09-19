"""Reading the two kinds of input.

A plain VCF, or a table annotated with ANNOVAR. In the annotated table the VCF
fields sit in the Otherinfo columns, which is how ANNOVAR carries them over.

A VCF may hold one sample or many. In single-cell work the second is common,
one column per cell barcode, and then every column is treated as its own cell.
A cell is only counted for a variant when its genotype or its allele depth says
it carries it, otherwise every cell would appear to carry everything.

Only biallelic single base substitutions are kept, because that is what the
models were trained on. Everything else is counted and skipped.
"""

import os

from vartrace.variant import Variant

BASES = {"A", "C", "G", "T"}

# chromosome 1 is a different length in the two assemblies, so one header line
# settles which one a VCF was called against
CHR1_LENGTH = {248956422: "GRCh38", 249250621: "GRCh37"}

# ANNOVAR column names, as used when the models were trained
DBSNP_COL = "avsnp150"
COSMIC_COL = "cosmic104"
REDI_COL = "rediportal"
GNOMAD_COL = "AF"


class UnknownFormat(ValueError):
    """The file is a table, but not one this tool knows how to read."""


def is_snv(ref, alt):
    # some callers write soft-masked positions in lower case, and those are
    # ordinary substitutions rather than something to throw away
    return ref.upper() in BASES and alt.upper() in BASES


def genome_build(path):
    """Which assembly a VCF was called against, or None if it cannot be told.

    The models and the Ensembl lookup both work on GRCh38. Variants given on
    GRCh37 sit at positions that mean something else, and the whole run comes
    out quietly wrong, so it is worth saying which one arrived. An ANNOVAR
    table carries no header and leaves the question open; None says so rather
    than guessing.
    """
    if looks_annotated(path):
        return None
    with open(path) as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if not line.startswith("##contig="):
                continue
            fields = dict(part.split("=", 1) for part in
                          line.strip().rstrip(">")[len("##contig=<"):].split(",")
                          if "=" in part)
            if fields.get("ID") in ("1", "chr1"):
                try:
                    return CHR1_LENGTH.get(int(fields.get("length", 0)))
                except ValueError:
                    return None
    return None


def cell_name(path):
    name = os.path.basename(path)
    for suffix in (".labeled.tsv", ".tsv", ".vcf.gz", ".vcf", ".txt"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def carries(fmt, value):
    """Does this sample actually have the variant."""
    fields = dict(zip(fmt.split(":"), value.split(":")))
    gt = fields.get("GT", "")
    if gt:
        return "1" in gt.replace("|", "/").split("/")
    ad = fields.get("AD", "")
    if "," in ad:
        try:
            return int(ad.split(",")[1]) > 0
        except ValueError:
            return False
    return False


def read_vcf(path):
    variants, skipped = [], 0
    samples = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                skipped += 1
                continue
            chrom, pos, _, ref, alt, qual, _, info, fmt = c[:9]
            if not is_snv(ref, alt):
                skipped += 1
                continue
            # upper case from here on, so the substitution features and the
            # Ensembl queries see one spelling and not two
            ref, alt = ref.upper(), alt.upper()
            for i, value in enumerate(c[9:]):
                if len(samples) > 1:
                    if not carries(fmt, value):
                        continue
                    cell = samples[i]
                else:
                    cell = cell_name(path)
                variants.append(Variant(cell=cell, chrom=chrom, pos=pos, ref=ref,
                                        alt=alt, qual=qual, info=info,
                                        format=fmt, sample=value))
                if len(samples) <= 1:
                    break
    return variants, skipped, max(len(samples), 1)


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
            ref, alt = ref.upper(), alt.upper()
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
    return variants, skipped, 1


def looks_annotated(path):
    with open(path) as fh:
        first = fh.readline()
    return "Otherinfo1" in first


def read(path):
    if looks_annotated(path):
        return read_annovar(path)
    with open(path) as fh:
        first = fh.readline()
    # a VCF starts with its header; a first line of plain columns means some
    # other annotator's table, which would otherwise be read as an empty VCF
    if first and not first.startswith("#"):
        raise UnknownFormat(
            f"{os.path.basename(path)} looks like a table, but not one from "
            "ANNOVAR. Tables from other annotators are not read. Export a VCF "
            "from them instead, which works whoever made it.")
    return read_vcf(path)


def read_all(paths):
    """Every variant from every file, and how many cells that was.

    One file per cell gives one cell each; a VCF with many samples gives as
    many cells as it has samples.
    """
    variants, skipped, cells = [], 0, 0
    for path in paths:
        got, lost, n = read(path)
        variants.extend(got)
        skipped += lost
        cells += n
    return variants, cells, skipped
