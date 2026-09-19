"""The 35 features the models expect, computed as in the training pipeline.

Anything that cannot be computed stays None. None is not zero. A variant with
no gnomAD column is not a variant absent from gnomAD, and the choice of model
depends on that difference.
"""

import math

import numpy as np

SUBS = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
        "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
READLEVEL = [("fs", "FS"), ("sor", "SOR"), ("mq", "MQ"), ("mqranksum", "MQRankSum"),
             ("readposranksum", "ReadPosRankSum"), ("baseqranksum", "BaseQRankSum"),
             ("qd", "QD")]

BASELINE = ["ref_count", "alt_count", "vaf", "dp", "qual", "gq"]
SEQCTX = ["sub_" + s.replace(">", "_") for s in SUBS] + ["is_transition"]
CROSSCELL = ["n_cells_variant", "frac_cells", "vaf_mean_cells", "vaf_sd_cells"]
DATABASE = ["in_dbsnp", "gnomad_af", "in_gnomad", "in_rediportal", "in_cosmic"]

ALL_FEATURES = (BASELINE + [name for name, _ in READLEVEL] + SEQCTX
                + CROSSCELL + DATABASE)

DB_GROUPS = {"dbsnp", "gnomad", "redi", "cosmic"}


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def is_empty(x):
    return x is None or x in (".", "", "NA")


def parse_info(info):
    out = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
    return out


def database_features(v):
    f = dict.fromkeys(DATABASE)
    if "dbsnp" in v.databases:
        f["in_dbsnp"] = 0 if is_empty(v.dbsnp) else 1
    if "gnomad" in v.databases:
        # ANNOVAR writes one AF column per gnomAD source and the training
        # pipeline took the highest of them
        af = 0.0
        for x in v.gnomad_af:
            value = to_float(x)
            if value is not None:
                af = max(af, value)
        f["gnomad_af"] = round(af, 6)
        f["in_gnomad"] = 1 if af > 0 else 0
    if "redi" in v.databases:
        f["in_rediportal"] = 0 if is_empty(v.rediportal) else 1
    if "cosmic" in v.databases:
        f["in_cosmic"] = 0 if is_empty(v.cosmic) else 1
    return f


def variant_features(v):
    f = {}
    sample = dict(zip(v.format.split(":"), v.sample.split(":")))

    ref_count = alt_count = None
    ad = sample.get("AD", "")
    if "," in ad:
        parts = ad.split(",")
        ref_count, alt_count = to_float(parts[0]), to_float(parts[1])
    f["ref_count"] = ref_count
    f["alt_count"] = alt_count

    vaf = sample.get("VAF")
    if is_empty(vaf):
        if ref_count is not None and alt_count is not None and ref_count + alt_count > 0:
            vaf = round(alt_count / (ref_count + alt_count), 4)
        else:
            vaf = None
    else:
        vaf = to_float(vaf)
    f["vaf"] = vaf
    f["dp"] = to_float(sample.get("DP"))
    f["gq"] = to_float(sample.get("GQ"))
    f["qual"] = to_float(v.qual)

    info = parse_info(v.info)
    for name, key in READLEVEL:
        f[name] = to_float(info.get(key))

    sub = v.ref + ">" + v.alt
    for s in SUBS:
        f["sub_" + s.replace(">", "_")] = 1 if s == sub else 0
    f["is_transition"] = 1 if (v.ref, v.alt) in TRANSITIONS else 0

    for name in CROSSCELL:
        f[name] = None
    f.update(database_features(v))

    v.features = f
    return f


def add_cross_cell(variants, n_cells):
    """Needs every cell at once, so it runs after all files are read.

    A variant seen in several cells gets the same four values in each of them.
    """
    seen = {}
    for v in variants:
        seen.setdefault(v.key, []).append(v.features["vaf"])

    for v in variants:
        vafs = seen[v.key]
        known = [x for x in vafs if x is not None]
        f = v.features
        f["n_cells_variant"] = len(vafs)
        f["frac_cells"] = round(len(vafs) / n_cells, 4) if n_cells else None
        if known:
            mean = sum(known) / len(known)
            if len(known) > 1:
                sd = math.sqrt(sum((x - mean) ** 2 for x in known) / len(known))
            else:
                sd = 0.0
            f["vaf_mean_cells"] = round(mean, 4)
            f["vaf_sd_cells"] = round(sd, 4)


def compute(variants, n_cells):
    for v in variants:
        variant_features(v)
    if n_cells > 1:
        add_cross_cell(variants, n_cells)


def matrix(variants, names):
    rows = []
    for v in variants:
        f = v.features
        rows.append([np.nan if f[n] is None else float(f[n]) for n in names])
    return np.array(rows, dtype=float)


def without_allele_counts(variants):
    """How many variants came without the allele fraction.

    It is the feature the models lean on most after the databases, and a file
    that does not carry AD or VAF cannot produce it. Scores from such a file
    would look ordinary and mean nothing, so the caller warns about it.
    """
    return sum(1 for v in variants if v.features["vaf"] is None)


def without_read_level(variants):
    """How many variants came without any of the read-level annotations.

    FS, SOR, MQ and the rest of them are written by GATK. A VCF from
    DeepVariant, bcftools or Clair3 carries none, and unlike the databases
    there is no model trained for their absence, so seven of the thirty-five
    features arrive empty and nothing in the numbers shows it.
    """
    names = [name for name, _ in READLEVEL]
    return sum(1 for v in variants if all(v.features[n] is None for n in names))


def available(variants):
    """Which feature groups this input can support."""
    f = variants[0].features
    groups = set()
    if f["n_cells_variant"] is not None:
        groups.add("cells")
    if f["in_dbsnp"] is not None:
        groups.add("dbsnp")
    if f["gnomad_af"] is not None:
        groups.add("gnomad")
    if f["in_rediportal"] is not None:
        groups.add("redi")
    if f["in_cosmic"] is not None:
        groups.add("cosmic")
    return groups
