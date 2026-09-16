"""The result table, the summary, and a record of what produced them."""

import datetime

from vartrace import __version__

COLUMNS = ["cell", "chrom", "pos", "ref", "alt", "score", "call",
           "model", "n_features"]


def rows(variants, scores, model, threshold):
    out = []
    for v, p in zip(variants, scores):
        out.append([
            v.cell, v.chrom, v.pos, v.ref, v.alt, f"{p:.4f}",
            "absent_from_genome" if p >= threshold else "in_genome",
            model.file, str(len(model.features)),
        ])
    return out


def write(path, table):
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in table:
            fh.write("\t".join(row) + "\n")


def provenance(model, why, threshold, annotation=None):
    """Enough to say later exactly what produced a given number."""
    lines = [
        f"vartrace {__version__}",
        f"date {datetime.date.today().isoformat()}",
        f"model {model.file}",
        f"model fingerprint {model.fingerprint}",
        f"features {len(model.features)}",
        f"reason {why}",
        f"threshold {threshold}",
    ]
    if annotation:
        lines.append(f"annotation {annotation['source']} release "
                     f"{annotation['release']} from {annotation['server']}")
        lines.append(f"annotation covered {annotation['with_a_record']} of "
                     f"{annotation['queried']} variants")
    else:
        lines.append("annotation taken from the input")
    return "\n".join(lines)


def summary(scores, threshold, model, why, n_cells, skipped):
    confident_absent = sum(1 for p in scores if p > 0.9)
    confident_present = sum(1 for p in scores if p < 0.1)
    middle = len(scores) - confident_absent - confident_present
    called = sum(1 for p in scores if p >= threshold)

    lines = [
        f"{len(scores)} variants from {n_cells} cell(s)",
        f"model: {model.file} ({len(model.features)} features)",
        f"why: {why}",
        f"threshold {threshold}: {called} called absent from the genome",
        f"confident absent (>0.9): {confident_absent}",
        f"confident present (<0.1): {confident_present}",
        f"uncertain: {middle}",
    ]
    if skipped:
        lines.append(f"skipped {skipped} rows that were not biallelic SNVs")
    return "\n".join(lines)
