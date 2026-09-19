"""Is annotation fetched from Ensembl as good as annotation done with ANNOVAR?

The models were trained on the database columns ANNOVAR writes. A user with no
ANNOVAR gets dbSNP and gnomAD from Ensembl instead, and until now nothing had
checked that the two agree. The labelled tables say which variants really are
absent from the genome, so the question can be settled on the result and not
only on how far the values drift.

The same variants are scored three ways:

    annovar_full   every database column ANNOVAR wrote
    annovar_two    only dbSNP and gnomAD, REDIportal and COSMIC withheld
    ensembl_two    dbSNP and gnomAD as Ensembl gave them, the other two withheld

annovar_two against ensembl_two is the comparison that matters. Both reach the
same model with the same feature groups present, so any difference between them
comes from the values alone. annovar_full is there to show what the two
withheld databases were worth in the first place.

The answers from Ensembl are read from the file fetch_ensembl.py wrote.

    PYTHONPATH=/data/vartrace /data/conda/envs/ml-cpu/bin/python \\
        /data/vartrace/tests/compare_annotation_sources.py \\
        /data/results/step9_labeled_bamfix /data/results/model_dataset_sig.csv \\
        /data/results/ensembl_answers.tsv

One caveat on the numbers. These variants are the ones the models were trained
on, so the figures here are higher than the cross-validated ones in the model
descriptions and are not comparable with them. They are comparable with each
other, which is the whole point.
"""

import csv
import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from vartrace import features, io, model as models

labeled_dir, training_csv, ensembl_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

POSITIVE = "RNA_only"


def group_of(cell):
    """Cross-cell features were computed per cell line during training."""
    if "MOLM13" in cell:
        return "molm13"
    if cell.startswith("ResolveOMEv2.1") or "NA12878" in cell:
        return "hg001"
    if cell.startswith("ResolveOMEv2.0-384"):
        return "hg002"
    return "other"


expected = pd.read_csv(training_csv).dropna(subset=["label"])
label_of = {}
for r in expected.itertuples(index=False):
    label_of[(r.cell, r.chr, str(r.pos), r.ref, r.alt)] = 1 if r.label == POSITIVE else 0
print(f"{len(label_of)} labelled rows, {sum(label_of.values())} of them {POSITIVE}")

from_ensembl = {}
with open(ensembl_tsv) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        from_ensembl[(row["chrom"], row["pos"], row["ref"], row["alt"])] = (
            row["rs"], row["gnomad_af"])
print(f"{len(from_ensembl)} variants answered for by Ensembl")

variants = []
for path in sorted(glob.glob(os.path.join(labeled_dir, "*.labeled.tsv"))):
    got, _, _ = io.read_annovar(path)
    variants.extend(v for v in got
                    if (v.cell, v.chrom, v.pos, v.ref, v.alt) in label_of)
print(f"{len(variants)} rows read back from the labelled tables")

by_line = {}
for v in variants:
    by_line.setdefault(group_of(v.cell), []).append(v)

# what ANNOVAR said, kept aside before anything is overwritten
as_annovar = {id(v): (set(v.databases), v.dbsnp, list(v.gnomad_af)) for v in variants}

loaded = models.load_all()
y = np.array([label_of[(v.cell, v.chrom, v.pos, v.ref, v.alt)] for v in variants])


def score_everything(prepare):
    """Apply prepare() to every variant, then score exactly as the tool would."""
    for v in variants:
        prepare(v)
    for line, group in by_line.items():
        for v in group:
            features.variant_features(v)
        features.add_cross_cell(group, len({c.cell for c in group}))

    chosen, why = models.choose(loaded, features.available(variants))
    order = {id(v): i for i, v in enumerate(variants)}
    out = np.empty(len(variants))
    for line, group in by_line.items():
        for v, p in zip(group, chosen.predict(group)):
            out[order[id(v)]] = p
    return chosen, why, out


def as_it_was(v):
    dbs, rs, af = as_annovar[id(v)]
    v.databases, v.dbsnp, v.gnomad_af = set(dbs), rs, list(af)


def two_from_annovar(v):
    as_it_was(v)
    v.databases -= {"redi", "cosmic"}


def two_from_ensembl(v):
    rs, af = from_ensembl.get((v.chrom, v.pos, v.ref, v.alt), (".", "0.0"))
    v.databases = {"dbsnp", "gnomad"}
    v.dbsnp, v.gnomad_af = rs, [af]


runs = {}
for name, prepare in (("annovar_full", as_it_was),
                      ("annovar_two", two_from_annovar),
                      ("ensembl_two", two_from_ensembl)):
    chosen, why, scores = score_everything(prepare)
    runs[name] = scores
    print(f"\n{name}")
    print(f"  model {chosen.file}")
    print(f"  {why}")
    print(f"  PR-AUC  {average_precision_score(y, scores):.4f}")
    print(f"  ROC-AUC {roc_auc_score(y, scores):.4f}")
    called = scores >= 0.5
    print(f"  at 0.5: precision {(y[called] == 1).mean():.4f}  "
          f"recall {called[y == 1].mean():.4f}")

print("\nhow far the two sources disagree on the values")
same_rs = differ = both_absent = 0
gaps = []
for v in variants:
    dbs, rs_annovar, af_annovar = as_annovar[id(v)]
    rs_ensembl, af_ensembl = from_ensembl.get(
        (v.chrom, v.pos, v.ref, v.alt), (".", "0.0"))
    in_annovar = 0 if rs_annovar in (".", "", "NA") else 1
    in_ensembl = 0 if rs_ensembl in (".", "", "NA") else 1
    if in_annovar == in_ensembl:
        same_rs += 1
        if not in_annovar:
            both_absent += 1
    else:
        differ += 1
    best = max([float(x) for x in af_annovar if x not in (".", "", "NA")] or [0.0])
    gaps.append(abs(best - float(af_ensembl)))

gaps = np.array(gaps)
print(f"  dbSNP agrees on {same_rs} of {len(variants)} "
      f"({100 * same_rs / len(variants):.1f} percent), "
      f"of which {both_absent} are absent from both")
print(f"  gnomAD frequency: median difference {np.median(gaps):.6f}, "
      f"mean {gaps.mean():.6f}")
for cut in (0.001, 0.01, 0.1):
    print(f"    differs by more than {cut}: {(gaps > cut).sum()} variants")

print("\nhow far the scores move")
shift = np.abs(runs["ensembl_two"] - runs["annovar_two"])
print(f"  ensembl_two against annovar_two: median {np.median(shift):.4f}, "
      f"mean {shift.mean():.4f}, largest {shift.max():.4f}")
flipped = (runs["annovar_two"] >= 0.5) != (runs["ensembl_two"] >= 0.5)
print(f"  the call at 0.5 changes for {flipped.sum()} of {len(variants)} "
      f"({100 * flipped.sum() / len(variants):.2f} percent)")
