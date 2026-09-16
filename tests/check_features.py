"""Check that features.py reproduces the training table exactly.

The models were trained on model_dataset_sig.csv, built by build_dataset.py.
If the tool computes a feature even slightly differently, the models receive
something other than what they learned on, and nothing else in the tool means
anything. So every value of every feature is compared, for every row.

Both files live on hpc00, so this runs there:

    PYTHONPATH=/data/vartrace /data/conda/envs/ml-cpu/bin/python \\
        /data/vartrace/tests/check_features.py \\
        /data/results/step9_labeled_bamfix /data/results/model_dataset_sig.csv

Two things are mirrored from build_dataset.py on purpose. Cross-cell features
were computed per cell line, not over all 45 cells, so the cells are grouped the
same way here. And they counted only labelled variants, so rows that are not in
the training table are dropped before the cross-cell step.
"""

import glob
import os
import sys

import numpy as np
import pandas as pd

from vartrace import features, io

labeled_dir, training_csv = sys.argv[1], sys.argv[2]


def group_of(cell):
    if "MOLM13" in cell:
        return "molm13"
    if cell.startswith("ResolveOMEv2.1") or "NA12878" in cell:
        return "hg001"
    if cell.startswith("ResolveOMEv2.0-384"):
        return "hg002"
    return "other"


expected = pd.read_csv(training_csv)
expected = expected.dropna(subset=["label"]).reset_index(drop=True)
print(f"training table: {len(expected)} rows, {len(expected.columns)} columns")

wanted = set(zip(expected["cell"], expected["chr"],
                 expected["pos"].astype(str), expected["ref"], expected["alt"]))

variants = []
for path in sorted(glob.glob(os.path.join(labeled_dir, "*.labeled.tsv"))):
    got, _ = io.read_annovar(path)
    variants.extend(v for v in got
                    if (v.cell, v.chrom, v.pos, v.ref, v.alt) in wanted)
print(f"read back: {len(variants)} variants")

by_line = {}
for v in variants:
    by_line.setdefault(group_of(v.cell), []).append(v)

for line, group in by_line.items():
    for v in group:
        features.variant_features(v)
    n_cells = len({v.cell for v in group})
    features.add_cross_cell(group, n_cells)
    print(f"  {line}: {len(group)} variants in {n_cells} cells")

mine = {}
for v in variants:
    mine[(v.cell, v.chrom, v.pos, v.ref, v.alt)] = v.features

missing = 0
diffs = {}
for row in expected.itertuples(index=False):
    key = (row.cell, row.chr, str(row.pos), row.ref, row.alt)
    f = mine.get(key)
    if f is None:
        missing += 1
        continue
    for name in features.ALL_FEATURES:
        theirs = getattr(row, name)
        ours = f[name]
        theirs_missing = theirs is None or (isinstance(theirs, float) and np.isnan(theirs))
        if ours is None or theirs_missing:
            same = (ours is None) and theirs_missing
        else:
            same = float(ours) == float(theirs)
        if not same:
            d = diffs.setdefault(name, [])
            if len(d) < 3:
                d.append((key, theirs, ours))
            diffs[name] = d

print()
if missing:
    print(f"NOT FOUND in the input: {missing} rows")
if diffs:
    print("MISMATCHES")
    for name, examples in diffs.items():
        print(f"  {name}")
        for key, theirs, ours in examples:
            print(f"    {key[1]}:{key[2]} {key[3]}>{key[4]} in {key[0]}")
            print(f"      training {theirs!r}   tool {ours!r}")
else:
    print("every feature of every row matches")

sys.exit(1 if (diffs or missing) else 0)
