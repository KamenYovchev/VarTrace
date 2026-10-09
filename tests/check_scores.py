"""End to end check that the tool's scores match the model run directly.

check_features.py showed that the features come out identical to the training
table. This goes one step further: it runs the whole tool on the labelled
tables and compares the scores with what the same model gives when the training
table is fed to it straight. That covers the parts features.py does not, namely
the order of the feature matrix, the imputer, and the choice of model.

Run per cell line, because that is how the cross-cell features were computed
during training.

    PYTHONPATH=. python tests/check_scores.py <labelled_dir> <model_dataset_sig.csv> [models_dir]
"""

import glob
import os
import sys

import numpy as np
import pandas as pd

from vartrace import features, io, model as models

labeled_dir, training_csv = sys.argv[1], sys.argv[2]
models_dir = sys.argv[3] if len(sys.argv) > 3 else models.MODELS_DIR


def group_of(cell):
    if "MOLM13" in cell:
        return "molm13"
    if cell.startswith("ResolveOMEv2.1") or "NA12878" in cell:
        return "hg001"
    if cell.startswith("ResolveOMEv2.0-384"):
        return "hg002"
    return "other"


expected = pd.read_csv(training_csv).dropna(subset=["label"]).reset_index(drop=True)
wanted = set(zip(expected["cell"], expected["chr"],
                 expected["pos"].astype(str), expected["ref"], expected["alt"]))

loaded = models.load_all(models_dir)

by_line = {}
for path in sorted(glob.glob(os.path.join(labeled_dir, "*.labeled.tsv"))):
    got, _, _ = io.read_annovar(path)
    got = [v for v in got if (v.cell, v.chrom, v.pos, v.ref, v.alt) in wanted]
    if got:
        by_line.setdefault(group_of(got[0].cell), []).append(got)

worst = 0.0
for line, cells in sorted(by_line.items()):
    variants = [v for cell in cells for v in cell]
    for v in variants:
        features.variant_features(v)
    features.add_cross_cell(variants, len(cells))

    model, why = models.choose(loaded, features.available(variants))
    ours = model.predict(variants)

    # the same model, fed the training table directly
    rows = expected[expected["cell"].isin({v.cell for v in variants})]
    lookup = {}
    for r in rows.itertuples(index=False):
        lookup[(r.cell, r.chr, str(r.pos), r.ref, r.alt)] = r
    X = np.array([[getattr(lookup[(v.cell, v.chrom, v.pos, v.ref, v.alt)], name)
                   for name in model.features] for v in variants], dtype=float)
    if model.imputer is not None:
        X = model.imputer.transform(X)
    theirs = model.estimator.predict_proba(X)[:, 1]

    gap = float(np.abs(ours - theirs).max())
    worst = max(worst, gap)
    print(f"{line}: {len(variants)} variants in {len(cells)} cells, "
          f"{model.file}, largest difference {gap:.2e}")

print()
if worst == 0.0:
    print("the tool gives exactly the same scores")
else:
    print(f"largest difference anywhere: {worst:.2e}")
sys.exit(0 if worst == 0.0 else 1)
