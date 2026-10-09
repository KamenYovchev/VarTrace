"""Loading the models and choosing between them.

The choice is made by what the input can support, not by anything the user has
to declare. What each model needs is read from its own description, and the
feature groups are listed here.
"""

import glob
import hashlib
import json
import os

import joblib

from vartrace import features

MODELS_DIR = os.path.join(os.path.dirname(__file__), "models")


class Model:
    def __init__(self, path):
        bundle = joblib.load(path)
        with open(path[: -len(".pkl")] + ".json") as fh:
            self.description = json.load(fh)
        self.file = os.path.basename(path)
        with open(path, "rb") as fh:
            self.fingerprint = hashlib.sha256(fh.read()).hexdigest()[:12]
        self.name = self.description.get("name", self.file)
        self.features = bundle["features"]
        self.estimator = bundle["model"]
        self.imputer = bundle.get("imputer")
        self.handles_missing = bundle.get("handles_missing", False)

    @property
    def groups(self):
        needed = set()
        if any(f in features.CROSSCELL for f in self.features):
            needed.add("cells")
        if "in_dbsnp" in self.features:
            needed.add("dbsnp")
        if "gnomad_af" in self.features:
            needed.add("gnomad")
        if "in_rediportal" in self.features:
            needed.add("redi")
        if "in_cosmic" in self.features:
            needed.add("cosmic")
        return needed

    def predict(self, variants):
        X = features.matrix(variants, self.features)
        if self.imputer is not None:
            X = self.imputer.transform(X)
        return self.estimator.predict_proba(X)[:, 1]

    def __repr__(self):
        return f"Model({self.name}, {len(self.features)} features)"


def load_all(directory=MODELS_DIR):
    return [Model(p) for p in sorted(glob.glob(os.path.join(directory, "*.pkl")))]


def choose(models, available):
    """Picks a model and says in words why that one.

    Full annotation or none at all is covered by a model trained for exactly
    that case. Anything in between goes to the model trained with feature
    groups withheld, which is the only one that can take a partial input.
    """
    dbs = available & features.DB_GROUPS
    cells = "cells" in available

    if dbs and dbs != features.DB_GROUPS:
        for m in models:
            if m.handles_missing:
                missing = sorted(features.DB_GROUPS - dbs)
                return m, (f"the input has {', '.join(sorted(dbs))} but not "
                           f"{', '.join(missing)}, so the model for partial "
                           f"annotation was used")
        raise RuntimeError("no model can take a partial annotation")

    wanted = set(dbs)
    if cells:
        wanted.add("cells")
    for m in models:
        if not m.handles_missing and m.groups == wanted:
            if not dbs and cells:
                why = "the input has no database annotation but several cells"
            elif not dbs:
                why = "the input has no database annotation and one cell"
            elif cells:
                why = "the input has all four databases and several cells"
            else:
                why = "the input has all four databases and one cell"
            return m, why + f", so {m.name} was used"
    raise RuntimeError(f"no model for {sorted(wanted)}")
