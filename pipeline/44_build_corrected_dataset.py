#!/usr/bin/env python3
"""Write a copy of the training table with the contradicted labels corrected.

Usage: 44_build_corrected_dataset.py <dataset_sig.csv> <step 40 outdir> <step 43 outdir> <out.csv>
"""
import os
import sys

import pandas as pd

csv_in, giab_dir, molm_dir, csv_out = sys.argv[1:5]
KEY = ["cell", "chr", "pos", "ref", "alt"]

df = pd.read_csv(csv_in)
g = pd.concat([pd.read_csv(os.path.join(giab_dir, "%s_rows.tsv" % l), sep="\t")
               for l in ["hg001", "hg002"]])
g["wrong"] = (g.label == "RNA_only") & g.in_region & g.giab_match
m = pd.read_csv(os.path.join(molm_dir, "molm13_rows.tsv"), sep="\t")
flags = pd.concat([g[KEY + ["wrong"]], m[KEY + ["wrong"]]])
assert not flags.duplicated(KEY).any()

n = len(df)
df = df.merge(flags, on=KEY, how="left")
df["wrong"] = df["wrong"].fillna(False).astype(bool)
assert len(df) == n

print("corrected:", int(df.wrong.sum()), df[df.wrong].group.value_counts().to_dict())
df.loc[df.wrong, "label"] = "DNA_origin"
df.drop(columns="wrong").to_csv(csv_out, index=False)
print(df.label.value_counts().to_dict(), "RNA_only share %.4f" % (df.label == "RNA_only").mean())
