"""Ask Ensembl about every labelled variant and keep the answers in a file.

The second half of this check, compare_annotation_sources.py, needs the models
and so runs under ml-cpu. That environment has no requests, and Ensembl should
not be asked twice anyway, so the fetching is done here on its own and written
out once.

Run with an environment that has requests:

    PYTHONPATH=. python tests/fetch_ensembl.py <model_dataset_sig.csv> <ensembl_answers.tsv>

The questions are put through vartrace.annotate itself rather than through a
copy of it, so what is measured afterwards is what the tool really does.
"""

import csv
import sys

from vartrace import annotate
from vartrace.variant import Variant

training_csv, out_path = sys.argv[1], sys.argv[2]

seen = {}
with open(training_csv) as fh:
    for row in csv.DictReader(fh):
        if not row.get("label"):
            continue
        key = (row["chr"], row["pos"], row["ref"], row["alt"])
        seen.setdefault(key, Variant(cell="", chrom=key[0], pos=key[1],
                                     ref=key[2], alt=key[3]))

variants = list(seen.values())
print(f"{len(variants)} unique labelled variants to ask about")


def say(done, total):
    print(f"  {done} of {total}", flush=True)


report = annotate.annotate(variants, progress=say)
print(f"Ensembl release {report['release']} had a record for "
      f"{report['with_a_record']} of {report['queried']}")

with open(out_path, "w") as fh:
    fh.write("chrom\tpos\tref\talt\trs\tgnomad_af\n")
    for v in variants:
        af = v.gnomad_af[0] if v.gnomad_af else "0.0"
        fh.write(f"{v.chrom}\t{v.pos}\t{v.ref}\t{v.alt}\t{v.dbsnp}\t{af}\n")

print(f"written to {out_path}")
