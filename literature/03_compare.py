#!/usr/bin/env python3
"""Compare the rule screening with the corpus used in the thesis.

A paper of the rebuilt pool counts as kept by the original screening if it is in
the corpus (matched by DOI, PMID or normalised title). Rule score 2 or 3 counts as kept by the rules.
"""
import re
import pandas as pd

pool = pd.read_csv("pool.tsv", sep="\t", dtype={"pmid": str})
corp = pd.read_csv("variant_origin_literature_corpus.csv", dtype={"pmid": str})
comp = pd.read_csv("benchmark_comparators_summary.csv", dtype={"PMID": str})


def norm(t):
    return re.sub(r"[^a-z0-9]", "", str(t).lower().replace("&amp;", "&"))[:80]


def ids(df, pmid="pmid"):
    out = []
    for _, r in df.iterrows():
        k = set()
        if isinstance(r.get("doi"), str):
            k.add("doi:" + r["doi"].lower())
        p = r.get(pmid)
        if isinstance(p, str) and p not in ("", "nan"):
            k.add("pmid:" + p.split(".")[0])
        k.add("ti:" + norm(r["title"]))
        out.append(k)
    return out


pool_ids = ids(pool)
index = {}
for i, k in enumerate(pool_ids):
    for x in k:
        index.setdefault(x, i)

corp["pool_row"] = [next((index[x] for x in k if x in index), None) for k in ids(corp)]
corp["in_pool"] = corp.pool_row.notna()
corp["rule_score"] = [pool.rule_score[int(i)] if pd.notna(i) else None for i in corp.pool_row]

lines = []
p = lines.append
p("Rebuilt pool: %d papers" % len(pool))
p("Rules keep (score 2 or 3): %d" % (pool.rule_score >= 2).sum())
for lab in ["directly relevant", "useful background"]:
    c = corp[corp.relevance_label == lab]
    inp = c[c.in_pool]
    kept = (inp.rule_score >= 2).sum()
    p("%s: %d in corpus, %d found in the rebuilt pool, rules keep %d of them (%.1f%%)"
      % (lab, len(c), len(inp), kept, 100 * kept / len(inp)))
    p("   rule scores: " + inp.rule_score.value_counts().sort_index().to_string().replace("\n", ", "))

# agreement on the whole pool: original keep = present in corpus.
# Every pool row is checked, so a preprint and its journal version both count.
corp_ids = set().union(*ids(corp))
pool["orig_keep"] = [bool(k & corp_ids) for k in pool_ids]
pool["rule_keep"] = pool.rule_score >= 2
ct = pd.crosstab(pool.orig_keep, pool.rule_keep)
p("\nCross table on the pool (rows original, columns rules):\n" + ct.to_string())
n = len(pool)
po = (pool.orig_keep == pool.rule_keep).mean()
pe = pool.orig_keep.mean() * pool.rule_keep.mean() + (1 - pool.orig_keep.mean()) * (1 - pool.rule_keep.mean())
p("Agreement %.3f, Cohen kappa %.3f" % (po, (po - pe) / (1 - pe)))
p("Precision of the rules against the original: %.3f" % (pool.orig_keep & pool.rule_keep).sum() / 1 if False else
  "Of the %d kept by rules, %d are in the corpus (%.1f%%)" % (pool.rule_keep.sum(), (pool.orig_keep & pool.rule_keep).sum(),
                                                             100 * (pool.orig_keep & pool.rule_keep).sum() / pool.rule_keep.sum()))

# the twelve comparator tools
p("\nComparator tools:")
pm = {str(x).split(".")[0] for x in pool.pmid.dropna()}
for _, r in comp.iterrows():
    pid = str(r["PMID"]).split(".")[0]
    rows = pool[pool.pmid.astype(str).str.split(".").str[0] == pid]
    if len(rows):
        p("  %-13s PMID %s  in pool, rule score %d" % (r["Tool"], pid, rows.rule_score.iloc[0]))
    else:
        p("  %-13s PMID %s  not in the rebuilt pool" % (r["Tool"], pid))

# directly relevant papers the rules would drop, for reading by hand
p("\nDirectly relevant papers the rules drop:")
miss = corp[(corp.relevance_label == "directly relevant") & corp.in_pool & (corp.rule_score < 2)]
for _, r in miss.iterrows():
    p("  [%d] %s (%s) %s" % (r.rule_score, r.section_title, int(r.year), r.title[:110]))
p("\nDirectly relevant papers not found in the rebuilt pool:")
for _, r in corp[(corp.relevance_label == "directly relevant") & ~corp.in_pool].iterrows():
    p("  %s (%s) %s" % (r.section_title, int(r.year), r.title[:110]))

open("results.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
