#!/usr/bin/env python3
"""Rebuild the candidate pool and screen it with fixed keyword rules.

Steps 1 and 2 repeat the original pipeline: deduplicate by DOI, PMID or title,
compute the same keyword topic score and keep records with score >= 3 and an
abstract longer than 100 characters.

Step 3 is new. It replaces the screening of the abstracts with rules that were
written before the comparison with the corpus and were not changed afterwards.
The rules see the same text the original screening saw: the title and the
first 900 characters of the abstract.

Output: pool.tsv with the rule score for every candidate.
"""
import json, re
import pandas as pd

raw = json.load(open("raw.json"))


def rows(recs, theme):
    out = []
    for r in recs:
        out.append({"theme": theme, "pmid": r.get("pmid"), "doi": r.get("doi"),
                    "title": (r.get("title") or "").strip().rstrip("."),
                    "year": r.get("pubYear"), "cited": r.get("citedByCount", 0),
                    "abstract": r.get("abstractText") or ""})
    return out


def key(r):
    if isinstance(r["doi"], str) and r["doi"]:
        return "doi:" + r["doi"].lower()
    if isinstance(r["pmid"], str) and r["pmid"]:
        return "pmid:" + r["pmid"]
    return "ti:" + re.sub(r"[^a-z0-9]", "", (r["title"] or "").lower())[:80]


def frame(recs):
    df = pd.DataFrame(recs)
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df["cited"] = pd.to_numeric(df["cited"], errors="coerce").fillna(0).astype(int)
    df["key"] = df.apply(key, axis=1)
    return df


# the original keyword topic score, unchanged
core_terms = {
 r"single[- ]cell rna": 3, r"scrna-?seq": 3, r"rna editing|a-to-i|adar": 3,
 r"somatic (?:mutation|variant)": 2, r"variant calling|genotyp": 2,
 r"mutational signature": 2, r"allelic dropout|allele drop": 3,
 r"transcription(?:al)? error|transcript fidelity|rna polymerase fidelity": 3,
 r"paired .{0,20}(?:dna|genome).{0,20}(?:rna|transcriptome)|g&t-seq|dr-seq|scwgs": 3,
 r"neoantigen|neoepitope": 2, r"strand bias|transcription-coupled repair|strand asymmetr": 2,
 r"gradient boosting|random forest|xgboost|machine learning|deep learning": 1,
 r"label noise|noisy label|positive-unlabeled|calibrat": 1,
 r"rna-dna difference|rdd|rna variant": 3, r"false positive|artifact|artefact": 1,
 r"acute myeloid|aml|leukemi": 1,
}


def topic_score(df):
    txt = (df["title"].fillna("") + " || " + df["abstract"].fillna("")).str.lower()
    s = pd.Series(0.0, index=df.index)
    for pat, w in core_terms.items():
        s += txt.str.contains(pat, regex=True) * w
    return s


# step 1: twelve themes plus named tools
main = []
for theme, recs in raw.items():
    if theme not in ("targeted", "must"):
        main += rows(recs, theme)
df = frame(main)
ded = df.sort_values("cited", ascending=False).drop_duplicates("key")
ded = ded[ded["year"].notna()].copy()
ded["topic_score"] = topic_score(ded)
pool1 = ded[(ded.topic_score >= 3) & (ded.abstract.str.len() > 100)]

# step 2: targeted searches, then the landmark searches, new records only
d2 = frame(rows(raw["targeted"], "targeted"))
d2 = d2[~d2.key.isin(ded.key)].drop_duplicates("key")
d2 = d2[d2.year.notna() & (d2.abstract.str.len() > 100)].copy()
d2["topic_score"] = topic_score(d2)
pool2 = d2[d2.topic_score >= 3]

d3 = frame(rows(raw["must"], "must"))
d3 = d3[~d3.key.isin(ded.key) & ~d3.key.isin(d2.key)].drop_duplicates("key")
d3 = d3[d3.abstract.str.len() > 100].copy()
d3["topic_score"] = topic_score(d3)
pool3 = d3[d3.topic_score >= 3]

print("records after deduplication:", len(ded), "+", len(d2), "+", len(d3), "=", len(ded) + len(d2) + len(d3))
print("pool:", len(pool1), "+", len(pool2), "+", len(pool3), "=", len(pool1) + len(pool2) + len(pool3))
pool = pd.concat([pool1, pool2, pool3], ignore_index=True)

# step 3: screening rules
# A paper counts as on topic if it matches at least one of the nine topics of the review.
SC = r"single[- ]cell|scrna|scdna|scwgs"
VAR = r"variant|mutation|snv|genotyp"
topics = {
 "A scRNA variant calling":       [SC, r"scrna|single[- ]cell rna|transcriptom", VAR],
 "B RNA editing":                 [r"rna editing|a-to-i|adar"],
 "C RNA-DNA differences":         [r"rna-dna difference|transcription(?:al)? error|transcriptional mutagenesis|polymerase fidelity|transcription fidelity"],
 "D paired DNA and RNA":          [SC, r"genome and transcriptome|dna and rna|g&t-seq|dr-seq|paired|simultaneous|joint"],
 "E scWGS artefacts and ADO":     [r"allelic dropout|allele dropout|amplification (?:artifact|artefact|bias|error)|whole[- ]genome amplification|single[- ]cell whole[- ]genome|scwgs|primary template"],
 "F signatures and strand":       [r"mutational signature|sigprofiler|transcription-coupled repair|strand bias|strand asymmetr|transcriptional strand"],
 "G ML for variants":             [r"machine learning|random forest|gradient boost|xgboost|deep learning|neural network|classifier|label noise|calibrat", VAR],
 "H expressed mutations":         [r"neoantigen|neoepitope|expressed (?:somatic )?mutation", r"rna|express|transcript"],
 "I single-cell clonal genotyping": [SC, r"leukemi|aml|clonal hematopoiesis|myeloid|lineage tracing|mitochondrial", VAR],
}
# The paper contributes a method, a tool, a resource or a benchmark.
METHOD = (r"\bwe (?:develop|developed|present|introduce|propose|describe|built|benchmark)"
          r"|\b(?:tool|software|algorithm|pipeline|framework|database|benchmark|resource)\b")
# Typical clinical or prognostic papers that the review rubric names as not relevant.
EXCLUDE = (r"prognos|nomogram|risk (?:model|score|signature)|gene signature|immune infiltrat"
           r"|overall survival|immunotherapy response")
# The paper is about the question of the thesis: a variant seen in RNA compared with DNA.
CORE = r"(?:rna|transcriptom).{0,120}(?:dna|genom)|(?:dna|genom).{0,120}(?:rna|transcriptom)"


def rule_score(title, abstract):
    t = (title + " || " + abstract[:900]).lower()
    hit = [name for name, pats in topics.items() if all(re.search(p, t) for p in pats)]
    if not hit:
        return 0, ""
    if re.search(EXCLUDE, t) and not re.search(VAR, title.lower()):
        return 0, "excluded"
    s = 1
    if re.search(METHOD, t):
        s += 1
    if len(hit) >= 2 or (re.search(CORE, t) and re.search(VAR, t)):
        s += 1
    return s, "; ".join(hit)


scored = pool.apply(lambda r: rule_score(r["title"], r["abstract"]), axis=1)
pool["rule_score"] = [s for s, _ in scored]
pool["rule_topics"] = [h for _, h in scored]
pool[["key", "pmid", "doi", "year", "cited", "topic_score", "rule_score", "rule_topics", "title"]] \
    .to_csv("pool.tsv", sep="\t", index=False)
print(pool.rule_score.value_counts().sort_index().to_string())
