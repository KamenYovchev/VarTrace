Literature selection for the thesis

The corpus of 332 papers (variant_origin_literature_corpus.csv) was built in
three steps: a Europe PMC search on 12 topics, a keyword filter, and a 0 to 3
relevance score of each abstract against a description of the project, keeping
scores 3 and 2. The 12 closest methods were extracted from their full texts into
benchmark_comparators_summary.csv.

The three scripts repeat the selection with plain keyword rules:

  python 01_fetch.py          runs the same 12 Europe PMC queries, writes raw.json
  python 02_screen_rules.py   removes duplicates, scores keywords and abstracts, writes pool.tsv
  python 03_compare.py        compares the result with the corpus, writes results.txt

Run them in this folder, in this order. 01_fetch.py needs internet access.
02 and 03 need pandas.
