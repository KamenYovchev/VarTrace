#!/usr/bin/env python3
"""Match our MOLM13 RNA variants against the MOLM13 DepMap somatic truth set.

Usage: 38_match_depmap_somatic.py <molm13_depmap_somatic.tsv> '<labeled_glob>'
"""
import sys
import glob

depmap_file, labeled_glob = sys.argv[1], sys.argv[2]

somatic = set()
with open(depmap_file) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) < 4:
            continue
        somatic.add((c[0], c[1], c[2], c[3]))

files = sorted(glob.glob(labeled_glob))
total = 0
matched = 0
by_label = {}
hit_positions = set()

for f in files:
    with open(f) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ic, ip, ir, ia = (header.index(x) for x in
                          ("Otherinfo1", "Otherinfo2", "Otherinfo4", "Otherinfo5"))
        il = header.index("label")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            total += 1
            chrom = c[ic] if c[ic].startswith("chr") else "chr" + c[ic]
            key = (chrom, c[ip], c[ir], c[ia])
            if key in somatic:
                matched += 1
                by_label[c[il]] = by_label.get(c[il], 0) + 1
                hit_positions.add(key)

print("MOLM13 cells scanned:            ", len(files))
print("total variant rows (all cells):  ", total)
print("rows matching DepMap somatic:    ", matched)
print("unique DepMap positions hit:     ", len(hit_positions), "of", len(somatic))
print("current labels of matched rows:  ", by_label)
