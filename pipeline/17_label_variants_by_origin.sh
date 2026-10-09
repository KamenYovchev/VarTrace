#!/bin/bash
# Label each RNA variant by origin using the DNA of the same cell.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

IN=${ROOT}/results/rna_annotated
DNADIR=$DATA
OUT=${ROOT}/results/labeled
mkdir -p "$OUT"
SUM=${ROOT}/results/label_summary.tsv
printf "group\tcell\tgermline\tsomatic\tnovel_DNA\tRNA_editing\tRNA_novel\tindeterminate\ttotal\n" > "$SUM"

for f in "$IN"/*.hg38_multianno.txt; do
  base=$(basename "$f" .hg38_multianno.txt)
  if [[ "$base" == *MOLM13* ]]; then
    grp=molm13
    dna="$DNADIR/molm13/dna_vcf/google-deepvariant_1.8.0_ptamodel/$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')_deepvariant.vcf.gz"
  elif [[ "$base" == ResolveOMEv2.1* ]]; then
    grp=hg001
    dna="$DNADIR/hg001/dna_vcf/google-deepvariant_1.8.0_ptamodel/$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')_deepvariant.vcf.gz"
  elif [[ "$base" == *NA12878* ]]; then
    # the 10 HG001 cells from ResolveOME v1.0, same donor, other file names
    grp=hg001
    dna="$DNADIR/hg001_v1/dna_vcf/$(echo "$base" | sed 's/^RNA/DNA/; s/NA12878-RNA-/NA12878-/; s/_S1_L001//')_deepvariant.vcf.gz"
  elif [[ "$base" == ResolveOMEv2.0-384* ]]; then
    # HG002 calls from steps 11 and 12, DeepVariant with the PTA model
    grp=hg002
    dna="$HG002_DIR/vcf_ptamodel/$(echo "$base" | sed 's/RNA-QC/WGS/')_deepvariant.vcf.gz"
  else
    echo "warn no mapping for $base"; continue
  fi
  if [ ! -f "$dna" ]; then echo "warn DNA not found: $dna"; continue; fi
  echo ">>> label $base ($grp) ($(date))"
  python "$(dirname "${BASH_SOURCE[0]}")"/17_label_variants_by_origin.py "$f" "$dna" "$grp" "$OUT/${base}.labeled.tsv" 2>> "$SUM"
done

echo ">>> labelling done ($(date))"
echo ""
echo "sanity check: total somatic per group (hg001/hg002 must be ~0)"
awk -F'\t' 'NR>1{s[$1]+=$4; g[$1]+=$3} END{for(k in s) printf "%s\tsomatic=%d\tgermline=%d\n", k, s[k], g[k]}' "$SUM"
echo ""
echo "full summary"
column -t "$SUM"
