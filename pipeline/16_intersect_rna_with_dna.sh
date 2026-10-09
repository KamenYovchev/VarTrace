#!/bin/bash
# Join each cell's RNA variants with its DNA VCF by chr, pos, ref and alt.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

IN=${ROOT}/results/rna_annotated
DNADIR=$DATA
OUT=${ROOT}/results/rna_dna_intersect
mkdir -p "$OUT"

for f in "$IN"/*.hg38_multianno.txt; do
  base=$(basename "$f" .hg38_multianno.txt)

  # map RNA cell name, then corresponding DNA VCF (naming differs per group)
  if [[ "$base" == *MOLM13* ]]; then
    dbase="$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')_deepvariant"
    dna="$DNADIR/molm13/dna_vcf/google-deepvariant_1.8.0_ptamodel/${dbase}.vcf.gz"
  elif [[ "$base" == ResolveOMEv2.1* ]]; then
    dbase="$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')_deepvariant"
    dna="$DNADIR/hg001/dna_vcf/google-deepvariant_1.8.0_ptamodel/${dbase}.vcf.gz"
  elif [[ "$base" == *NA12878* ]]; then
    dbase="$(echo "$base" | sed 's/^RNA/DNA/; s/NA12878-RNA-/NA12878-/; s/_S1_L001//')_deepvariant"
    dna="$DNADIR/hg001_v1/dna_vcf/${dbase}.vcf.gz"
  elif [[ "$base" == ResolveOMEv2.0-384* ]]; then
    dbase="$(echo "$base" | sed 's/RNA-QC/WGS/')_deepvariant"
    dna="$HG002_DIR/vcf_ptamodel/${dbase}.vcf.gz"
  else
    echo ">>> warn no DNA mapping for $base"; continue
  fi

  if [ ! -f "$dna" ]; then echo ">>> warn DNA not found for $base: $dna"; continue; fi

  echo ">>> Intersect $base  x  $(basename "$dna") ($(date))"
  python "$(dirname "${BASH_SOURCE[0]}")"/16_intersect_rna_with_dna.py "$f" "$dna" "$OUT/${base}.intersect.tsv"
done
echo ">>> All intersect done ($(date))"
