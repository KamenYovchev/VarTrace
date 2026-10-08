#!/bin/bash
# Label each RNA variant by origin using the DNA of the same cell.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

IN=${ROOT}/results/rna_annotated
DNADIR=${ROOT}/singlecell_data
OUT=${ROOT}/results/step9_labeled
mkdir -p "$OUT"
SUM=${ROOT}/results/step9_summary.tsv
printf "group\tcell\tgermline\tsomatic\tnovel_DNA\tRNA_editing\tRNA_novel\tindeterminate\ttotal\n" > "$SUM"

for f in "$IN"/*.hg38_multianno.txt; do
  base=$(basename "$f" .hg38_multianno.txt)
  if [[ "$base" == *MOLM13* ]]; then
    grp=molm13
    dna="$dnadir/molm13/dna_vcf/google-deepvariant_1.8.0_ptamodel/$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')_deepvariant.vcf.gz"
  elif [[ "$base" == ResolveOMEv2.1* ]]; then
    grp=hg001
    dna="$dnadir/hg001/dna_vcf/google-deepvariant_1.8.0_ptamodel/$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')_deepvariant.vcf.gz"
  elif [[ "$base" == *NA12878* ]]; then
    # 10 additional HG001 cells (ResolveOME v1.0) added 2026-08-15 - same donor NA12878,
    # so same group; only the file-name pattern differs.
    grp=hg001
    dna="$dnadir/hg001_v1/dna_vcf/$(echo "$base" | sed 's/^RNA/DNA/; s/NA12878-RNA-/NA12878-/; s/_S1_L001//')_deepvariant.vcf.gz"
  elif [[ "$base" == ResolveOMEv2.0-384* ]]; then
    # hg002 re-processed 2026-08-15/16 with BioSkryb's own pipeline (Parabricks 4.4.0 +
    # their GATK-bundle reference + DeepVariant 1.8.0 with their PTA model), so it now
    # matches molm13/hg001 instead of the original Sentieon DNAscope calls.
    grp=hg002
    dna="$dnadir/hg002/dna_vcf_ptamodel/$(echo "$base" | sed 's/RNA-QC/WGS/')_deepvariant.vcf.gz"
  else
    echo "warn no mapping for $base"; continue
  fi
  if [ ! -f "$dna" ]; then echo "warn DNA not found: $dna"; continue; fi
  echo ">>> Step9 $base ($grp) ($(date))"
  python "$(dirname "${BASH_SOURCE[0]}")"/17_label_variants_by_origin.py "$f" "$dna" "$grp" "$OUT/${base}.labeled.tsv" 2>> "$SUM"
done

echo ">>> all STEP9 done ($(date))"
echo ""
echo "sanity check: total somatic per group (hg001/hg002 must be ~0)"
awk -F'\t' 'NR>1{s[$1]+=$4; g[$1]+=$3} END{for(k in s) printf "%s\tsomatic=%d\tgermline=%d\n", k, s[k], g[k]}' "$SUM"
echo ""
echo "full summary"
column -t "$SUM"
