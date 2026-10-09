#!/bin/bash
# Quality filter on the RNA variants.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate variant-calling

IN=${ROOT}/results/rna_mainchr
OUT=${ROOT}/results/rna_qc_filtered
mkdir -p "$OUT"
SUM=${ROOT}/results/rna_qc_filter_summary.tsv
printf "cell\tbefore\tafter\n" > "$SUM"

# format fields: DP, AD[0:1]=alt depth, GQ ; INFO fields: Qd, mq, fs, sor
EXPR='FMT/DP>=3 && FMT/AD[0:1]>=3 && FMT/GQ>=30 && INFO/QD>=2 && INFO/MQ>=50 && INFO/FS<=30 && INFO/SOR<=3'

for f in "$IN"/*.vcf.gz; do
  base=$(basename "$f" .mainchr.vcf.gz)
  out="$OUT/${base}.qc.vcf.gz"
  # skip cells already filtered
  if [ -s "$out" ]; then echo "skip (done): $base"; continue; fi
  bcftools view -f PASS "$f" -Ou 2>/dev/null \
    | bcftools view -i "$EXPR" -Oz -o "$out" 2>/dev/null
  bcftools index -t "$out" 2>/dev/null
  b=$(bcftools view -H "$f"   | wc -l)
  a=$(bcftools view -H "$out" | wc -l)
  printf "%s\t%s\t%s\n" "$base" "$b" "$a" >> "$SUM"
  echo "$base: $b -> $a"
done
echo ""
echo ">>> done. Filtered VCFs in $OUT. Summary: $SUM"
column -t "$SUM"
