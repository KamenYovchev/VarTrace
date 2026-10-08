#!/bin/bash
# Annotate the DNA VCFs with ANNOVAR.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

ANNOVAR=${ROOT}/annovar
HUMANDB=$ANNOVAR/humandb
IN=${ROOT}/results/dna_step6_mainchr
OUT=${ROOT}/results/dna_annotated
mkdir -p "$OUT"

PROTO="refGene,cytoBand,avsnp150,gnomad30_genome,gnomad211_exome,dbnsfp47a,clinvar_20240611,cosmic104,rediportal"
OP="g,r,f,f,f,f,f,f,f"

for f in "$IN"/*.mainchr.vcf.gz; do
  base=$(basename "$f" .mainchr.vcf.gz)
  echo ">>> Annotating DNA $base ($(date))"
  perl "$ANNOVAR/convert2annovar.pl" -format vcf4 -includeinfo "$f" 2>/dev/null \
    | sed 's/^chr//' > "$OUT/$base.avinput"
  perl "$ANNOVAR/table_annovar.pl" "$OUT/$base.avinput" "$HUMANDB/" \
    -buildver hg38 -out "$OUT/$base" -remove \
    -protocol "$PROTO" -operation "$OP" \
    -nastring . -otherinfo > "$out/$base.annovar.log" 2>&1
  if [ -f "$OUT/${base}.hg38_multianno.txt" ]; then
    n=$(( $(wc -l < "$OUT/${base}.hg38_multianno.txt") - 1 ))
    echo ">>> done DNA $base : $n annotated variants ($(date))"
  else
    echo ">>> Warn failed dna $base"
  fi
done
echo ">>> All dna annotation done ($(date))"
