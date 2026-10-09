#!/bin/bash
# Annotate the filtered RNA VCFs with ANNOVAR.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

ANNOVAR=${ROOT}/annovar
HUMANDB=$ANNOVAR/humandb
IN=${ROOT}/results/rna_qc_filtered
OUT=${ROOT}/results/rna_annotated
mkdir -p "$OUT"

PROTO="refGene,cytoBand,avsnp150,gnomad30_genome,gnomad211_exome,dbnsfp47a,clinvar_20240611,cosmic104,rediportal"
OP="g,r,f,f,f,f,f,f,f"

for f in "$IN"/*.qc.vcf.gz; do
  base=$(basename "$f" .qc.vcf.gz)
  # annotation takes hours, so skip cells already annotated
  if [ -s "$OUT/${base}.hg38_multianno.txt" ]; then
    echo ">>> skip (done): $base"; continue
  fi
  echo ">>> Annotating $base ($(date))"
  # VCF, then avinput; strip leading 'chr' from chrom (col 1) to match ANNOVAR dbs
  perl "$ANNOVAR/convert2annovar.pl" -format vcf4 -includeinfo "$f" 2>/dev/null \
    | sed 's/^chr//' > "$OUT/$base.avinput"
  # annotate with the nine databases in PROTO
  perl "$ANNOVAR/table_annovar.pl" "$OUT/$base.avinput" "$HUMANDB/" \
    -buildver hg38 -out "$OUT/$base" -remove \
    -protocol "$PROTO" -operation "$OP" \
    -nastring . -otherinfo > "$OUT/$base.annovar.log" 2>&1
  if [ -f "$OUT/${base}.hg38_multianno.txt" ]; then
    n=$(( $(wc -l < "$OUT/${base}.hg38_multianno.txt") - 1 ))
    echo ">>> done $base : $n annotated variants ($(date))"
  else
    echo ">>> warn failed $base (see $OUT/$base.annovar.log)"
  fi
done
echo ">>> All annotation done ($(date))"
