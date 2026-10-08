#!/bin/bash
# Recheck variants without a DNA record using the depth in the DNA BAM.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

: "${BIOSKRYB_KEY:?set BIOSKRYB_KEY}"
: "${BIOSKRYB_SECRET:?set BIOSKRYB_SECRET}"

SAMTOOLS=${ENVS}/variant-calling/bin/samtools
RCLONE=${ENVS}/rclone/bin/rclone
R=":s3,provider=Wasabi,access_key_id=${BIOSKRYB_KEY},secret_access_key=${BIOSKRYB_SECRET},endpoint=${BIOSKRYB_ENDPOINT},region=${BIOSKRYB_REGION},force_path_style=true:"
BASE="$BIOSKRYB_DATA"

# (re)start the http server at the illumina base
pkill -f "rclone serve http" 2>/dev/null
sleep 1
nohup "$RCLONE" serve http "${R}${BASE}/" --addr localhost:8080 --read-only \
  > /tmp/rclone_serve.log 2>&1 &
sleep 4

IN=${ROOT}/results/step9_labeled
OUT=${ROOT}/results/step9_labeled_bamfix
mkdir -p "$OUT"
SUM=${ROOT}/results/step9b_summary.tsv
printf "cell\tgermline\tsomatic\tnovel_DNA\tRNA_editing\tRNA_novel\tindeterminate\ttotal\n" > "$SUM"

# HG002 BAMs served over http by rclone, see below
HG002_BAM_URL="${HG002_BAM_URL:-http://localhost:8081}"

for f in "$IN"/*.labeled.tsv; do
  base=$(basename "$f" .labeled.tsv)
  if [[ "$base" == *MOLM13* ]]; then
    dna=$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')
    url="http://localhost:8080/ResolveOME_v1.0/manual/molm13/dna/bam/parabricks_4.4.0/${dna}.bam"
  elif [[ "$base" == ResolveOMEv2.1* ]]; then
    dna=$(echo "$base" | sed 's/RNA/DNA/; s/_S1_L001//')
    url="http://localhost:8080/ResolveOME_v2.1/manual/hg001/dna/bam/parabricks_4.4.0/${dna}.bam"
  elif [[ "$base" == *NA12878* ]]; then
    # 10 additional HG001 cells (ResolveOME v1.0) added 2026-08-15
    dna=$(echo "$base" | sed 's/^RNA/DNA/; s/NA12878-RNA-/NA12878-/; s/_S1_L001//')
    url="http://localhost:8080/ResolveOME_v1.0/manual/hg001/dna/bam/parabricks_4.4.0/${dna}.bam"
  elif [[ "$base" == ResolveOMEv2.0-384* ]]; then
    # HG002 BAMs from step 11.
    dna=$(echo "$base" | sed 's/RNA-QC/WGS/')
    url="${HG002_BAM_URL}/${dna}.bam"
  else
    url="NONE"
  fi
  echo ">>> step9b $base ($(date))"
  python "$(dirname "${BASH_SOURCE[0]}")"/18_rescue_undetermined_from_bam.py "$f" "$OUT/${base}.labeled.tsv" "$url" "$SAMTOOLS" 2>> "$SUM"
done

pkill -f "rclone serve http" 2>/dev/null
echo ">>> STEP9b done ($(date))"
echo ""
echo "summary after BAM rescue"
column -t "$SUM"
echo ""
echo "indeterminate: before vs after (RNA-originating gained)"
echo "before (step9):"; awk -F'\t' 'NR>1{i+=$8; e+=$5; nv+=$6} end{printf "  indeterminate=%d  RNA_editing=%d  RNA_novel=%d\n", i,e,nv}' ${ROOT}/results/step9_summary.tsv
echo "after  (step9b):"; awk -F'\t' 'NR>1{i+=$7; e+=$5; nv+=$6} end{printf "  indeterminate=%d  RNA_editing=%d  RNA_novel=%d\n", i,e,nv}' "$sum"
