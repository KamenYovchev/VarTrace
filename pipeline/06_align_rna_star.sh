#!/bin/bash
# Build the STAR index if missing, then two-pass alignment of every RNA cell.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate variant-calling

IDX=${ROOT}/reference/grch38/star_index
FA=${ROOT}/reference/grch38/GRCh38.primary_assembly.genome.fa
GTF=${ROOT}/reference/grch38/gencode.v49.primary_assembly.annotation.gtf
TRIM=${ROOT}/results/trimmed
OUT=${ROOT}/results/aligned

# 1) Build STAR index only if it is not already present
if [ ! -f "$IDX/SAindex" ]; then
  echo ">>> Building STAR index ($(date))"
  mkdir -p "$IDX"
  STAR --runMode genomeGenerate --genomeDir "$IDX" \
    --genomeFastaFiles "$FA" --sjdbGTFfile "$GTF" \
    --sjdbOverhang 100 --runThreadN 8 || { echo ">>> index failed"; exit 1; }
fi

# 2) Align all cells (2-pass).
mkdir -p "$OUT"
find "$TRIM" -name "*_R1.trim.fastq.gz" | sort | while read r1; do
  r2="${r1/_R1.trim/_R2.trim}"
  base=$(basename "$r1" | sed 's/_R1.trim.fastq.gz//')
  grp=$(echo "$r1" | grep -oE 'molm13|hg002|hg001')
  # resume: skip cells that already have a finished BAM (added 2026-08-15 so new cells
  # can be aligned without overwriting the already-validated old ones)
  if [ -s "$OUT/$grp/$base/Aligned.sortedByCoord.out.bam" ]; then
    echo ">>> skip (done): $grp / $base"
    continue
  fi
  mkdir -p "$OUT/$grp/$base"
  echo ">>> STAR align: $grp / $base ($(date))"
  if STAR --runThreadN 8 --genomeDir "$IDX" \
       --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
       --twopassMode Basic --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$base SM:$base PL:ILLUMINA LB:$base \
       --quantMode GeneCounts \
       --outFileNamePrefix "$OUT/$grp/$base/"; then
    samtools index "$OUT/$grp/$base/Aligned.sortedByCoord.out.bam"
  else
    echo ">>> warn: failed $grp/$base ($(date))"
  fi
done
echo ">>> all done ($(date))"
