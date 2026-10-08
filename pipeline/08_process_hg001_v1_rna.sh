#!/bin/bash
# Trim, align and call the 10 HG001 v1.0 cells (steps 5 to 7 for these cells).
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh

RAW=${ROOT}/singlecell_data/hg001_v1/rna_fastq
TRIM=${ROOT}/results/trimmed/hg001
IDX=${ROOT}/reference/grch38/star_index
REF=${ROOT}/reference/grch38/GRCh38.primary_assembly.genome.fa
ALIGN=${ROOT}/results/aligned/hg001
OUT=${ROOT}/results/vcf/hg001
TMP=${ROOT}/results/gatk_tmp
mkdir -p "$TRIM" "$ALIGN" "$OUT" "$TMP"

# 1) fastp trimming
conda activate read-qc
for r1 in "$RAW"/*_R1_001.fastq.gz; do
  base=$(basename "$r1" _R1_001.fastq.gz)
  r2="${r1/_R1_001/_R2_001}"
  o1="$TRIM/${base}_R1.trim.fastq.gz"
  o2="$TRIM/${base}_R2.trim.fastq.gz"
  if [ -s "$o1" ] && [ -s "$o2" ]; then
    echo ">>> skip trim (done): $base"; continue
  fi
  echo ">>> fastp: $base ($(date))"
  fastp -i "$r1" -I "$r2" -o "$o1" -O "$o2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 25 \
        --thread 8 \
        --html "$TRIM/${base}.fastp.html" --json "$TRIM/${base}.fastp.json" \
    || { echo ">>> warn fastp failed $base"; continue; }
done

# 2) STAR 2-pass alignment
conda activate variant-calling
for r1 in "$TRIM"/RNA-OME-Manuscript-NA12878-*_R1.trim.fastq.gz; do
  [ -e "$r1" ] || continue
  r2="${r1/_R1.trim/_R2.trim}"
  base=$(basename "$r1" | sed 's/_R1.trim.fastq.gz//')
  celldir="$ALIGN/$base"
  if [ -s "$celldir/Aligned.sortedByCoord.out.bam" ]; then
    echo ">>> skip align (done): $base"; continue
  fi
  mkdir -p "$celldir"
  echo ">>> STAR align: hg001 / $base ($(date))"
  if STAR --runThreadN 8 --genomeDir "$IDX" \
       --readFilesIn "$r1" "$r2" --readFilesCommand zcat \
       --twopassMode Basic --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$base SM:$base PL:ILLUMINA LB:$base \
       --quantMode GeneCounts \
       --outFileNamePrefix "$celldir/"; then
    samtools index "$celldir/Aligned.sortedByCoord.out.bam"
  else
    echo ">>> warn: STAR failed $base ($(date))"
  fi
done

# 3) GATK RNA-seq variant calling
for bam in "$ALIGN"/RNA-OME-Manuscript-NA12878-*/Aligned.sortedByCoord.out.bam; do
  [ -e "$bam" ] || continue
  base=$(basename "$(dirname "$bam")")
  cellout="$OUT/$base"
  mkdir -p "$cellout"
  if [ -f "$cellout/${base}.filtered.vcf.gz" ]; then
    echo ">>> skip gatk (done): $base"; continue
  fi
  echo ">>> GATK start: hg001 / $base ($(date))"

  gatk MarkDuplicates -I "$bam" -O "$cellout/md.bam" \
    -M "$cellout/dup_metrics.txt" --CREATE_INDEX true --TMP_DIR "$TMP" \
    || { echo ">>> warn MarkDuplicates failed $base"; continue; }

  gatk SplitNCigarReads -R "$REF" -I "$cellout/md.bam" \
    -O "$cellout/split.bam" --tmp-dir "$TMP" \
    || { echo ">>> warn SplitNCigarReads failed $base"; continue; }

  gatk HaplotypeCaller -R "$REF" -I "$cellout/split.bam" \
    -O "$cellout/${base}.raw.vcf.gz" \
    --dont-use-soft-clipped-bases true \
    --standard-min-confidence-threshold-for-calling 20 \
    --native-pair-hmm-threads 4 --tmp-dir "$TMP" \
    || { echo ">>> warn HaplotypeCaller failed $base"; continue; }

  gatk VariantFiltration -R "$REF" -V "$cellout/${base}.raw.vcf.gz" \
    -O "$cellout/${base}.filtered.vcf.gz" \
    -window 35 -cluster 3 \
    --filter-name FS -filter "FS > 30.0" \
    --filter-name QD -filter "QD < 2.0" \
    --tmp-dir "$TMP" \
    || { echo ">>> warn VariantFiltration failed $base"; continue; }

  rm -f "$cellout/md.bam" "$cellout/md.bai" "$cellout/md.bam.bai" \
        "$cellout/split.bam" "$cellout/split.bai" "$cellout/split.bam.bai"
  echo ">>> done $base ($(date))"
done

echo "all done ($(date))"
echo "trimmed: $(ls "$TRIM"/RNA-ome-Manuscript-NA12878-*.trim.fastq.gz 2>/dev/null | wc -l) / 20"
echo "aligned: $(ls -d "$ALIGN"/RNA-ome-Manuscript-NA12878-*/ 2>/dev/null | wc -l) / 10"
echo "vcf    : $(ls "$OUT"/RNA-ome-Manuscript-NA12878-*/*.filtered.vcf.gz 2>/dev/null | wc -l) / 10"
