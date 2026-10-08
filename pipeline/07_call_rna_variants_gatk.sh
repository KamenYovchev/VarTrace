#!/bin/bash
# Call RNA variants with GATK for every cell.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate variant-calling

REF=${ROOT}/reference/grch38/GRCh38.primary_assembly.genome.fa
ALIGN=${ROOT}/results/aligned
OUT=${ROOT}/results/vcf
TMP=${ROOT}/results/gatk_tmp
mkdir -p "$OUT" "$TMP"

# 0) One-time reference prep: .fai index + .dict dictionary (GATK requires both)
if [ ! -f "${REF}.fai" ]; then
  echo ">>> samtools faidx ($(date))"
  samtools faidx "$REF" || { echo ">>> faidx failed"; exit 1; }
fi
DICT="${REF%.fa}.dict"
if [ ! -f "$DICT" ]; then
  echo ">>> CreateSequenceDictionary ($(date))"
  gatk CreateSequenceDictionary -R "$REF" -O "$dict" || { echo ">>> dict failed"; exit 1; }
fi

# 1) Process every aligned cell.
find "$ALIGN" -name "Aligned.sortedByCoord.out.bam" | sort | while read bam; do
  base=$(basename "$(dirname "$bam")")
  grp=$(echo "$bam" | grep -oE 'molm13|hg002|hg001')
  cellout="$OUT/$grp/$base"
  mkdir -p "$cellout"

  # resume: skip if final VCF already present
  if [ -f "$cellout/${base}.filtered.vcf.gz" ]; then
    echo ">>> skip (done): $grp/$base"
    continue
  fi

  echo ">>> GATK start: $grp/$base ($(date))"

  # 1a) MarkDuplicates
  gatk MarkDuplicates \
    -I "$bam" \
    -O "$cellout/md.bam" \
    -M "$cellout/dup_metrics.txt" \
    --CREATE_INDEX true --TMP_DIR "$TMP" \
    || { echo ">>> warn MarkDuplicates failed $grp/$base"; continue; }

  # 1b) SplitNCigarReads (splits spliced reads at N cigar / intron gaps)
  gatk SplitNCigarReads \
    -R "$REF" \
    -I "$cellout/md.bam" \
    -O "$cellout/split.bam" \
    --tmp-dir "$TMP" \
    || { echo ">>> warn SplitNCigarReads failed $grp/$base"; continue; }

  # 1c) HaplotypeCaller (RNA-seq settings)
  gatk HaplotypeCaller \
    -R "$REF" \
    -I "$cellout/split.bam" \
    -O "$cellout/${base}.raw.vcf.gz" \
    --dont-use-soft-clipped-bases true \
    --standard-min-confidence-threshold-for-calling 20 \
    --native-pair-hmm-threads 4 \
    --tmp-dir "$TMP" \
    || { echo ">>> warn HaplotypeCaller failed $grp/$base"; continue; }

  # 1d) VariantFiltration (GATK RNA-seq hard filters: clustered SNPs, FS, QD)
  gatk VariantFiltration \
    -R "$REF" \
    -V "$cellout/${base}.raw.vcf.gz" \
    -O "$cellout/${base}.filtered.vcf.gz" \
    -window 35 -cluster 3 \
    --filter-name FS -filter "FS > 30.0" \
    --filter-name QD -filter "QD < 2.0" \
    --tmp-dir "$TMP" \
    || { echo ">>> warn VariantFiltration failed $grp/$base"; continue; }

  # cleanup large intermediates (keep raw + filtered VCF)
  rm -f "$cellout/md.bam" "$cellout/md.bai" "$cellout/md.bam.bai" \
        "$cellout/split.bam" "$cellout/split.bai" "$cellout/split.bam.bai"

  echo ">>> done $grp/$base ($(date))"
done
echo ">>> All gatk done ($(date))"
