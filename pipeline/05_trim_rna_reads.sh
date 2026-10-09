#!/bin/bash
# Trim the single-cell RNA reads with fastp and write the quality reports.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate read-qc

# all cells under $DATA, including the HG001 v1.0 cells; set DATA to a
# subfolder to run only some of them
OUT=${ROOT}/results/trimmed

find "$DATA" -path "*rna_fastq*" -name "*_R1_001.fastq.gz" | sort | while read r1; do
  r2="${r1/_R1_001/_R2_001}"
  base=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')
  grp=$(echo "$r1" | grep -oE 'molm13|hg002|hg001')
  mkdir -p "$OUT/$grp"
  echo ">>> $grp / $base"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$OUT/$grp/${base}_R1.trim.fastq.gz" \
    -O "$OUT/$grp/${base}_R2.trim.fastq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 25 \
    --thread 8 \
    --html "$OUT/$grp/${base}.fastp.html" \
    --json "$OUT/$grp/${base}.fastp.json"
done

# quality reports: FastQC on the raw reads, MultiQC over FastQC and fastp
QC=${ROOT}/results/qc
mkdir -p "$QC/fastqc_raw"
find "$DATA" -path "*rna_fastq*" -name "*.fastq.gz" | sort \
  | xargs fastqc -t 8 -o "$QC/fastqc_raw"
multiqc -f "$QC/fastqc_raw" -o "$QC" -n multiqc_raw
multiqc -f "$OUT" -o "$QC" -n multiqc_trimmed
