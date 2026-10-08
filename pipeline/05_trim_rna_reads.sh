#!/bin/bash
# Trim the single-cell RNA reads with fastp.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate read-qc

# data can be overridden to restrict the run to a subset, e.g.
DATA=${DATA:-${ROOT}/singlecell_data}
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
