#!/bin/bash
# Realign and call the 10 HG002 cells close to BioSkryb's pipeline:
# Parabricks 4.4.0 fq2bam with their reference, then DeepVariant 1.8.0.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u
BASE=$HG002_DIR
REF=$BASE/ref_bioskryb/genome.fa
PB_SIF=$BASE/parabricks_4.4.0.sif
DV_SIF=$BASE/deepvariant_1.8.0-gpu.sif
SHARDS=${SHARDS:-36}   # make_examples threads per lane
mkdir -p "$BASE/bam" "$BASE/vcf" "$BASE/logs" "$BASE/tmp"

SAMPLES=(
  ResolveOMEv2.0-384-WGS-B4_L001
  ResolveOMEv2.0-384-WGS-D3_L001
  ResolveOMEv2.0-384-WGS-D6_L001
  ResolveOMEv2.0-384-WGS-F22_L001
  ResolveOMEv2.0-384-WGS-G22_L001
  ResolveOMEv2.0-384-WGS-H5_L001
  ResolveOMEv2.0-384-WGS-K3_L001
  ResolveOMEv2.0-384-WGS-M4_L001
  ResolveOMEv2.0-384-WGS-O11_L001
  ResolveOMEv2.0-384-WGS-O7_L001
)

run_sample () {
  local s=$1 gpu=$2
  local bam=$BASE/bam/${s}.bam
  local vcf=$BASE/vcf/${s}_deepvariant.vcf.gz
  if [ ! -s "$bam" ]; then
    echo ">>> [$s] fq2bam on GPU $gpu ($(date))"
    APPTAINERENV_CUDA_DEVICE_ORDER=PCI_BUS_ID \
    APPTAINERENV_CUDA_VISIBLE_DEVICES=$gpu \
    apptainer exec --nv --bind "$BASE":"$BASE" "$PB_SIF" \
      pbrun fq2bam \
      --ref "$REF" \
      --in-fq "$BASE/fastq/${s}_R1_001.fastq.gz" "$BASE/fastq/${s}_R2_001.fastq.gz" \
      --out-bam "$bam" \
      --num-gpus 1 \
      >> "$BASE/logs/${s}.fq2bam.log" 2>&1 || { echo "!!! [$s] fq2bam failed"; return 1; }
  fi
  if [ ! -s "$vcf" ]; then
    echo ">>> [$s] deepvariant 1.8.0 on GPU $gpu ($(date))"
    mkdir -p "$BASE/tmp/$s"
    APPTAINERENV_CUDA_DEVICE_ORDER=PCI_BUS_ID \
    APPTAINERENV_CUDA_VISIBLE_DEVICES=$gpu \
    apptainer exec --nv --bind "$BASE":"$BASE" "$DV_SIF" \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --ref="$REF" \
      --reads="$bam" \
      --output_vcf="$vcf" \
      --num_shards=$SHARDS \
      --intermediate_results_dir="$BASE/tmp/$s" \
      >> "$BASE/logs/${s}.deepvariant.log" 2>&1 || { echo "!!! [$s] deepvariant failed"; return 1; }
    rm -rf "$BASE/tmp/$s"
  fi
  echo ">>> [$s] done ($(date))"
}

# One lane per GPU from GPUS in config.sh, samples spread over the lanes
read -ra GPU_LIST <<< "$GPUS"
NL=${#GPU_LIST[@]}
for lane in $(seq 0 $((NL - 1))); do
  (
    gpu=${GPU_LIST[$lane]}
    i=$lane
    while [ $i -lt ${#SAMPLES[@]} ]; do
      run_sample "${SAMPLES[$i]}" "$gpu"
      i=$((i + NL))
    done
  ) &
done
wait
echo "All samples processed ($(date))"
echo "VCFs ready: $(ls "$BASE"/vcf/*.vcf.gz 2>/dev/null | wc -l) / 10"
