#!/bin/bash
# Call the HG002 BAMs again with BioSkryb's PTA model for DeepVariant,
# the same model used for the MOLM13 and HG001 VCFs.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u
BASE=$HG002_DIR
REF=$BASE/ref_bioskryb/genome.fa
MODELDIR=$BASE/ptamodel/model
POPVCF="$BASE/ptamodel/population_vcfs/cohort-chr*.release_missing2ref.no_calls.vcf.gz"
CKPT=sc_af.model.ckpt
DV_SIF=$BASE/deepvariant_1.8.0-gpu.sif
SHARDS=36
BATCH=1024   # matches BioSkryb's nextflow.config: call_variants_batch_size = 1024
OUT=$BASE/vcf_ptamodel
mkdir -p "$OUT" "$BASE/logs" "$BASE/tmp_pta"

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
  local vcf=$OUT/${s}_deepvariant.vcf.gz
  local wd=$BASE/tmp_pta/$s
  local log=$BASE/logs/${s}.ptamodel.log
  [ -s "$vcf" ] && { echo ">>> skip (done): $s"; return 0; }
  [ -s "$bam" ] || { echo "!!! [$s] BAM missing"; return 1; }
  mkdir -p "$wd"
  echo ">>> [$s] make_examples ($(date))"
  APPTAINERENV_CUDA_DEVICE_ORDER=PCI_BUS_ID APPTAINERENV_CUDA_VISIBLE_DEVICES=$gpu \
  apptainer exec --nv --bind "$BASE":"$BASE" "$DV_SIF" bash -c "
      seq 0 $((SHARDS-1)) | parallel -q --halt 2 --line-buffer \
        /opt/deepvariant/bin/make_examples \
          --mode calling \
          --ref '$REF' \
          --reads '$bam' \
          --examples '$wd/${s}.tfrecord@${SHARDS}.gz' \
          --checkpoint '$MODELDIR/$CKPT' \
          --population_vcfs='$POPVCF' \
          --task {}" >> "$log" 2>&1 \
    || { echo "!!! [$s] make_examples failed"; return 1; }

  echo ">>> [$s] call_variants on GPU $gpu ($(date))"
  APPTAINERENV_CUDA_DEVICE_ORDER=PCI_BUS_ID APPTAINERENV_CUDA_VISIBLE_DEVICES=$gpu \
  apptainer exec --nv --bind "$BASE":"$BASE" "$DV_SIF" \
      /opt/deepvariant/bin/call_variants \
        --examples="$wd/${s}.tfrecord@${SHARDS}.gz" \
        --outfile="$wd/${s}_variants_output.tfrecord.gz" \
        --checkpoint="$MODELDIR/$CKPT" \
        --batch_size=$BATCH >> "$log" 2>&1 \
    || { echo "!!! [$s] call_variants failed"; return 1; }

  echo ">>> [$s] postprocess_variants ($(date))"
  apptainer exec --bind "$BASE":"$BASE" "$DV_SIF" \
      /opt/deepvariant/bin/postprocess_variants \
        --ref="$REF" \
        --infile="$wd/${s}_variants_output.tfrecord.gz" \
        --outfile="$vcf" \
        --cpus=8 >> "$log" 2>&1 \
    || { echo "!!! [$s] postprocess failed"; return 1; }

  rm -rf "$wd"
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
echo "PTA-Model calling done ($(date))"
echo "VCFs: $(ls "$OUT"/*.vcf.gz 2>/dev/null | wc -l) / 10  -> $OUT"
