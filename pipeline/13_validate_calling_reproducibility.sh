#!/bin/bash
# Check that our calling reproduces BioSkryb's VCF on one HG001 cell.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u

# Bucket keys come from the environment.
: "${BIOSKRYB_KEY:?set BIOSKRYB_KEY}"
: "${BIOSKRYB_SECRET:?set BIOSKRYB_SECRET}"
BASE=$HG002_DIR
VAL=$BASE/validation
REF=$BASE/ref_bioskryb/genome.fa
MODELDIR=$BASE/ptamodel/model
POPVCF="$BASE/ptamodel/population_vcfs/cohort-chr*.release_missing2ref.no_calls.vcf.gz"
CKPT=sc_af.model.ckpt
DV_SIF=$BASE/deepvariant_1.8.0-gpu.sif
SHARDS=36
BATCH=1024
CELL=DNA-OME-Manuscript-NA12878-SC1
REMOTE=":s3,provider=Wasabi,access_key_id=${BIOSKRYB_KEY},secret_access_key=${BIOSKRYB_SECRET},endpoint=${BIOSKRYB_ENDPOINT},region=${BIOSKRYB_REGION},force_path_style=true:"
SRC="$BIOSKRYB_DATA/ResolveOME_v1.0/manual/hg001/dna"
mkdir -p "$VAL/bam" "$VAL/theirs" "$VAL/ours" "$VAL/tmp" "$BASE/logs"

# 1.
if [ ! -s "$VAL/bam/${CELL}.bam" ]; then
  echo ">>> downloading their BAM ($(date))"
  "$BASE/rclone" copy "${REMOTE}${SRC}/bam/parabricks_4.4.0/${CELL}.bam"     "$VAL/bam/" -P
  "$BASE/rclone" copy "${REMOTE}${SRC}/bam/parabricks_4.4.0/${CELL}.bam.bai" "$VAL/bam/" -P
fi
if [ ! -s "$VAL/theirs/${CELL}_deepvariant.vcf.gz" ]; then
  echo ">>> downloading their VCF ($(date))"
  "$BASE/rclone" copy "${REMOTE}${SRC}/vcf/google-deepvariant_1.8.0_ptamodel/${CELL}_deepvariant.vcf.gz"     "$VAL/theirs/" -P
  "$BASE/rclone" copy "${REMOTE}${SRC}/vcf/google-deepvariant_1.8.0_ptamodel/${CELL}_deepvariant.vcf.gz.tbi" "$VAL/theirs/" -P
fi

# 2.
OURVCF=$VAL/ours/${CELL}_deepvariant.vcf.gz
LOG=$BASE/logs/validation_${CELL}.log
if [ ! -s "$OURVCF" ]; then
  echo ">>> make_examples ($(date))"
  APPTAINERENV_CUDA_DEVICE_ORDER=PCI_BUS_ID APPTAINERENV_CUDA_VISIBLE_DEVICES=${GPUS%% *} \
  apptainer exec --nv --bind "$BASE":"$BASE" "$DV_SIF" bash -c "
      seq 0 $((SHARDS-1)) | parallel -q --halt 2 --line-buffer \
        /opt/deepvariant/bin/make_examples \
          --mode calling --ref '$REF' --reads '$VAL/bam/${CELL}.bam' \
          --examples '$VAL/tmp/${CELL}.tfrecord@${SHARDS}.gz' \
          --checkpoint '$MODELDIR/$CKPT' \
          --population_vcfs='$POPVCF' \
          --task {}" >> "$LOG" 2>&1 || { echo "!!! make_examples FAILED"; exit 1; }

  echo ">>> call_variants ($(date))"
  APPTAINERENV_CUDA_DEVICE_ORDER=PCI_BUS_ID APPTAINERENV_CUDA_VISIBLE_DEVICES=${GPUS%% *} \
  apptainer exec --nv --bind "$BASE":"$BASE" "$DV_SIF" \
      /opt/deepvariant/bin/call_variants \
        --examples="$VAL/tmp/${CELL}.tfrecord@${SHARDS}.gz" \
        --outfile="$VAL/tmp/${CELL}_variants_output.tfrecord.gz" \
        --checkpoint="$MODELDIR/$CKPT" --batch_size=$BATCH >> "$LOG" 2>&1 \
    || { echo "!!! call_variants failed"; exit 1; }

  echo ">>> postprocess_variants ($(date))"
  apptainer exec --bind "$BASE":"$BASE" "$DV_SIF" \
      /opt/deepvariant/bin/postprocess_variants \
        --ref="$REF" --infile="$VAL/tmp/${CELL}_variants_output.tfrecord.gz" \
        --outfile="$OURVCF" --cpus=8 >> "$LOG" 2>&1 \
    || { echo "!!! postprocess failed"; exit 1; }
  rm -rf "$VAL/tmp"
fi

# 3.
ENV=$BASE/env
[ -x "$ENV/bin/bcftools" ] || mamba install -y -p "$ENV" -c bioconda -c conda-forge bcftools
BC=$ENV/bin/bcftools
THEIRS=$VAL/theirs/${CELL}_deepvariant.vcf.gz

echo ""
echo "Vcf header comparison"
diff <(zcat "$THEIRS"  | grep "^##" | grep -v "^##contig") \
     <(zcat "$OURVCF"  | grep "^##" | grep -v "^##contig") \
  && echo "headers identical (apart from contigs)" || echo "^ header differences above"

echo ""
echo "variant concordance"
# Real variants only (drop RefCall rows), then intersect on chr:pos:ref:alt
mkdir -p "$VAL/cmp"
$BC view -f PASS "$THEIRS" -Oz -o "$VAL/cmp/theirs.pass.vcf.gz" 2>/dev/null; $BC index -t -f "$VAL/cmp/theirs.pass.vcf.gz"
$BC view -f PASS "$OURVCF" -Oz -o "$VAL/cmp/ours.pass.vcf.gz"   2>/dev/null; $BC index -t -f "$VAL/cmp/ours.pass.vcf.gz"
$BC isec -p "$VAL/cmp/isec" "$VAL/cmp/theirs.pass.vcf.gz" "$VAL/cmp/ours.pass.vcf.gz" 2>/dev/null

only_theirs=$(grep -vc "^#" "$VAL/cmp/isec/0000.vcf" 2>/dev/null || echo 0)
only_ours=$(  grep -vc "^#" "$VAL/cmp/isec/0001.vcf" 2>/dev/null || echo 0)
shared=$(     grep -vc "^#" "$VAL/cmp/isec/0002.vcf" 2>/dev/null || echo 0)
t_tot=$((shared + only_theirs)); o_tot=$((shared + only_ours))
echo "PASS variants in theirs : $t_tot"
echo "PASS variants in ours   : $o_tot"
echo "shared                  : $shared"
echo "only in theirs          : $only_theirs"
echo "only in ours            : $only_ours"
if [ "$t_tot" -gt 0 ] && [ "$o_tot" -gt 0 ]; then
  awk -v s="$shared" -v t="$t_tot" -v o="$o_tot" 'BEGIN{
    printf "recall  (share of theirs we reproduce): %.4f\n", s/t;
    printf "precision (share of ours that is theirs): %.4f\n", s/o;
    printf "F1: %.4f\n", 2*(s/t)*(s/o)/((s/t)+(s/o));
  }'
fi
echo ""
echo "outputs in $val"
