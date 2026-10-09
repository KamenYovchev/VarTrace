#!/bin/bash
# Align and call the MOLM13 bulk sample from step 41.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u
BASE=$RESULTS/molm13_bulk
REF=${REF:-$HG002_DIR/ref_bioskryb/genome.fa}
SIF=${SIF:-$HG002_DIR}
PB="apptainer exec --nv --bind $ROOT:$ROOT $SIF/parabricks_4.4.0.sif pbrun"
DV="apptainer exec --nv --bind $ROOT:$ROOT $SIF/deepvariant_1.8.0-gpu.sif"
BWA="apptainer exec --bind $ROOT:$ROOT $SIF/bwa.sif bwa"            # bwa:0.7.17--he4a0461_11
SAM="apptainer exec --bind $ROOT:$ROOT $SIF/samtools.sif samtools"  # samtools:1.21--h50ea8bc_0
S=SRR17524983
RG='@RG\tID:SRR17524983.1.1\tLB:lib1\tPL:bar\tSM:sample\tPU:SRR17524983.1.1'
THREADS=${THREADS:-80}
mkdir -p "$BASE/bam" "$BASE/vcf" "$BASE/logs" "$BASE/tmp"
cd "$BASE"

if [ ! -s bam/$S.aligned.bam ]; then
  set -o pipefail
  $BWA mem -t "$THREADS" -K 10000000 -R "$RG" "$REF" fastq/${S}_1.fastq.gz fastq/${S}_2.fastq.gz 2> logs/$S.bwa.log \
    | $SAM view -@ 8 -b -o bam/$S.aligned.tmp.bam - || { echo "!!! bwa mem"; exit 1; }
  set +o pipefail
  mv bam/$S.aligned.tmp.bam bam/$S.aligned.bam
fi

if [ ! -s bam/$S.bam ]; then
  $PB bamsort --ref "$REF" --in-bam bam/$S.aligned.bam --out-bam bam/$S.qn.bam --sort-order queryname \
    > logs/$S.bamsort.log 2>&1 || { echo "!!! bamsort"; exit 1; }
  $PB markdup --ref "$REF" --in-bam bam/$S.qn.bam --out-bam bam/$S.bam \
    > logs/$S.markdup.log 2>&1 || { echo "!!! markdup"; exit 1; }
  $SAM index -@ 16 bam/$S.bam || exit 1
  rm -f bam/$S.qn.bam
fi

if [ ! -s vcf/${S}_deepvariant.vcf.gz ]; then
  $DV /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="$REF" \
    --reads=bam/$S.bam \
    --output_vcf=vcf/${S}_deepvariant.vcf.gz \
    --output_gvcf=vcf/${S}_deepvariant.g.vcf.gz \
    --num_shards=36 \
    --intermediate_results_dir=tmp/$S \
    > logs/$S.deepvariant.log 2>&1 || { echo "!!! deepvariant"; exit 1; }
  rm -rf tmp/$S
fi
echo "done: $BASE/vcf/${S}_deepvariant.vcf.gz"
