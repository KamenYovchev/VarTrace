#!/bin/bash
# Download the 10 extra HG001 cells (ResolveOME v1.0).
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u

# Bucket keys come from the environment.
: "${BIOSKRYB_KEY:?set BIOSKRYB_KEY}"
: "${BIOSKRYB_SECRET:?set BIOSKRYB_SECRET}"
REMOTE=":s3,provider=Wasabi,access_key_id=${BIOSKRYB_KEY},secret_access_key=${BIOSKRYB_SECRET},endpoint=${BIOSKRYB_ENDPOINT},region=${BIOSKRYB_REGION},force_path_style=true:"
SRC="$BIOSKRYB_DATA/ResolveOME_v1.0/manual/hg001"
DEST=$DATA/hg001_v1

mkdir -p "$DEST/dna_vcf" "$DEST/rna_fastq"

echo ">>> DNA VCFs (DeepVariant 1.8.0 ptamodel) ..."
"$RCLONE" copy "${REMOTE}${SRC}/dna/vcf/google-deepvariant_1.8.0_ptamodel/" "$DEST/dna_vcf/" -P

echo ">>> RNA fastq ..."
"$RCLONE" copy "${REMOTE}${SRC}/rna/fastq/" "$DEST/rna_fastq/" -P

echo ""
echo "check"
echo "DNA VCF  : $(ls "$DEST"/dna_vcf/*.vcf.gz 2>/dev/null | wc -l) / 10"
echo "RNA fastq: $(ls "$DEST"/rna_fastq/*.fastq.gz 2>/dev/null | wc -l) / 20  (10 cells, paired)"
du -sh "$DEST"
echo ""
echo "note: DNA BAMs stay in the bucket at"
echo "  ${SRC}/dna/bam/parabricks_4.4.0/"
echo "and are streamed on demand (rclone serve http + samtools), never downloaded."
