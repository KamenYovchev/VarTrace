#!/bin/bash
# Download the DNA VCFs of the three cell lines from the public BioSkryb bucket.
#   export BIOSKRYB_KEY=...
#   export BIOSKRYB_SECRET=...
#   bash 01_download_dna_vcfs.sh
set -euo pipefail

: "${BIOSKRYB_KEY:?set BIOSKRYB_KEY}"
: "${BIOSKRYB_SECRET:?set BIOSKRYB_SECRET}"

REMOTE=":s3,provider=Wasabi,access_key_id=${BIOSKRYB_KEY},secret_access_key=${BIOSKRYB_SECRET},endpoint=${BIOSKRYB_ENDPOINT},region=${BIOSKRYB_REGION},force_path_style=true:"
BASE="$BIOSKRYB_DATA"

mkdir -p dna_vcf/molm13 dna_vcf/hg002 dna_vcf/hg001

rclone copy "${REMOTE}${BASE}/ResolveOME_v1.0/manual/molm13/dna/vcf/google-deepvariant_1.8.0_ptamodel/" dna_vcf/molm13/ -P
rclone copy "${REMOTE}${BASE}/ResolveOME_v2.0/automated/hg002/dna/vcf/sentieon-dnascope_202503.02_bioskryb129/" dna_vcf/hg002/ -P
rclone copy "${REMOTE}${BASE}/ResolveOME_v2.1/manual/hg001/dna/vcf/google-deepvariant_1.8.0_ptamodel/" dna_vcf/hg001/ -P
