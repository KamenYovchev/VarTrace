#!/bin/bash
# Download the DNA VCFs of the three cell lines from the public BioSkryb bucket.
#   export BIOSKRYB_KEY=...
#   export BIOSKRYB_SECRET=...
#   bash 01_download_dna_vcfs.sh
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail

: "${BIOSKRYB_KEY:?set BIOSKRYB_KEY}"
: "${BIOSKRYB_SECRET:?set BIOSKRYB_SECRET}"

REMOTE=":s3,provider=Wasabi,access_key_id=${BIOSKRYB_KEY},secret_access_key=${BIOSKRYB_SECRET},endpoint=${BIOSKRYB_ENDPOINT},region=${BIOSKRYB_REGION},force_path_style=true:"
BASE="$BIOSKRYB_DATA"

# MOLM13 and HG001 come called with DeepVariant and the PTA model, as used for
# the labels. HG002 comes called with Sentieon DNAscope and only lists variant
# positions, so it is called again in steps 11 and 12; these files are kept for
# the comparison.
PTA=google-deepvariant_1.8.0_ptamodel
mkdir -p "$DATA/molm13/dna_vcf/$PTA" "$DATA/hg001/dna_vcf/$PTA" "$DATA/hg002/dna_vcf"

"$RCLONE" copy "${REMOTE}${BASE}/ResolveOME_v1.0/manual/molm13/dna/vcf/$PTA/" "$DATA/molm13/dna_vcf/$PTA/" -P
"$RCLONE" copy "${REMOTE}${BASE}/ResolveOME_v2.1/manual/hg001/dna/vcf/$PTA/" "$DATA/hg001/dna_vcf/$PTA/" -P
"$RCLONE" copy "${REMOTE}${BASE}/ResolveOME_v2.0/automated/hg002/dna/vcf/sentieon-dnascope_202503.02_bioskryb129/" "$DATA/hg002/dna_vcf/" -P
