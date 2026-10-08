#!/bin/bash
# Download COSMIC v104 GRCh38 files for ANNOVAR. Needs a COSMIC account:
#   export AUTH=$(echo -n 'email:password' | base64)
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

mkdir -p ${ROOT}/cosmic
cd ${ROOT}/cosmic

: "${AUTH:?Set: export AUTH=\$(echo -n 'email:password' | base64)}"
API="https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"

for path in "grch38/cosmic/v104/VCF/Cosmic_GenomeScreensMutant_Vcf_v104_GRCh38.tar" \
            "grch38/cosmic/v104/TSV/Cosmic_GenomeScreensMutant_Tsv_v104_GRCh38.tar"; do
  f=$(basename "$path")
  echo ">>> Getting URL for $f ($(date))"
  URL=$(curl -s -H "Authorization: Basic $AUTH" \
        "$API?path=$path&bucket=downloads" | grep -o 'https[^"]*')
  if [ -z "$URL" ]; then
    echo ">>> failed to get URL for $f (check credentials / auth)"; continue
  fi
  echo ">>> Downloading $f ($(date))"
  if curl -s "$URL" --output "$f"; then
    echo ">>> OK $f ($(date)) size=$(du -h "$f" | cut -f1)"
  else
    echo ">>> download failed $f"
  fi
done
echo ">>> COSMIC done ($(date))"
