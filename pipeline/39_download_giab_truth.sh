#!/bin/bash
# Download the GIAB v4.2.1 benchmark for HG001 and HG002 (GRCh38, chr1-22).
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u
OUT=$RESULTS/giab_label_noise/ref
mkdir -p "$OUT"
cd "$OUT"

B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
for f in \
  NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
  NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
  AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
  AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
do
  [ -s "$(basename "$f")" ] || curl -s -f -O "$B/$f" || { echo "!!! $f"; exit 1; }
done
ls -la
