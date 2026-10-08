#!/bin/bash
# Download bulk WGS of MOLM13 (SRR17524983) from ENA and check the md5 sums.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

set -u
OUT=$RESULTS/molm13_bulk/fastq
mkdir -p "$OUT"
cd "$OUT"

B=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR175/083/SRR17524983
printf "cc0b075bf7f8db9feef804bf11972ed0  SRR17524983_1.fastq.gz\n56f592e9de2267c6251a4bb7240dc2da  SRR17524983_2.fastq.gz\n" > md5.txt
for i in 1 2; do
  wget -c -q "$B/SRR17524983_$i.fastq.gz" || { echo "!!! SRR17524983_$i"; exit 1; }
done
md5sum -c md5.txt
