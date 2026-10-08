#!/bin/bash
# Download the ANNOVAR hg38 databases.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

ANNOVAR=${ROOT}/annovar
HUMANDB=$ANNOVAR/humandb
mkdir -p "$HUMANDB"
cd "$ANNOVAR"

# RefGene = gene annotation; avsnp150 = dbSNP; gnomad30_genome = population freq
# (germline); dbnsfp47a = 18-predictor panel (~40GB, slow); clinvar = clinical
# significance.
for db in refGene avsnp150 gnomad30_genome dbnsfp47a clinvar_20240611; do
  echo ">>> start $db ($(date))"
  if perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar "$db" "$HUMANDB/"; then
    echo ">>> OK $db ($(date))"
  else
    echo ">>> warn failed $db ($(date))"
  fi
done
echo ">>> all DB downloads done ($(date))"
