#!/bin/bash
# Download the ANNOVAR hg38 databases.
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate annotation

ANNOVAR=${ROOT}/annovar
HUMANDB=$ANNOVAR/humandb
mkdir -p "$HUMANDB"
cd "$ANNOVAR"

# refGene and cytoBand: gene and band; avsnp150: dbSNP; gnomad30_genome and
# gnomad211_exome: population frequency; dbnsfp47a: predictors (about 40 GB);
# clinvar: clinical significance. cosmic104 and rediportal are not on the
# ANNOVAR server, see README.md.
for db in refGene cytoBand avsnp150 gnomad30_genome gnomad211_exome dbnsfp47a clinvar_20240611; do
  echo ">>> start $db ($(date))"
  if perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar "$db" "$HUMANDB/"; then
    echo ">>> OK $db ($(date))"
  else
    echo ">>> warn failed $db ($(date))"
  fi
done
echo ">>> all DB downloads done ($(date))"
