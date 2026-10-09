# Paths and environments for the pipeline.

ROOT="${ROOT:-/data}"

DATA="$ROOT/singlecell_data"      # raw reads and DNA VCFs
RESULTS="$ROOT/results"           # everything the pipeline writes
REFERENCE="$ROOT/reference"       # GRCh38 and GENCODE
ANNOVAR="$ROOT/annovar"           # ANNOVAR and its databases
COSMIC="$ROOT/cosmic"             # COSMIC v104, needs a licence

# Recreate the environments from environments/ with
#   conda env create -f environments/variant-calling.yml
CONDA_BASE="${CONDA_BASE:-/opt/miniforge3}"
ENVS="${ENVS:-$ROOT/conda/envs}"

# rclone for the bucket downloads and for serving BAMs over http
RCLONE="${RCLONE:-$ENVS/rclone/bin/rclone}"

# Public BioSkryb demo data. The keys are not here, set BIOSKRYB_KEY and
# BIOSKRYB_SECRET in the environment.
BIOSKRYB_ENDPOINT="${BIOSKRYB_ENDPOINT:-s3.us-east-1.wasabisys.com}"
BIOSKRYB_REGION="${BIOSKRYB_REGION:-us-east-1}"
BIOSKRYB_DATA="${BIOSKRYB_DATA:-bioskryb-public-data/demo-data/illumina}"

# GPU ids for the steps that need them (11, 12, 13); one lane per GPU
GPUS="${GPUS:-1 2 3 4}"

# Where the HG002 reprocessing and its containers live
HG002_DIR="${HG002_DIR:-$ROOT/hg002_reprocess}"

export ROOT DATA RESULTS REFERENCE ANNOVAR COSMIC CONDA_BASE ENVS
export BIOSKRYB_ENDPOINT BIOSKRYB_REGION BIOSKRYB_DATA HG002_DIR GPUS RCLONE
