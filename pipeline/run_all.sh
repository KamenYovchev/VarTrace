#!/bin/bash
# Driver: runs the pipeline in order.
#   for e in environments/*.yml; do conda env create -f "$e"; done
#   export BIOSKRYB_KEY=...  BIOSKRYB_SECRET=...
#   bash run_all.sh            run everything
#   bash run_all.sh 14 20      run steps 14 to 20 only
#   bash run_all.sh --list     only show the steps
set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
source config.sh

FROM="${1:-1}"
TO="${2:-99}"
LIST=0
if [ "${1:-}" = "--list" ]; then LIST=1; FROM=1; TO=99; fi

GPU_STEPS="11 12 13"

run_step () {
  local n="$1" script="$2" what="$3"
  local num=$((10#$n))
  [ "$num" -lt "$((10#$FROM))" ] && return 0
  [ "$num" -gt "$((10#$TO))" ] && return 0
  if [ "$LIST" = "1" ]; then
    printf "%2s  %-42s %s\n" "$n" "$script" "$what"
    return 0
  fi
  if [[ " $GPU_STEPS " == *" $num "* ]] && [ "${WITH_GPU:-0}" != "1" ]; then
    echo ">>> $n skipped (needs a GPU, run with WITH_GPU=1)"
    return 0
  fi
  echo ">>> $n  $what  ($(date '+%H:%M:%S'))"
  bash "$script" || { echo "!!! $n failed: $script"; exit 1; }
}

run_step 01 01_download_dna_vcfs.sh                "download the DNA VCFs"
run_step 02 02_download_hg001_v1_cells.sh          "download the ten HG001 v1.0 cells"
run_step 03 03_download_annovar_databases.sh       "download the ANNOVAR databases"
run_step 04 04_download_cosmic.sh                  "download COSMIC v104"
run_step 05 05_trim_rna_reads.sh                   "trim RNA reads with fastp"
run_step 06 06_align_rna_star.sh                   "two-pass STAR alignment"
run_step 07 07_call_rna_variants_gatk.sh           "call variants with GATK"
run_step 08 08_process_hg001_v1_rna.sh             "same three steps for the HG001 v1.0 cells"
run_step 09 09_keep_biallelic_snv_main_chrom.sh    "biallelic SNVs, main chromosomes"
run_step 10 10_filter_rna_by_quality.sh            "quality filter on RNA variants"
run_step 11 11_align_hg002_dna_parabricks.sh       "align HG002 DNA"
run_step 12 12_call_hg002_dna_deepvariant.sh       "call HG002 DNA"
run_step 13 13_validate_calling_reproducibility.sh "reproducibility check"
run_step 14 14_annotate_rna_annovar.sh             "annotate RNA"
run_step 15 15_annotate_dna_annovar.sh             "annotate DNA"
run_step 16 16_intersect_rna_with_dna.sh           "intersect RNA with DNA"
run_step 17 17_label_variants_by_origin.sh         "label by origin"
run_step 18 18_rescue_undetermined_from_bam.sh     "rescue indeterminate by coverage"

# The Python steps take their paths as arguments.
if [ "$LIST" = "1" ] || { [ "$((10#$FROM))" -le 19 ] && [ "$((10#$TO))" -ge 19 ]; }; then
  cat <<'TXT'

The steps after labelling are run one at a time, each with its own arguments.
Environment ml-cpu, except 19 and 20 which need signatures.

 19  19_make_sigprofiler_input.py            input for SigProfiler
 20  20_fit_cosmic_signatures.py             SBS96 spectra and COSMIC fitting
 21  21_build_feature_table.py               table with the 35 features and the label
 22  22_add_signature_features.py            adds the 86 signature features
 23  23_compare_algorithms_early.py          early algorithm comparison
 24  24_evaluate_final_model.py              cross-validation, figures, model
 25  25_compare_algorithms_and_lines.py      XGBoost, held-out line, neural net
 26  26_audit_leakage_and_ablation.py        leakage audit and ablation by group
 27  27_train_four_shipped_models.py         the four models for the tool
 28  28_train_partial_annotation_model.py    the fifth model
 29  29_calibrate_partial_model.py           calibration of the fifth model
 30  30_analyse_errors_and_thresholds.py     errors, operating points, rule baselines
 31  31_split_rna_class_by_rediportal.py     real editing inside the RNA class
 32  32_measure_cosmic_contribution.py       what the COSMIC feature adds
 33  33_measure_database_tiers.py            behaviour with different databases
 34  34_learning_curve_over_cells.py         learning curve over cells
 35  35_check_variants_outside_databases.py  variants outside all databases
 36  36_score_indeterminate_variants.py      scores for the indeterminate variants
 37  37_train_missing_aware_model.py         model trained with missing groups
 38  38_match_depmap_somatic.py              comparison with DepMap

Label check against DNA not amplified from a single cell:

 39  39_download_giab_truth.sh               GIAB v4.2.1 for HG001 and HG002
 40  40_check_labels_giab.py                 HG001 and HG002 labels against GIAB
 41  41_download_molm13_bulk.sh              MOLM13 bulk genome, SRR17524983
 42  42_process_molm13_bulk.sh               alignment and calling, needs a GPU
 43  43_check_labels_molm13_bulk.py          MOLM13 labels against the bulk sample
 44  44_build_corrected_dataset.py           table with corrected labels
 45  45_measure_effect_of_label_noise.py     effect of the noise on the numbers
 46  46_make_thesis_figures.py               figures of the results chapter

The VarTrace models from version 0.3.0 on are steps 27, 28 and 29 run on the table from step 44.
Arguments are at the top of each file.
TXT
fi

[ "$LIST" = "1" ] || echo ">>> done"
