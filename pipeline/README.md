# Analysis pipeline

These are the scripts behind the MSc thesis "Direct Inference of Expressed SNV
Origin from Single-Cell RNA Sequencing Using Machine Learning". They go from the
raw reads to the numbers and figures in the thesis.

VarTrace, the tool, is separate. It takes processed variants and returns a
score. This folder is the processing that produces those variants and trains
the models.

## Running

    source config.sh          paths and settings, change them for another machine
    bash run_all.sh --list    show the steps without running them
    bash run_all.sh           run steps 1 to 18 in order
    bash run_all.sh 14 18     run only some of them

The steps after 18 are Python scripts that take their paths as arguments, so
they are run one at a time. The order is printed by run_all.sh --list and the
arguments are at the top of each file.

## Settings

Everything that depends on the machine is in config.sh, and every value can be
overridden from the environment:

    ROOT                 working directory, /data by default
    CONDA_BASE, ENVS     where conda and the environments are
    BIOSKRYB_ENDPOINT    address of the public BioSkryb data
    BIOSKRYB_REGION
    BIOSKRYB_DATA
    GPUS                 GPU ids for steps 11 to 13, one lane per GPU
    HG002_DIR            folder for the HG002 reprocessing and its containers
    RCLONE               rclone binary, from the rclone environment by default

Access keys are never stored in the scripts. Set them before running:

    export BIOSKRYB_KEY=...  BIOSKRYB_SECRET=...           steps 1, 2, 13, 18
    export AUTH=$(echo -n 'email:password' | base64)       step 4, COSMIC account

## Environments

Recreate them from the files in environments:

    conda env create -f environments/variant-calling.yml

The versions in these files are the ones used for the thesis.

    read-qc           fastp, FastQC, MultiQC
    variant-calling   STAR, GATK, bcftools, samtools
    annotation        ANNOVAR
    signatures        SigProfilerMatrixGenerator, SigProfilerAssignment
    ml-cpu            scikit-learn, xgboost, pandas, shap
    rclone            data download

Steps 11 to 13 and 42 run through Apptainer on a machine with a GPU, with the
Parabricks 4.4.0 and DeepVariant 1.8.0 images in HG002_DIR.

The reference is GRCh38, primary assembly, with GENCODE 49. The HG002 recall
uses GRCh38 from the GATK bundle, because the header of the BioSkryb files has
3366 contigs and only that version matches it.

## Steps

### 1. Data

    01_download_dna_vcfs.sh                    DNA VCFs of the three cell lines
    02_download_hg001_v1_cells.sh              the 10 extra HG001 v1.0 cells

The DNA BAMs (about 540 GB) are not downloaded. They are read remotely in step 18.

### 2. Reference databases

    03_download_annovar_databases.sh           refGene, avsnp150, gnomAD, dbNSFP, ClinVar
    04_download_cosmic.sh                      COSMIC v104, needs a licence

REDIportal has no automatic download and is fetched by hand. COSMIC and
REDIportal are converted to ANNOVAR tables by hand.

### 3. RNA reads

    05_trim_rna_reads.sh                       trimming with fastp, FastQC and MultiQC reports
    06_align_rna_star.sh                       STAR index if missing, two-pass alignment
    07_call_rna_variants_gatk.sh               MarkDuplicates, SplitNCigarReads, HaplotypeCaller, VariantFiltration
    08_process_hg001_v1_rna.sh                 the same three steps for the HG001 v1.0 cells

### 4. Filtering

    09_keep_biallelic_snv_main_chrom.sh        biallelic SNVs, then chr1-22, X and Y
    10_filter_rna_by_quality.sh                quality filter on the RNA variants

The quality filter keeps records with FILTER PASS, depth at least 3, at least 3
reads for the alternative allele, GQ at least 30, QD at least 2, MQ at least 50,
FS at most 30 and SOR at most 3.

### 5. HG002 DNA recall

    11_align_hg002_dna_parabricks.sh           alignment with Parabricks fq2bam
    12_call_hg002_dna_deepvariant.sh           calling with the PTA model
    13_validate_calling_reproducibility.sh     check that the chain reproduces BioSkryb's VCF

The delivered HG002 files only list positions with a variant, so most RNA
variants could not be labelled. The reproducibility check uses an HG001 cell,
DNA-OME-Manuscript-NA12878-SC1, because no aligned reads are published for HG002.

### 6. Annotation

    14_annotate_rna_annovar.sh                 RNA variants
    15_annotate_dna_annovar.sh                 DNA variants

The chr prefix is removed for ANNOVAR, its databases use 1, 2 and so on.

### 7. RNA against DNA

    16_intersect_rna_with_dna.sh / .py         left join by chr, pos, ref and alt

### 8. Labelling

    17_label_variants_by_origin.sh / .py       label by the DNA of the same cell
    18_rescue_undetermined_from_bam.sh / .py   recheck variants without a DNA record in the BAM

A variant found in the DNA (same alternative allele, non-reference genotype) is
DNA_origin. A variant not in the DNA, where the DNA covers the position with at
least five reads, is RNA_only. Reads are not counted for the alternative allele;
the rule relies on DeepVariant not having called it. Anything else is
indeterminate and left out.

Step 18 is needed because callers only write positions with a variant. It reads
the depth straight from the DNA BAMs through their index, a few kilobytes per
position. For HG002 the BAMs from step 11 in $HG002_DIR/bam are served locally
over http by the same script; set HG002_BAM_URL to use another address.

### 9. Mutational signatures

    19_make_sigprofiler_input.py               SigProfiler input per cell and per group
    20_fit_cosmic_signatures.py                SBS96 spectra and COSMIC v3.4 fitting

### 10. Training table

    21_build_feature_table.py                  the 35 features in five groups and the label
    22_add_signature_features.py               adds the 86 signature features

The cross-cell features are computed per cell line, so a variant seen in two
lines does not carry information between them.

### 11. Training and evaluation

    23_compare_algorithms_early.py             feature groups added one at a time
    24_evaluate_final_model.py                 cross-validation by variant and the final model
    25_compare_algorithms_and_lines.py         XGBoost, held-out cell line, neural net
    26_audit_leakage_and_ablation.py           leakage audit, full metrics, ablation by group
    27_train_four_shipped_models.py            the four models VarTrace ships
    28_train_partial_annotation_model.py       the fifth model, for partial annotation
    29_calibrate_partial_model.py              calibration of the fifth model

Cross-validation is split by variant, not by row. Calibration is isotonic and
nested in the same split.

### 12. Further analyses

    30_analyse_errors_and_thresholds.py        errors, operating points, rule baselines
    31_split_rna_class_by_rediportal.py        real editing inside the RNA class
    32_measure_cosmic_contribution.py          what the COSMIC feature adds
    33_measure_database_tiers.py               behaviour with different databases
    34_learning_curve_over_cells.py            would more cells help
    35_check_variants_outside_databases.py     variants outside all databases
    36_score_indeterminate_variants.py         scores for the indeterminate variants
    37_train_missing_aware_model.py            model trained with feature groups withheld
    38_match_depmap_somatic.py                 MOLM13 against DepMap

### 13. Label check

    39_download_giab_truth.sh                  GIAB v4.2.1 for HG001 and HG002
    40_check_labels_giab.py                    HG001 and HG002 labels against GIAB
    41_download_molm13_bulk.sh                 bulk WGS of MOLM13, SRR17524983
    42_process_molm13_bulk.sh                  BWA-MEM, Parabricks markdup, DeepVariant WGS model
    43_check_labels_molm13_bulk.py             MOLM13 labels against the bulk sample
    44_build_corrected_dataset.py              table with the corrected labels
    45_measure_effect_of_label_noise.py        effect of the wrong labels on ROC-AUC and calibration

Amplifying the DNA of one cell sometimes loses one copy, and an inherited
variant then looks RNA_only. These steps check the labels against DNA that was
not amplified from a single cell. The thesis analysis uses the original labels.
The VarTrace models from version 0.3.0 are steps 27 to 29 run on the table from
step 44. The models from the original labels are under the tag thesis-2026-08
of the VarTrace repository.

### 14. Figures

    46_make_thesis_figures.py                  figures of the results chapter

The learning curve figure comes from step 34. Plot labels are in Bulgarian
because the figures are in the thesis.

## Not included

Things that are placed by hand before step 5, because no script downloads them:

    RNA reads of MOLM13, HG001 v2.1 and HG002     $DATA/<line>/rna_fastq/
    HG002 DNA reads                               $HG002_DIR/fastq/
    GRCh38 primary assembly and GENCODE 49        $REFERENCE/
    GRCh38 from the GATK bundle (3366 contigs)    $HG002_DIR/ref_bioskryb/
    DeepVariant PTA model and population VCFs     $HG002_DIR/ptamodel/
    containers: Parabricks 4.4.0, DeepVariant 1.8.0, BWA, samtools

The reads are the public BioSkryb ResolveOME demo data. Step 2 downloads the
RNA reads of the 10 HG001 v1.0 cells.

The reference databases: step 3 fetches the ANNOVAR ones. COSMIC v104 needs a
licence (step 4) and has to be converted to the ANNOVAR format hg38_cosmic104
by hand, as does REDIportal (hg38_rediportal). The trained models ship with
VarTrace.
