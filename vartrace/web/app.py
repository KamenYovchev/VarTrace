"""The web interface. Started with: vartrace serve

Uploaded files are written to a temporary directory, read, and deleted as soon
as the scores are computed. Nothing is kept on the server.
"""

import os
import tempfile
import zipfile
from io import BytesIO

import numpy as np
import pandas as pd
import streamlit as st

from vartrace import __version__, cli, features, io, model as models, report
from vartrace.web import charts

ABSENT, PRESENT = charts.ABSENT, charts.PRESENT


@st.cache_resource(show_spinner=False)
def loaded_models():
    """Kept across runs, since the five of them are seven megabytes."""
    return models.load_all()


def which_model(variants, use_ensembl):
    """What will be used, and what Ensembl would change about that.

    Worth saying before anything is scored, because the choice follows from
    the input rather than from anything the person asked for, and because the
    difference between the strongest and the weakest of the five is large.
    """
    available = features.available(variants)
    chosen, why = models.choose(loaded_models(), available)
    now = {"name": chosen.name, "features": len(chosen.features), "why": why}

    later = None
    if use_ensembl and not (available & features.DB_GROUPS):
        after, why_after = models.choose(loaded_models(),
                                         available | {"dbsnp", "gnomad"})
        later = {"name": after.name, "features": len(after.features),
                 "why": why_after}
    return now, later


@st.cache_data(show_spinner=False)
def inspect(files, use_ensembl):
    """Read the files and check them, without asking Ensembl or scoring.

    Kept apart from the scoring so the person sees what arrived and agrees to
    it before anything slow or irreversible happens, and so a refusal costs
    nothing.
    """
    tmp = tempfile.mkdtemp(prefix="vartrace_")
    paths = []
    try:
        for name, data in files:
            path = os.path.join(tmp, os.path.basename(name))
            with open(path, "wb") as fh:
                fh.write(data)
            paths.append(path)

        try:
            variants, n_cells, skipped = io.read_all(paths)
        except io.UnknownFormat as problem:
            return {"refused": str(problem)}
        if not variants:
            return None

        features.compute(variants, n_cells)
        now, later = which_model(variants, use_ensembl)
        return {
            "problems": cli.check(paths, variants),
            "model_now": now,
            "model_later": later,
            "facts": {
                "files": len(paths),
                "kinds": sorted({"ANNOVAR table" if io.looks_annotated(p)
                                 else "VCF" for p in paths}),
                "cells": n_cells,
                "variants": len(variants),
                "skipped": skipped,
                "databases": sorted(set().union(*(v.databases
                                                  for v in variants))),
                "builds": sorted({b for b in (io.genome_build(p)
                                              for p in paths) if b}),
                "blind": features.without_allele_counts(variants),
                "plain": features.without_read_level(variants),
            },
        }
    finally:
        for path in paths:
            os.remove(path)
        os.rmdir(tmp)


@st.cache_data(show_spinner=False)
def compute(files, use_ensembl, anyway):
    """Files are (name, bytes). Cached so moving the slider does not rerun it."""
    tmp = tempfile.mkdtemp(prefix="vartrace_")
    paths = []
    try:
        for name, data in files:
            path = os.path.join(tmp, os.path.basename(name))
            with open(path, "wb") as fh:
                fh.write(data)
            paths.append(path)

        try:
            variants, n_cells, skipped = io.read_all(paths)
        except io.UnknownFormat as problem:
            return {"refused": str(problem)}
        if not variants:
            return None
        messages = []
        try:
            model, why, scores, annotation = cli.score(
                variants, n_cells, models.MODELS_DIR, use_ensembl,
                log=messages.append, paths=paths, anyway=anyway)
        except cli.Refused as refusal:
            return {"problems": refusal.problems}

        table = pd.DataFrame([{
            "cell": v.cell,
            "chrom": v.chrom,
            "pos": int(v.pos),
            "ref": v.ref,
            "alt": v.alt,
            "substitution": f"{v.ref}>{v.alt}",
            "vaf": v.features["vaf"],
            "dp": v.features["dp"],
            "score": float(p),
        } for v, p in zip(variants, scores)])
        return {
            "table": table, "n_cells": n_cells, "skipped": skipped, "why": why,
            "messages": messages, "model": model.file,
            "provenance": report.provenance(model, why, "see below", annotation),
            "annotation": annotation,
        }
    finally:
        for path in paths:
            os.remove(path)
        os.rmdir(tmp)


PRETTY = {"dbsnp": "dbSNP", "gnomad": "gnomAD", "redi": "REDIportal",
          "cosmic": "COSMIC"}

# one entry per thing the input has to satisfy, each opened from the page
# rather than read as a block of prose nobody finishes
REQUIREMENTS = [
    {
        "name": "GRCh38",
        "checked": True,
        "status": "Checked. VarTrace refuses without it.",
        "body": [
            "Your variants must sit on the GRCh38 assembly. VarTrace reads "
            "which one it is from the contig lengths in the VCF header: "
            "chromosome 1 is 248,956,422 bases in GRCh38 and 249,250,621 in "
            "GRCh37.",
            "An ANNOVAR table carries no header, so there the assembly cannot "
            "be told and GRCh38 is assumed.",
            "On GRCh37 every position means a different place in the genome, "
            "and both the models and the Ensembl lookup would be answering "
            "about somewhere else. Lift the variants over with CrossMap or "
            "the Ensembl Assembly Converter first.",
        ],
    },
    {
        "name": "single-cell RNA",
        "checked": False,
        "status": "Not checked. This one is down to you.",
        "body": [
            "The variants must come from RNA sequencing of single cells. Not "
            "bulk RNA, and not DNA.",
            "Nothing in a VCF records which it was, so VarTrace cannot tell "
            "and will score bulk data without complaint. Four of the "
            "thirty-five features describe how a variant behaves across "
            "cells, and on bulk data they mean something else entirely.",
        ],
    },
    {
        "name": "GATK-called",
        "checked": True,
        "status": "Checked. VarTrace refuses when they are missing.",
        "body": [
            "The INFO field must carry FS, SOR, MQ, MQRankSum, "
            "ReadPosRankSum, BaseQRankSum and QD. Seven of the thirty-five "
            "features are read straight from there.",
            "GATK writes them. DeepVariant, bcftools and Clair3 do not. "
            "DRAGEN has a GATK-compatible mode; whether yours wrote them is "
            "answered by looking, which is what VarTrace does rather than "
            "trusting the name of the caller.",
            "Unlike the database features, there is no model trained without "
            "these, which is why their absence stops the run instead of "
            "choosing a different model.",
        ],
    },
    {
        "name": "AD or VAF",
        "checked": True,
        "status": "Checked. VarTrace refuses without it.",
        "body": [
            "The sample column must carry allele depths (AD) or an allele "
            "fraction (VAF). One of them gives the variant allele fraction, "
            "which is the strongest single feature after the databases.",
            "Some filtering and conversion steps quietly drop AD. If a file "
            "that once had it arrives without it, that is usually where it "
            "went.",
        ],
    },
    {
        "name": "DP >= 5",
        "checked": True,
        "status": "Checked. VarTrace refuses when too many are shallower.",
        "body": [
            "The training data has no variant under 4 reads and few at 4 "
            "(651 of 12,916 rows), because the quality filter asks for depth "
            "3 and at least 3 reads with the alternative allele. Below 5 the "
            "models have seen little, so filter yours like this:",
        ],
        "command": "bcftools view -e 'FMT/DP<5' in.vcf -o out.vcf",
        "after": "A few shallow variants are tolerated. The run stops when "
                 "more than a fifth of the input is below 5.",
    },
]


def requirement(need, with_name=False):
    """One requirement, either opened from a chip or listed in full."""
    status = f"{need['name']}. {need['status']}" if with_name else need["status"]
    (st.success if need["checked"] else st.warning)(status)
    for paragraph in need["body"]:
        st.write(paragraph)
    if need.get("command"):
        st.code(need["command"], language="bash")
    if need.get("after"):
        st.write(need["after"])

# the numbers behind the tool, each opening on what it actually covers
FEATURE_GROUPS = pd.DataFrame([
    ("allele counts", 6, "reference and alternative depths, the allele "
                         "fraction, total depth, variant quality, genotype "
                         "quality"),
    ("read-level", 7, "FS, SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum, "
                      "QD, all from INFO"),
    ("substitution", 13, "one for each of the twelve base changes, plus "
                         "whether it is a transition"),
    ("cross-cell", 4, "how many cells carry it, what share that is, and the "
                      "mean and spread of the allele fraction across them"),
    ("database", 5, "dbSNP, gnomAD, REDIportal and COSMIC"),
], columns=["group", "features", "what they are"])

# Written as lines rather than as a table, because the sidebar is too narrow
# for four columns and a table that has to be scrolled sideways is not read.
ACCEPTED = """
1. A VCF, from any caller and any annotator, or none at all. One sample, or
   many with a column for each cell barcode.
2. A table from ANNOVAR, the output of table_annovar.pl with its Otherinfo
   columns. Always a single cell, so one file for each.
"""

PREPARING = """\
# the same steps as the training data, see pipeline/ in the repository

# 1. trim and align the single-cell RNA reads to GRCh38, two-pass
fastp -i cell_01_R1.fq.gz -I cell_01_R2.fq.gz -o r1.fq.gz -O r2.fq.gz \\
  --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 25
STAR --genomeDir GRCh38 --readFilesIn r1.fq.gz r2.fq.gz --readFilesCommand zcat \\
  --twopassMode Basic --outSAMtype BAM SortedByCoordinate \\
  --outSAMattrRGline ID:cell_01 SM:cell_01 PL:ILLUMINA

# 2. call with GATK, which writes the seven read-level fields
gatk MarkDuplicates -I Aligned.sortedByCoord.out.bam -O md.bam -M md.txt
gatk SplitNCigarReads -R GRCh38.fa -I md.bam -O split.bam
gatk HaplotypeCaller -R GRCh38.fa -I split.bam -O raw.vcf.gz \\
  --dont-use-soft-clipped-bases true \\
  --standard-min-confidence-threshold-for-calling 20
gatk VariantFiltration -R GRCh38.fa -V raw.vcf.gz -O filt.vcf.gz \\
  -window 35 -cluster 3 --filter-name FS -filter "FS > 30.0" \\
  --filter-name QD -filter "QD < 2.0"

# 3. biallelic substitutions on chr1-22, X, Y, then the quality filter
bcftools view -v snps -m2 -M2 filt.vcf.gz \\
  | bcftools view -f PASS -i 'FMT/DP>=3 && FMT/AD[0:1]>=3 && FMT/GQ>=30 &&
      INFO/QD>=2 && INFO/MQ>=50 && INFO/FS<=30 && INFO/SOR<=3' \\
  -o cell_01.ready.vcf

# 4. optional: annotate with ANNOVAR, through convert2annovar so the
#    Otherinfo columns hold the VCF fields
convert2annovar.pl -format vcf4 -includeinfo cell_01.ready.vcf \\
  | sed 's/^chr//' > cell_01.avinput
table_annovar.pl cell_01.avinput humandb/ -buildver hg38 -out cell_01 \\
  -protocol avsnp150,gnomad30_genome,gnomad211_exome,cosmic104,rediportal \\
  -operation f,f,f,f,f -nastring . -otherinfo

# 5. upload cell_01.ready.vcf, or cell_01.hg38_multianno.txt
"""

ROUTING = """
ANNOVAR table, all four databases
- several cells: database-informed, with cross-cell, 35 features
- one cell: database-informed, no cross-cell, 31 features

ANNOVAR table, only some of them
- either way: partial annotation, 35 features

Any VCF, Ensembl switched on
- either way: partial annotation, 35 features

Any VCF, Ensembl switched off
- several cells: data-only, with cross-cell, 30 features
- one cell: data-only, no cross-cell, 26 features
"""

MODEL_TABLE = pd.DataFrame([
    ("database-informed, with cross-cell", 35, "0.9874"),
    ("database-informed, no cross-cell", 31, "0.9834"),
    ("partial annotation", 35, "0.9676 to 0.9864"),
    ("data-only, with cross-cell", 30, "0.8996"),
    ("data-only, no cross-cell", 26, "0.7797"),
], columns=["model", "features", "PR-AUC"])

STATS = [
    {
        "number": "35", "label": "features",
        "lead": "Thirty-five numbers are computed for every variant, in five "
                "groups.",
        "table": FEATURE_GROUPS,
        "body": [
            "A feature that cannot be computed stays empty rather than "
            "becoming zero. A variant with no gnomAD column is not a variant "
            "absent from gnomAD, and which model is chosen turns on that "
            "difference.",
        ],
    },
    {
        "number": "5", "label": "models",
        "lead": "One is chosen automatically from what your input supports, "
                "and the choice is stated in words with every run.",
        "table": MODEL_TABLE,
        "body": [
            "The first four are gradient boosting. The fifth is XGBoost "
            "trained with feature groups withheld at random, which is why it "
            "alone tolerates a partial input.",
            "All five are calibrated with isotonic regression, so a score of "
            "0.87 means 0.87 and not merely more than 0.80.",
            "Which one your own files will land on is in the sidebar, under "
            "Which model runs on what.",
        ],
    },
    {
        "number": "45", "label": "cells trained on",
        "lead": "45 single cells from three cell lines: HG001 and HG002, two "
                "healthy reference genomes, and MOLM13, a line from acute "
                "myeloid leukaemia.",
        "body": [
            "12,916 labelled rows over 6,317 distinct variants, of which 38.2 "
            "percent arose at the RNA level.",
            "The truth comes from sequencing the DNA of the same cell. A "
            "variant found in the DNA too is DNA_origin. One absent from the "
            "DNA, where the DNA was covered at least five times and no "
            "alternative allele was called there, is RNA_only. Where the coverage did not reach, the "
            "variant was left out rather than guessed at.",
            "Amplifying the DNA of one cell sometimes loses one of the two "
            "copies, and an inherited variant on the lost copy then looks "
            "RNA_only. Such labels were found by checking against DNA that "
            "was not amplified from a single cell: the Genome in a Bottle "
            "truth sets for HG001 and HG002 and bulk genome sequencing of "
            "MOLM13. 461 of them were corrected to DNA_origin before "
            "training.",
        ],
    },
    {
        "number": "0.99", "label": "PR-AUC",
        "lead": "0.9874, cross-validated with the folds grouped by variant, "
                "so no variant appears in both training and testing. Measured "
                "with all four databases present and more than one cell, "
                "which is the best case.",
        "body": [
            "PR-AUC rather than ROC-AUC because the classes are uneven. 38.2 "
            "percent of the data arose at the RNA level, and that prevalence "
            "is what a model knowing nothing would score.",
            "The matching precision is 0.9711 at a threshold of 0.5 and "
            "0.9879 at 0.867, the high confidence setting. That threshold "
            "was chosen on the original labels, where it gave 0.99.",
            "These figures are from the three cell lines above. What the "
            "models would score on primary tissue or another technology is "
            "not known.",
        ],
    },
]


def statistic(stat):
    """One of the headline numbers, opened from the page."""
    st.write(stat["lead"])
    if stat.get("table") is not None:
        st.dataframe(stat["table"], hide_index=True)
    for paragraph in stat["body"]:
        st.write(paragraph)


def which_model_and_why(now, later):
    """Name the model that will run, and say what decided it."""
    with st.container(border=True):
        if later:
            st.write(f"With nothing but the files, the model would be "
                     f"{now['name']} ({now['features']} features), because "
                     f"{now['why'].split(', so ')[0]}.")
            st.write(f"After Ensembl answers, it becomes {later['name']} "
                     f"({later['features']} features), because the input will "
                     f"then have dbSNP and gnomAD.")
        else:
            st.write(f"The model will be {now['name']} ({now['features']} "
                     f"features), because {now['why'].split(', so ')[0]}.")
        st.caption("The five models differ in which features they were "
                   "trained to expect. The one that matches your input is "
                   "chosen for you, and the reason is kept with the results.")


def what_was_read(f):
    """Say what the input turned out to be, before any number is shown.

    A file read as one cell when it holds three hundred, or scores computed
    without allele fractions, both give results that look perfectly ordinary.
    Stating it here is what makes such a mistake visible.
    """
    files = f"{f['files']} file" + ("s" if f["files"] != 1 else "")
    cells = f"{f['cells']} cell" + ("s" if f["cells"] != 1 else "")
    lines = [f"Read {files} ({', '.join(f['kinds'])}) as {cells}, "
             f"{f['variants']} variants."]

    if "GRCh37" in f["builds"]:
        lines.append("GRCh37, from the contig lengths in the header.")
    elif f["builds"] == ["GRCh38"]:
        lines.append("GRCh38, from the contig lengths in the header.")
    else:
        lines.append("The assembly could not be told from the files, so "
                     "GRCh38 is assumed.")

    if f["databases"]:
        lines.append("Database annotation in the files: "
                     + ", ".join(PRETTY[d] for d in f["databases"]) + ".")
    else:
        lines.append("No database annotation in the files. dbSNP and gnomAD "
                     "will be looked up at Ensembl.")

    lines.append("Allele fractions present for every variant." if not f["blind"]
                 else f"No allele fraction for {f['blind']} of {f['variants']}.")
    lines.append("Read-level annotations present for every variant."
                 if not f["plain"] else
                 f"No read-level annotations for {f['plain']} of {f['variants']}.")

    if f["skipped"]:
        lines.append(f"{f['skipped']} rows were not biallelic single base "
                     "substitutions and were skipped.")

    with st.container(border=True):
        for line in lines:
            st.write(line)


LOGO = os.path.join(os.path.dirname(__file__), "logo.svg")

st.set_page_config(page_title="VarTrace", page_icon=LOGO, layout="wide")

# Flat and plain on purpose. Tools in this field are read all day next to a
# terminal and a genome browser, so the page stays out of the way: a solid bar
# at the top, square corners, one accent colour, and no decoration anywhere.
st.markdown("""
<style>
.block-container { padding-top: 3rem; }
/* Streamlit's captions are pale enough to read as disabled text; these are
   real instructions and have to be legible */
[data-testid="stCaptionContainer"],
[data-testid="stCaptionContainer"] * { color: #3E4A54 !important; }
[data-testid="stFileUploaderDropzone"] {
  min-height: 6.5rem;
  border: 1px dashed #B6C0BA;
  border-radius: 8px;
  background: #FAFBFA;
}
/* a shadow this faint is not decoration; it is what separates a panel from
   the page without drawing another line */
[data-testid="stVerticalBlockBorderWrapper"] {
  border-radius: 10px;
  box-shadow: 0 1px 2px rgba(16, 24, 40, 0.04), 0 1px 3px rgba(16, 24, 40, 0.05);
}
[data-testid="stExpander"] details { border-radius: 8px; }
/* the one expander that belongs with the requirement buttons is coloured like
   them, so it reads as part of the same row rather than as a stray panel */
.st-key-needs summary {
  background: #24443A !important; border-color: #24443A !important;
  border-radius: 8px !important;
}
.st-key-needs summary,
.st-key-needs summary * { color: #FFFFFF !important; fill: #FFFFFF !important; }
.st-key-needs summary:hover {
  background: #2F6B45 !important; border-color: #2F6B45 !important;
}
/* the questions in the sidebar are headings you open, so they carry the same
   colour as everything else on the page that does that */
[data-testid="stSidebar"] summary {
  background: #24443A !important; border-color: #24443A !important;
  border-radius: 8px !important;
}
[data-testid="stSidebar"] summary,
[data-testid="stSidebar"] summary * {
  color: #FFFFFF !important; fill: #FFFFFF !important; font-weight: 600;
}
[data-testid="stSidebar"] summary:hover {
  background: #2F6B45 !important; border-color: #2F6B45 !important;
}
[data-testid="stSidebar"] [data-testid="stCaptionContainer"] * {
  color: #3E4A54 !important;
}
.vt-bar {
  display: flex; align-items: center; gap: 0.85rem;
  margin-bottom: 1.4rem; padding: 0.9rem 1.3rem;
  border-radius: 10px;
  background: #24443A; color: #FFFFFF;
  box-shadow: 0 1px 3px rgba(16, 24, 40, 0.10);
}
.vt-bar svg { width: 30px; height: 30px; flex: none; }
/* the mark itself is the same dark green as the step numbers, which is right
   for the browser tab but invisible on the bar, so there it is lightened */
.vt-bar svg path { stroke: #DCEAE1 !important; }
.vt-bar-title {
  font-size: 1.12rem; font-weight: 700; letter-spacing: 0.03em;
}
.vt-bar-sub { font-size: 0.83rem; color: #D6E3DB; }
.vt-bar-right { margin-left: auto; font-size: 0.76rem; color: #A9C2B4; }
.vt-need {
  display: inline-block; padding: 0.18rem 0.5rem; margin: 0 0.35rem 0.3rem 0;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 0.76rem; color: #34423A; border-radius: 5px;
  background: #F1F4F2; border: 1px solid #D4DAD6;
}
/* the requirements and the headline numbers are popovers, so that each opens
   on what it means; these two rules give the buttons back the plain look they
   had as text, since a row of grey buttons says nothing about what it is */
[class*="st-key-need-"] button, [class*="st-key-stat-"] button {
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace !important;
  color: #FFFFFF !important; background: #24443A !important;
  border: 1px solid #24443A !important; border-radius: 6px !important;
  min-height: 0 !important; font-variant-numeric: tabular-nums;
}
[class*="st-key-need-"] button {
  font-size: 0.76rem !important; padding: 0.24rem 0.62rem !important;
}
[class*="st-key-stat-"] button {
  font-size: 1.25rem !important; font-weight: 700 !important;
  padding: 0.25rem 0.8rem !important;
}
[class*="st-key-need-"] button:hover, [class*="st-key-stat-"] button:hover {
  background: #2F6B45 !important; border-color: #2F6B45 !important;
}
/* a small i, so it is plain that these are not labels but open on something */
[class*="st-key-need-"] button::after,
[class*="st-key-stat-"] button::after {
  content: "i"; display: inline-block; margin-left: 0.5rem;
  width: 1em; height: 1em; line-height: 0.95em; text-align: center;
  font-size: 0.72em; font-style: italic; font-family: Georgia, serif;
  border: 1px solid rgba(255, 255, 255, 0.6); border-radius: 50%;
  vertical-align: middle; opacity: 0.85;
}
/* the chevron is a Material icon rather than an svg, and it turns a label into
   something that looks like a half-open dropdown */
[class*="st-key-need-"] button [data-testid="stIconMaterial"],
[class*="st-key-stat-"] button [data-testid="stIconMaterial"],
[class*="st-key-need-"] button svg,
[class*="st-key-stat-"] button svg { display: none !important; }
.vt-step { display: flex; gap: 0.75rem; margin-bottom: 1rem; }
.vt-step-n {
  flex: none; width: 1.6rem; height: 1.6rem; background: #24443A;
  color: #FFFFFF; font-size: 0.78rem; font-weight: 700; border-radius: 7px;
  display: flex; align-items: center; justify-content: center;
}
.vt-step-t { font-weight: 600; font-size: 0.9rem; line-height: 1.5rem; }
.vt-step-d { font-size: 0.83rem; color: #4A555E; margin-top: 0.12rem; }
.vt-stats { display: flex; gap: 2.2rem; margin-top: 0.4rem; }
.vt-stat-n {
  font-size: 1.45rem; font-weight: 700; color: #24443A; line-height: 1.2;
  font-variant-numeric: tabular-nums;
}
.vt-stat-l {
  font-size: 0.72rem; color: #4A555E; text-transform: uppercase;
  letter-spacing: 0.06em;
}
</style>
""", unsafe_allow_html=True)

with st.sidebar:
    # one question per heading, all shut to begin with, because a wall of
    # explanation beside the thing you came to do is read by nobody
    with st.expander("What question does this answer?"):
        st.write(
            "A variant seen in RNA may be there in the genome too, or it may "
            "have arisen only at the RNA level, through ADAR editing, an error "
            "of reverse transcription, or sequencing itself. Telling the two "
            "apart normally needs the DNA of the same cell sequenced as well. "
            "VarTrace answers from the RNA alone.")

    with st.expander("What does a score mean?"):
        st.write(
            "The probability that the variant is absent from the genome of "
            "that cell. It is calibrated, not a ranking: of a hundred variants "
            "scoring 0.87, about 87 really are absent. The threshold beside "
            "the results only turns scores into calls and never changes them.")

    with st.expander("How does it decide?"):
        st.write(
            "Thirty-five features per variant, in five groups: the allele "
            "counts, the read-level annotations, the substitution type, how "
            "the variant behaves across cells, and what four public databases "
            "say about the position.")
        st.write(
            "Five models come with the package and one is chosen from what "
            "your input can support. The choice is stated in words with every "
            "run, and recorded in the run's provenance.")

    with st.expander("What files can I upload?"):
        st.write("Two kinds, and nothing else.")
        st.markdown(ACCEPTED)
        st.write(
            "Upload as many as you have. One file for each cell, or a single "
            "VCF with a column per barcode, or one file if that is all there "
            "is. VarTrace counts the cells from the files themselves.")
        st.write(
            "A table from VEP, SnpEff or VarSeq is refused, because their "
            "columns are not the ones VarTrace reads and such tables usually "
            "drop the fields it needs. Export a VCF from those instead, which "
            "works whoever made it.")

    with st.expander("Which model runs on what?"):
        st.write(
            "You choose nothing. The model follows from what your files "
            "carry, and the reason is written into the results.")
        st.markdown(ROUTING)
        st.write(
            "All four databases means an ANNOVAR table carrying avsnp150, AF, "
            "cosmic104 and rediportal. Any VCF counts as carrying none, "
            "whoever annotated it, because only ANNOVAR's columns are read.")
        st.write(
            "Several cells means more than one in total, whether from many "
            "files or from one VCF with barcode columns. With a single cell "
            "the four cross-cell features cannot be computed at all.")
        st.caption("What each model is worth is under the number 5 on the "
                   "start page.")

    with st.expander("What was it trained on?"):
        st.write(
            "45 single cells from three cell lines, HG001, HG002 and MOLM13, "
            "with the DNA of the same cells sequenced as the truth. 12,916 "
            "labelled rows over 6,317 variants, of which 38.2 percent arose "
            "at the RNA level. Cross-validated with the folds grouped by "
            "variant, so no variant appears in both training and testing.")

    with st.expander("Where is it untested?"):
        st.write(
            "Primary tissue, solid tumours, other sequencing technologies and "
            "other variant callers. Two of the three lines are healthy "
            "reference genomes and the third is a leukaemia line. If your data "
            "fall outside that, the scores may be off in ways nothing here "
            "can show.")

    with st.expander("Where do the database features come from?"):
        st.write(
            "From ANNOVAR's columns, or from Ensembl. Nothing else is read. A "
            "VCF annotated with VEP, SnpEff or VarSeq is treated as the plain "
            "VCF it still is, and dbSNP and gnomAD are fetched from Ensembl "
            "instead.")
        st.write(
            "The two sources were compared on the data behind the models: "
            "using Ensembl rather than ANNOVAR cost 0.0009 PR-AUC and changed "
            "the call for 0.87 percent of variants.")

    st.markdown("[Project source](https://codeberg.org/KamenYovchev/vartrace)")
    st.caption("Code under PolyForm Noncommercial 1.0.0, models under "
               "CC BY-NC 4.0. Uploaded files are deleted as soon as the "
               "scores are computed.")

# The input sits in the middle of the page rather than in the sidebar, because
# it is the one thing anyone opening this has come to do. Everything that
# explains it stays below the fold, short, and out of the way once files are in.
with open(LOGO) as fh:
    mark = fh.read().split("?>")[-1].strip()

st.markdown(f"""
<div class="vt-bar">
  {mark}
  <div>
    <div class="vt-bar-title">VarTrace</div>
    <div class="vt-bar-sub">Variants seen in single-cell RNA that are absent
    from the genome of the cell</div>
  </div>
  <div class="vt-bar-right">v{__version__}</div>
</div>
""", unsafe_allow_html=True)

# Choosing the input and reading the results are two different jobs, so they
# are two screens. Once a run is agreed to, the files are kept in the session
# and the upload card is gone until the run is ended.
started = st.session_state.get("agreed", False)

if started:
    files = st.session_state["files"]
    use_ensembl = st.session_state["use_ensembl"]
    st.subheader("Results")

middle, aside = st.columns([3, 2], gap="large")
if not started:
    with middle:
        with st.container(border=True):
            # one field, no choice to make: the number of cells is counted
            # from the files themselves, so asking about it beforehand only
            # invited a wrong answer
            uploaded = st.file_uploader("Your variants",
                                        type=["vcf", "tsv", "txt"],
                                        accept_multiple_files=True,
                                        key="many")
            st.caption("One file for each cell, selected together, or a "
                       "single VCF with a column per cell barcode. Plain VCF "
                       "or a table from ANNOVAR. VarTrace counts the cells "
                       "itself.")

            # asked here rather than beside the results, because it decides
            # what the run does and cannot be changed after the fact
            use_ensembl = st.toggle(
                "Fetch dbSNP and gnomAD from Ensembl", value=True,
                help="Positions and alleles are the only thing that leaves "
                     "this machine. Turning it off works offline, at the cost "
                     "of a weaker model.")
            # worth saying here rather than in the documentation, because an
            # annotated VCF looks like a file that would not need this, and it
            # is the one case people get wrong
            st.caption("Only an ANNOVAR table counts as already annotated. A "
                       "VCF is sent to Ensembl even if VEP, SnpEff or VarSeq "
                       "annotated it, because their fields are not read.")
            if not use_ensembl:
                st.warning("Switched off, a VCF is scored with no database "
                           "features at all. That is the weakest of the five "
                           "models: 0.90 PR-AUC with several cells and 0.78 "
                           "with one, against 0.99 when the databases are "
                           "there. An ANNOVAR table is unaffected.")

if not started and not uploaded:
    with middle:
        st.caption("Required. Click any of them.")
        with st.container(horizontal=True):
            for n, need in enumerate(REQUIREMENTS):
                with st.popover(need["name"], key=f"need-{n}"):
                    requirement(need)
        with st.container(key="needs"), st.expander(
                "What your files must be, and how to get there"):
            st.subheader("The five requirements in full")
            for need in REQUIREMENTS:
                with st.container(border=True):
                    requirement(need, with_name=True)

            st.subheader("Preparing the files")
            st.write("If you are starting from reads, this is the whole of it.")
            st.code(PREPARING, language="bash")
            st.write("Step 4 is optional. Without it VarTrace fetches dbSNP "
                     "and gnomAD from Ensembl, which was measured to be worth "
                     "the same. With it, REDIportal and COSMIC come along too, "
                     "which are worth very little more.")

            st.subheader("Annotation is optional")
            st.write("If your files carry no database columns, VarTrace "
                     "fetches dbSNP and gnomAD from Ensembl by itself. You do "
                     "not have to annotate anything.")
            st.write("If you have annotated, only ANNOVAR tables are read. "
                     "VEP, SnpEff and VarSeq write their annotation in a form "
                     "VarTrace does not parse. Such a VCF still works, but "
                     "its annotation is skipped and Ensembl is used instead.")
            st.write("A table from any annotator other than ANNOVAR cannot be "
                     "read at all. Export a VCF from it instead.")

            st.subheader("Which variants are scored")
            st.write("Only biallelic single base substitutions, which is what "
                     "the models were trained on. Indels and multiallelic "
                     "sites are counted and skipped, and the count is "
                     "reported, so you always know how many of your variants "
                     "were scored.")

    with aside:
        st.markdown("""
<div class="vt-step">
  <div class="vt-step-n">1</div>
  <div><div class="vt-step-t">Your files are read and checked</div>
  <div class="vt-step-d">On this machine. VarTrace reports what it found,
  down to the assembly and the missing fields, and refuses rather than score
  an input it was not trained for.</div></div>
</div>
<div class="vt-step">
  <div class="vt-step-n">2</div>
  <div><div class="vt-step-t">You confirm, then it scores</div>
  <div class="vt-step-d">Missing database features are fetched from Ensembl
  first, if you left that on. One of five models is chosen from what the input
  supports, and the choice is stated in words.</div></div>
</div>
<div class="vt-step">
  <div class="vt-step-n">3</div>
  <div><div class="vt-step-t">A probability for every variant</div>
  <div class="vt-step-d">With the substitution spectrum, the allele fraction
  against the score, a per-cell breakdown, links out to IGV, and a record of
  exactly what produced the numbers.</div></div>
</div>

""", unsafe_allow_html=True)

        for n, (column, stat) in enumerate(zip(st.columns(4), STATS)):
            with column:
                with st.popover(stat["number"], key=f"stat-{n}"):
                    statistic(stat)
                st.markdown(f'<div class="vt-stat-l">{stat["label"]}</div>',
                            unsafe_allow_html=True)

    st.stop()

if not started:
    files = tuple((f.name, f.getvalue()) for f in uploaded)

with st.spinner("Checking the files"):
    checked = inspect(files, use_ensembl)

if checked is None:
    st.error("No usable variants in those files.")
    st.stop()
if "refused" in checked:
    st.error(checked["refused"])
    st.stop()

if not started:
    st.subheader("What arrived")
what_was_read(checked["facts"])
which_model_and_why(checked["model_now"], checked["model_later"])


problems = checked["problems"]
anyway = st.session_state.get("anyway", False)

if problems and not anyway:
    st.subheader("These files are not what the models were trained on")
    st.caption("Scored as they are, the numbers would look perfectly ordinary "
               "and mean something else. Each problem below says what to do "
               "about it.")
    for problem in problems:
        with st.container(border=True):
            st.error(problem["what"])
            st.write("What to do: " + problem["do"])
    st.checkbox("I have read this and want the scores anyway", key="anyway")
    st.stop()

if not started:
    if problems:
        st.warning("Scoring against the above, at your own reading of it.")
    else:
        st.success("These files are what the models expect. Nothing is "
                   "missing and the assembly is right.")
    if st.button("Score these variants", type="primary"):
        st.session_state["files"] = files
        st.session_state["use_ensembl"] = use_ensembl
        st.session_state["agreed"] = True
        st.rerun()
    st.stop()

with st.spinner("Annotating and scoring"):
    result = compute(files, use_ensembl, anyway)

if result is None or "problems" in result:
    st.error("The files could not be scored.")
    st.stop()

with st.container(border=True):
    left, right = st.columns([2, 3])
    threshold = left.slider("Threshold", 0.0, 1.0, 0.5, 0.01)
    right.write(f"A variant scoring {threshold:.2f} or above is called absent "
                "from the genome, and anything below it is called present. "
                "Moving this changes only the calls, never the scores.")
    right.write("At 0.5 the two mistakes count the same. At 0.867, of the "
                "variants called absent about 99 in 100 really are, at the "
                "price of finding fewer of them.")

table = result["table"].copy()
table["call"] = np.where(table["score"] >= threshold, ABSENT, PRESENT)
for message in result["messages"]:
    st.info(message)

if result["annotation"]:
    a = result["annotation"]
    st.caption(f"dbSNP and gnomAD came from Ensembl release {a['release']}, "
               f"which had a record for {a['with_a_record']} of "
               f"{a['queried']} variants.")

overview, variants_tab, details, about = st.tabs(
    ["Overview", "Variants", "Details", "About this run"])

with overview:
    a, b, c, d = st.columns(4)
    a.metric("Variants", len(table))
    b.metric("Cells", result["n_cells"])
    c.metric(f"Absent from the genome at {threshold:g}",
             int((table["call"] == ABSENT).sum()))
    d.metric("Uncertain (0.1 to 0.9)",
             int(((table["score"] > 0.1) & (table["score"] < 0.9)).sum()))
    st.caption(result["why"])
    if result["skipped"]:
        st.caption(f"{result['skipped']} rows were not biallelic single base "
                   "substitutions and were skipped")

    left, right = st.columns(2)
    with left:
        st.subheader("Scores")
        st.altair_chart(charts.score_histogram(table), use_container_width=True)
    with right:
        st.subheader("Substitutions")
        st.altair_chart(charts.spectrum(table), use_container_width=True)
        st.caption("Each class as a share of itself, so the two are comparable "
                   "even when one is much larger.")

with variants_tab:
    f1, f2, f3 = st.columns([2, 2, 1])
    chroms = sorted(table["chrom"].unique(), key=charts.chromosome_order)
    keep_chroms = f1.multiselect("Chromosome", chroms, default=chroms)
    keep_call = f2.selectbox("Call", ["all", ABSENT, PRESENT])
    depths = table["dp"].dropna()
    min_depth = f3.number_input("Minimum depth", min_value=0,
                                value=0, step=1,
                                help=None if depths.empty else
                                f"depth runs from {int(depths.min())} to "
                                f"{int(depths.max())}")

    shown = table[table["chrom"].isin(keep_chroms)]
    if keep_call != "all":
        shown = shown[shown["call"] == keep_call]
    if min_depth:
        shown = shown[shown["dp"].fillna(0) >= min_depth]

    shown = shown.sort_values("score", ascending=False)
    shown["IGV"] = [f"https://igv.org/app/?locus={c}:{p}"
                    for c, p in zip(shown["chrom"], shown["pos"])]
    st.write(f"{len(shown)} of {len(table)} variants")
    st.dataframe(
        shown, hide_index=True, use_container_width=True,
        column_config={
            "score": st.column_config.ProgressColumn(
                "score", min_value=0.0, max_value=1.0, format="%.4f"),
            "IGV": st.column_config.LinkColumn("IGV", display_text="open"),
        })
    st.download_button("Download this table",
                       shown.drop(columns="IGV").to_csv(sep="\t", index=False),
                       file_name="vartrace_results.tsv")

with details:
    st.subheader("Allele fraction against score")
    st.altair_chart(charts.vaf_against_score(table, threshold),
                    use_container_width=True)
    st.caption("Errors concentrate near an allele fraction of 0.5, where a "
               "germline heterozygote and an RNA level event look alike.")

    st.subheader("By cell")
    per_cell = table.groupby("cell").agg(
        variants=("score", "size"),
        absent=("call", lambda s: int((s == ABSENT).sum())),
        median_score=("score", "median"),
    ).reset_index()
    per_cell["share absent"] = (per_cell["absent"] / per_cell["variants"]).round(3)
    st.dataframe(per_cell, hide_index=True, use_container_width=True)

with about:
    st.code(result["provenance"].replace("see below", str(threshold)))
    st.caption("Everything needed to say later exactly what produced these "
               "numbers.")

st.divider()
_, saving, ending = st.columns([3, 1, 1])

# the table as scored plus the record of what produced it, in one file, so the
# two cannot be separated later and leave a table nobody can account for
bundle = BytesIO()
with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as archive:
    archive.writestr("vartrace_results.tsv",
                     table.to_csv(sep="\t", index=False))
    archive.writestr("vartrace_results.info.txt",
                     result["provenance"].replace("see below", str(threshold)))

saving.download_button("Save the results", bundle.getvalue(),
                       file_name="vartrace_results.zip", mime="application/zip",
                       type="primary", use_container_width=True,
                       help="The scored table and the record of what produced "
                            "it, together")

if ending.button("End this run", type="primary", use_container_width=True,
                 help="Clear the scores and go back to the files"):
    for name in ("agreed", "anyway", "files", "use_ensembl", "many", "single"):
        st.session_state.pop(name, None)
    st.rerun()
