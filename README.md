# VarTrace

VarTrace estimates, for each variant observed in single-cell RNA sequencing, the
probability that the variant is not present in the genome of that cell. It works
from the RNA data alone, without paired DNA sequencing.

Status: under development. The tool is not yet ready for use.


## Why

A variant seen in RNA may be present in the genome, or it may arise only at the RNA
level. Telling the two apart normally requires sequencing the DNA of the same cell,
which is expensive and often not possible. The common shortcut, checking whether the
variant is listed in public databases, misclassifies a large share of cases in both
directions.

VarTrace answers the question from the RNA data and returns a calibrated
probability rather than a yes or no.


## What it does

Input: variants already called from single-cell RNA sequencing, either as a plain
VCF or as a table annotated with ANNOVAR. One file per cell, one or several cells.

Output: a table with one row per variant, a score between 0 and 1, the class at the
chosen threshold, and a short statement of which model was used and why.

If the input carries no database annotation, VarTrace can look up dbSNP and gnomAD
through the public Ensembl service. This needs an internet connection and can be
switched off.

The score answers one question: is the variant present in the genome. It does not
say which biological process produced a variant that is absent from the genome.


## Installation

Not yet available. Planned:

    pip install git+https://codeberg.org/KamenYovchev/vartrace.git

or, to keep it in its own environment:

    pipx install git+https://codeberg.org/KamenYovchev/vartrace.git

Releases on PyPI and Bioconda will follow.


## Usage

Planned:

    vartrace cell_01.vcf
    vartrace cells/*.vcf --threshold 0.9 --output results.tsv


## Scope and limitations

The models were trained on single cells from three cell lines, sequenced with one
technology. They have not been validated on primary tissue or on other protocols.


## Licence

Code: PolyForm Noncommercial License 1.0.0, see LICENSE.
Trained models: Creative Commons Attribution-NonCommercial 4.0, see LICENSE-MODELS.

Free for research and any other non-commercial use. Commercial use is not permitted.


## Citation

To follow.
