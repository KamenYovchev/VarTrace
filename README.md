# VarTrace

VarTrace estimates, for each variant observed in single-cell RNA sequencing, the
probability that the variant is not present in the genome of that cell. It works
from the RNA data alone, without paired DNA sequencing.


## Why

A variant seen in RNA may be present in the genome, or it may arise only at the
RNA level. Telling the two apart normally requires sequencing the DNA of the same
cell, which is expensive and often not possible. The common shortcut, checking
whether the variant is listed in public databases, misclassifies a large share of
cases in both directions.

VarTrace answers the question from the RNA data and returns a calibrated
probability rather than a yes or no. A score of 0.87 means that of a hundred
variants scored that way, about 87 really are absent from the genome.


## Installing

    pipx install "git+https://codeberg.org/KamenYovchev/vartrace.git"

For the web interface as well:

    pipx install "vartrace[web] @ git+https://codeberg.org/KamenYovchev/vartrace.git"

The trained models come with the package, so nothing else has to be downloaded.


## Using it

One file per cell, or one file holding many cells:

    vartrace cell_01.vcf
    vartrace cells/*.tsv --threshold 0.867 --output results.tsv

A web interface, which opens in the browser and listens only on your own
machine:

    vartrace serve

Options:

    --threshold     0.5 by default. 0.867 is the high confidence mode, where
                    precision is about 0.99 at the cost of some recall.
    --no-annotate   do not ask Ensembl for dbSNP and gnomAD when the input
                    carries no annotation of its own.
    --output        where to write the table. A second file with the same name
                    and .info.txt records what produced it.


## Input

A plain VCF, as it comes from the variant caller, or a table annotated with
ANNOVAR, where the VCF fields sit in the Otherinfo columns.

A VCF may hold one sample or many. One column per cell barcode is read as many
cells, and a cell is counted for a variant only when its genotype or its allele
depth says it carries it.

Only biallelic single base substitutions are used, because that is what the
models were trained on. Anything else is counted and skipped, and the count is
reported.

If the input carries no database annotation, VarTrace looks up dbSNP and gnomAD
through the public Ensembl service. That needs an internet connection and can be
switched off. Of the five database features these two carry nearly all of the
weight, which is why the other two are used when present but never fetched.


## Output

A table with one row per variant: the cell, the position, the alleles, the
score, the call at the chosen threshold, and which model was used.

Alongside it a short record of what produced those numbers, including the model
file and its fingerprint, and where the annotation came from.

The tool also states in words why it chose the model it did, for example that
the input had dbSNP and gnomAD but not REDIportal or COSMIC.


## The models

Five models come with the package and one is chosen automatically from what the
input supports. No database annotation at all, or all four databases, is served
by a model trained for exactly that case, with or without the cross-cell
features depending on whether more than one cell was given. Anything in between
goes to a model trained with feature groups withheld at random, which is the
only one that can take a partial input.

Each model carries a description in JSON next to it: the algorithm, its
settings, the features it expects in order, what it was trained on, and its
validation metrics.


## Scope and limitations

The models were trained on 45 single cells from three cell lines, sequenced with
one technology, with paired DNA sequencing as the ground truth. They have not
been validated on primary tissue or on other protocols.

The score answers one question: is the variant present in the genome. It does
not say which biological process produced a variant that is absent from it.

The models were trained on gnomAD as ANNOVAR supplies it. When VarTrace fetches
the frequencies from Ensembl instead, the release may differ; which one answered
is recorded in the output.


## Licence

Code: PolyForm Noncommercial License 1.0.0, see LICENSE.
Trained models: Creative Commons Attribution-NonCommercial 4.0, see
LICENSE-MODELS.

Free for research and any other non-commercial use. Commercial use is not
permitted.


## Citation

To follow.
