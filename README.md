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

One line, and it does not matter whether you have Python:

    curl -LsSf https://codeberg.org/KamenYovchev/vartrace/raw/branch/main/install.sh | sh

Then vartrace opens the interface in your browser. The installer uses
whatever you already have and fetches a Python of its own only if you have
none. Read it first if you would rather: it is a short file and it is in this
repository.

On Windows, run the same line inside WSL.

If you keep your own Python environments, the ordinary route works too:

    pipx install vartrace

or, for the latest state of the source:

    pipx install "git+https://codeberg.org/KamenYovchev/vartrace.git"

The trained models come with the package and the web interface is included, so
nothing else has to be downloaded or installed separately.


## Using it

Run it with nothing after it and the web interface opens in the browser. It
listens only on your own machine, so the files you load stay on it.

    vartrace

Give it files instead and it works from the command line, one file per cell or
one file holding many cells:

    vartrace cell_01.vcf
    vartrace cells/*.tsv --threshold 0.867 --output results.tsv

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


## Where the database features come from

Only two things reach the models: the columns ANNOVAR writes, or the answers
Ensembl gives. Nothing else is read, and that is deliberate.

A VCF is a standard, so the nine fixed columns survive whoever touched the file
last. An annotator adds its own fields to INFO and leaves the rest alone. VEP
writes CSQ, SnpEff writes ANN, VarSeq writes its own; VarTrace reads none of
them. Such a file is read as the plain VCF it still is, its annotation passed
over, and dbSNP and gnomAD fetched from Ensembl instead. The result is the same
as for a file that was never annotated at all.

A table is not a standard, and every tool invents its own columns. The one from
ANNOVAR is read because the VCF fields survive in its Otherinfo columns. A
table from any other annotator is refused with a message saying to export a VCF
instead, which is one step in all of them and keeps more, not less: most such
tables drop the read-level INFO fields and the sample column, so the variants
would fail the checks below anyway.

Fetching from Ensembl needs an internet connection and can be switched off. Of
the five database features, dbSNP and gnomAD carry nearly all of the weight,
which is why REDIportal and COSMIC are used when present but never fetched.

The two sources were compared on the labelled data behind the models, all 12,916
rows. Taking dbSNP and gnomAD from Ensembl rather than from ANNOVAR cost 0.0009
PR-AUC and changed the call at 0.5 for 0.87 percent of variants, even though the
gnomAD frequencies themselves differed by more than 0.1 for 464 of them. That
measurement is on the cell lines the models were trained on and says the two
sources are interchangeable; it does not say what either would score on data
from somewhere else.


## What is checked before anything is scored

VarTrace reads the files, says what it found, and stops rather than score an
input that would produce ordinary looking numbers meaning something else. It
checks the assembly from the length of chr1 in the header and refuses GRCh37.
It also checks whether the read-level INFO fields are there, whether there is an
allele fraction, and whether the depth reaches 5; these three refuse when more
than a fifth of the variants fail. Each refusal says what to do about it, and
--anyway overrides them all.

Four things leave no trace in a VCF and cannot be checked: whether the reads
came from single cells or from bulk, the sequencing chemistry, the tissue, and
the aligner. Those are stated in the instructions instead.


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


## Analysis pipeline

The scripts behind the thesis, from the raw reads to the trained models, are in
pipeline/. See pipeline/README.md for the order and the settings.


## Licence

Code: PolyForm Noncommercial License 1.0.0, see LICENSE.
Trained models: Creative Commons Attribution-NonCommercial 4.0, see
LICENSE-MODELS.

Free for research and any other non-commercial use. Commercial use is not
permitted.


## Citation

To follow.
