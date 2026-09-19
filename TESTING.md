# Trying VarTrace on your own data

Thank you for agreeing to test this. It is part of a master's thesis at the
Faculty of Mathematics and Informatics, Sofia University, and honest criticism
is worth more to it than polite approval. If it refuses your files, says
something unclear, or gives you a number you do not believe, that is exactly
what I need to hear.


## What it does

For each variant seen in single-cell RNA sequencing, VarTrace estimates the
probability that the variant is not present in the genome of that cell. It
works from the RNA data alone, without paired DNA sequencing.

A score of 0.87 means that of a hundred variants scored that way, about 87
really are absent from the genome. It is a calibrated probability, not a
ranking.


## Installing and running

Paste this into a terminal:

    curl -LsSf https://codeberg.org/KamenYovchev/vartrace/raw/branch/main/install.sh | sh

On Windows, do it inside WSL. If that is a problem for you, say so in your
answers below.

Then type:

    vartrace

and the interface opens in your browser. You do not need Python for any of
this; the installer brings one if you have none. It is a short script and it
is in the repository, so read it first if you would rather.

Everything runs on your own machine and the files you load never leave it.

From the command line instead:

    vartrace cell_01.vcf
    vartrace cells/*.vcf --output results.tsv


## What your data has to be

The models were trained on 45 single cells from three cell lines, with paired
DNA sequencing as the truth. The closer your data sit to that, the more the
numbers mean.

Required:

    human variants on GRCh38
    from single-cell RNA sequencing, not bulk RNA and not DNA
    a VCF, or a table from ANNOVAR
    AD or VAF in the sample column, so there is an allele fraction
    coverage of at least 5 reads

A VCF from any caller and any annotator is fine. Tables from annotators other
than ANNOVAR are not read, so export a VCF from those instead. Everything below
5 reads was removed from the training data, so remove it from yours.

Strongly wanted: called with GATK, or with something that writes the same INFO
fields, FS, SOR, MQ, MQRankSum, ReadPosRankSum, BaseQRankSum and QD. Seven of
the thirty-five features come from there, and unlike the database features, no
model was trained without them.

Not needed: database annotation. If your files carry none, VarTrace looks up
dbSNP and gnomAD at Ensembl by itself.

Annotation from anything other than ANNOVAR is passed over rather than read. A
VCF from VEP, SnpEff or VarSeq is treated as the plain VCF it still is, and the
databases come from Ensembl. So you lose nothing by annotating, and gain
nothing either, unless you used ANNOVAR.

The two sources were compared on the data behind the models: taking dbSNP and
gnomAD from Ensembl rather than from ANNOVAR cost 0.0009 PR-AUC and changed the
call for 0.87 percent of variants. That was measured on the cell lines the
models were trained on, and says the two sources are interchangeable. It does
not say what either would score on your data.


## Which variants are scored

One base swapped for one other base, with a single alternative at that
position:

    REF   ALT
    G     C      scored
    G     C,T    skipped, two alternatives at one position
    GA    G      skipped, a deletion
    G     GA     skipped, an insertion

The training data were filtered to biallelic SNVs with bcftools view -v snps
-m2 -M2, so the models have seen nothing else. The question itself is about
substitutions: ADAR editing turns A into I, which sequences as G, and reverse
transcriptase errors are substitutions too. Indels arise differently.

Skipped rows are counted and the count is reported, so you always know how many
of the variants you submitted were actually scored.


## What VarTrace checks, and what it cannot

Before scoring anything it reads your files, says what it found, and stops if
something would make the numbers mean something else. It checks the assembly,
the read-level annotations, the allele fractions and the depth, and for each
problem it says what to do about it.

Four things it cannot check, and which are down to you:

    whether the reads came from single cells or from bulk
    which sequencing chemistry was used
    what tissue or cell type it was
    which aligner produced the BAM files

Nothing in a VCF records any of these. If your data differ from the training
data in one of them, the scores may be off in ways neither of us can see, and
that is worth saying in your answer below.


## Where it has not been tested

Primary tissue. The training data are cell lines, which are clean and uniform
in a way tissue is not.

Solid tumours. Two of the three lines are healthy reference genomes and the
third is a leukaemia line.

Other sequencing technologies. All 45 cells came from one.

Other callers. See the note about the read-level features above.

If your data fall into one of these, that is the most useful test of all. It is
also the one most likely to produce a wrong answer, so please do not take the
numbers on trust.


## What I would like to hear back

Short answers are fine. Anything you write will be described in the thesis,
anonymously unless you say otherwise.

1. Who you are and what you work on, in a sentence.

2. What you ran it on: organism, assembly, technology, variant caller, how
   many cells, roughly how many variants.

3. Did it install and start? How long did that take?

4. Did it refuse your files? If so, what did it say, and could you tell from
   the message what to do?

5. Do the scores look plausible? If you have any variants whose origin you
   already know, did VarTrace agree?

6. What was confusing or badly explained, in the interface or in the output?

7. What is missing before you would use this in real work?

8. Would you run it again on your next dataset? Why, or why not?

Send answers to kamen.yovchev@gmail.com.


## If something breaks

Please send the exact message, what you ran, and roughly what the input looked
like. You do not need to send the data itself.

Source: https://codeberg.org/KamenYovchev/vartrace
