#!/usr/bin/env python3
"""Run again the Europe PMC queries used to build the literature corpus.

The queries, page counts and sort order are copied from the original search
(2026-08-18). The citation counts change with time and the results are
sorted by them, so the returned set is close to the original but not identical.

Output: raw.json, one list of records per query group.
"""
import json, time, urllib.parse, urllib.request

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def get(params, tries=4):
    url = EPMC + "?" + urllib.parse.urlencode(params)
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception:
            time.sleep(3 * (t + 1))
    return None


def epmc(query, pages, pagesize):
    out, cursor = [], "*"
    for _ in range(pages):
        j = get({"query": query, "format": "json", "pageSize": pagesize, "resultType": "core",
                 "cursorMark": cursor, "sort": "CITED desc"})
        if j is None:
            print("FAILED", query[:80])
            break
        res = j.get("resultList", {}).get("result", [])
        out += res
        nc = j.get("nextCursorMark")
        if not nc or nc == cursor or len(res) < pagesize:
            break
        cursor = nc
        time.sleep(0.5)
    return out


themes = {
 "scRNA_variant_calling": '(TITLE_ABS:"single-cell RNA" OR TITLE_ABS:"scRNA-seq") AND (TITLE_ABS:"variant calling" OR TITLE_ABS:"somatic mutation" OR TITLE_ABS:"somatic variant" OR TITLE_ABS:"genotyping") AND FIRST_PDATE:[2015 TO 2026]',
 "rna_editing": '(TITLE_ABS:"RNA editing" OR TITLE_ABS:"A-to-I editing" OR TITLE_ABS:"ADAR") AND (TITLE_ABS:"single cell" OR TITLE_ABS:"single-cell" OR TITLE_ABS:"transcriptome-wide" OR TITLE_ABS:"detection" OR TITLE_ABS:"machine learning") AND FIRST_PDATE:[2015 TO 2026]',
 "paired_dna_rna_singlecell": '(TITLE_ABS:"single-cell" ) AND (TITLE_ABS:"paired" OR TITLE_ABS:"joint" OR TITLE_ABS:"simultaneous" OR TITLE_ABS:"multi-omic" OR TITLE_ABS:"multiomic") AND (TITLE_ABS:"genome and transcriptome" OR TITLE_ABS:"DNA and RNA" OR TITLE_ABS:"scWGS" OR TITLE_ABS:"whole-genome sequencing") AND FIRST_PDATE:[2014 TO 2026]',
 "allelic_dropout_artifacts": '(TITLE_ABS:"allelic dropout" OR TITLE_ABS:"allele dropout" OR TITLE_ABS:"amplification artifact" OR TITLE_ABS:"allele-specific expression" OR TITLE_ABS:"allelic imbalance") AND (TITLE_ABS:"single cell" OR TITLE_ABS:"single-cell") AND FIRST_PDATE:[2014 TO 2026]',
 "noisy_labels_ml": '(TITLE_ABS:"label noise" OR TITLE_ABS:"noisy labels" OR TITLE_ABS:"positive-unlabeled" OR TITLE_ABS:"weak supervision" OR TITLE_ABS:"uncertain labels" OR TITLE_ABS:"calibration") AND (TITLE_ABS:"genomics" OR TITLE_ABS:"variant" OR TITLE_ABS:"single-cell" OR TITLE_ABS:"bioinformatics" OR TITLE_ABS:"machine learning") AND FIRST_PDATE:[2016 TO 2026]',
 "mutational_signatures": '(TITLE_ABS:"mutational signature" OR TITLE_ABS:"mutational signatures" OR TITLE_ABS:"SigProfiler" OR TITLE_ABS:"COSMIC signatures") AND FIRST_PDATE:[2015 TO 2026]',
 "transcription_errors": '(TITLE_ABS:"transcription error" OR TITLE_ABS:"transcriptional error" OR TITLE_ABS:"transcription fidelity" OR TITLE_ABS:"RNA polymerase fidelity" OR TITLE_ABS:"transcriptional mutagenesis") AND FIRST_PDATE:[2010 TO 2026]',
 "neoantigen_expressed": '(TITLE_ABS:"neoantigen" OR TITLE_ABS:"neoepitope") AND (TITLE_ABS:"expressed" OR TITLE_ABS:"RNA-seq" OR TITLE_ABS:"single-cell" OR TITLE_ABS:"transcript") AND FIRST_PDATE:[2016 TO 2026]',
 "ml_variant_filtering": '(TITLE_ABS:"deep learning" OR TITLE_ABS:"machine learning" OR TITLE_ABS:"random forest" OR TITLE_ABS:"gradient boosting" OR TITLE_ABS:"neural network") AND (TITLE_ABS:"variant calling" OR TITLE_ABS:"somatic variant" OR TITLE_ABS:"false positive variant" OR TITLE_ABS:"variant filtering") AND FIRST_PDATE:[2016 TO 2026]',
 "strand_bias_transcribed": '(TITLE_ABS:"transcriptional strand bias" OR TITLE_ABS:"transcription-coupled repair" OR TITLE_ABS:"strand-specific mutation" OR TITLE_ABS:"strand asymmetry") AND FIRST_PDATE:[2012 TO 2026]',
 "clonal_hematopoiesis_singlecell": '(TITLE_ABS:"acute myeloid leukemia" OR TITLE_ABS:"clonal hematopoiesis" OR TITLE_ABS:"leukemia") AND (TITLE_ABS:"single-cell" ) AND (TITLE_ABS:"mutation" OR TITLE_ABS:"genotyping" OR TITLE_ABS:"clone") AND FIRST_PDATE:[2017 TO 2026]',
 "somatic_mosaicism_normal": '(TITLE_ABS:"somatic mutation" OR TITLE_ABS:"somatic mosaicism") AND (TITLE_ABS:"single cell" OR TITLE_ABS:"single-cell") AND (TITLE_ABS:"normal" OR TITLE_ABS:"neuron" OR TITLE_ABS:"aging" OR TITLE_ABS:"whole-genome") AND FIRST_PDATE:[2015 TO 2026]',
}

tools_named = [
 'TITLE_ABS:"SComatic"', 'TITLE_ABS:"Monopogen"', 'TITLE_ABS:"cellsnp-lite"', 'TITLE_ABS:"Vireo"',
 'TITLE_ABS:"scAllele"', 'TITLE_ABS:"RNA-MuTect"', 'TITLE_ABS:"SNPiR"', 'TITLE_ABS:"RESA" AND TITLE_ABS:"single-cell"',
 'TITLE_ABS:"Trinity CTAT"', 'TITLE_ABS:"SCmut"', 'TITLE_ABS:"RedNet" OR TITLE_ABS:"RED-ML" OR TITLE_ABS:"DeepRed"',
 'TITLE_ABS:"REDIportal"', 'TITLE_ABS:"MutaSeq" OR TITLE_ABS:"GoT" AND TITLE_ABS:"genotyping of transcriptomes"',
 'TITLE_ABS:"DeepVariant"', 'TITLE_ABS:"Mutect2" OR TITLE_ABS:"MuTect"', 'TITLE_ABS:"SigProfiler"',
 'TITLE_ABS:"scSNV" OR TITLE_ABS:"SCReadCounts"', 'TITLE_ABS:"XClone" OR TITLE_ABS:"numbat" AND TITLE_ABS:"single-cell"',
 'TITLE_ABS:"MAESTER"', 'TITLE_ABS:"SIGURD" OR TITLE_ABS:"mitochondrial lineage tracing"',
 'TITLE_ABS:"scMutRates" OR TITLE_ABS:"SCAN2"', 'TITLE_ABS:"LiRA" AND TITLE_ABS:"somatic"',
 'TITLE_ABS:"Parabricks"', 'TITLE_ABS:"scDNA-seq" AND TITLE_ABS:"variant"',
]

targeted = [
 'TITLE_ABS:"RNA-MuTect" OR (TITLE:"somatic mutations" AND TITLE:"RNA" AND AUTH:"Yizhak")',
 '(TITLE_ABS:"G&T-seq") OR (TITLE_ABS:"genome and transcriptome sequencing" AND TITLE_ABS:"single-cell")',
 '(TITLE_ABS:"genotyping of transcriptomes") OR (TITLE_ABS:"somatic mutations" AND TITLE_ABS:"clonal" AND TITLE_ABS:"single-cell RNA")',
 'TITLE_ABS:"LiRA" OR (TITLE_ABS:"linked-read" AND TITLE_ABS:"somatic" AND TITLE_ABS:"single cell")',
 'TITLE_ABS:"scAllele" OR TITLE_ABS:"allele-specific" AND TITLE_ABS:"variant" AND TITLE_ABS:"RNA-seq"',
 'TITLE_ABS:"REDItools" OR TITLE_ABS:"SPRINT" AND TITLE_ABS:"RNA editing" OR TITLE_ABS:"JACUSA"',
 '(TITLE_ABS:"RNA editing" AND TITLE_ABS:"single-cell") AND FIRST_PDATE:[2019 TO 2026]',
 '(TITLE_ABS:"neoantigen" AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"RNA editing" OR TITLE_ABS:"expression")) AND FIRST_PDATE:[2018 TO 2026]',
 '(TITLE_ABS:"mutational signature" AND (TITLE_ABS:"single-cell" OR TITLE_ABS:"strand" OR TITLE_ABS:"transcribed strand" OR TITLE_ABS:"RNA")) AND FIRST_PDATE:[2016 TO 2026]',
 '(TITLE_ABS:"single-cell whole-genome" OR TITLE_ABS:"single cell whole genome") AND (TITLE_ABS:"somatic" OR TITLE_ABS:"mutation" OR TITLE_ABS:"amplification")',
 '(TITLE_ABS:"RNA-DNA difference" OR TITLE_ABS:"RNA-DNA differences" OR TITLE_ABS:"RDD")',
 '(TITLE_ABS:"calibrat" AND TITLE_ABS:"probab") AND (TITLE_ABS:"variant" OR TITLE_ABS:"genomic" OR TITLE_ABS:"single-cell")',
 '(TITLE_ABS:"mitochondrial" AND TITLE_ABS:"single-cell" AND (TITLE_ABS:"lineage" OR TITLE_ABS:"heteroplasmy" OR TITLE_ABS:"variant"))',
 '(TITLE_ABS:"expressed" AND TITLE_ABS:"somatic mutations" AND (TITLE_ABS:"tumor" OR TITLE_ABS:"cancer"))',
]

must = [
 'TITLE:"RNA sequence analysis reveals macroscopic somatic clonal expansion across normal tissues"',
 'TITLE:"Genome-wide somatic variant calling using localized colocalization" OR TITLE_ABS:"SCAN2"',
 'TITLE:"Accurate identification of single-nucleotide variants in whole-genome-amplified single cells"',
 'TITLE_ABS:"cellsnp-lite"',
 'TITLE_ABS:"Vireo" AND TITLE_ABS:"donor deconvolution"',
 'TITLE:"Aging and neurodegeneration are associated with increased mutations in single human neurons"',
 'TITLE_ABS:"Somatic mutation landscapes at single-molecule resolution" OR TITLE_ABS:"NanoSeq"',
 'TITLE_ABS:"duplex sequencing" AND TITLE_ABS:"somatic mutation"',
 'TITLE:"Deciphering signatures of mutational processes operative in human cancer"',
 'TITLE_ABS:"transcription-coupled repair" AND TITLE_ABS:"mutation" AND FIRST_PDATE:[2018 TO 2026]',
 'TITLE_ABS:"A-to-I RNA editing" AND TITLE_ABS:"cancer" AND FIRST_PDATE:[2015 TO 2026]',
 'TITLE_ABS:"RNA editing" AND TITLE_ABS:"neoantigen"',
 'TITLE_ABS:"variant allele frequency" AND TITLE_ABS:"single-cell RNA" AND TITLE_ABS:"detection"',
 'TITLE_ABS:"error-corrected" AND TITLE_ABS:"RNA sequencing" AND TITLE_ABS:"transcription errors"',
]

raw = {}
for k, q in themes.items():
    raw[k] = epmc(q, pages=2, pagesize=100)
    print(k, len(raw[k]))
raw["tools_named"] = [r for q in tools_named for r in epmc(q, 1, 50)]
raw["targeted"] = [r for q in targeted for r in epmc(q, 1, 100)]
raw["must"] = [r for q in must for r in epmc(q, 1, 40)]
for k in ("tools_named", "targeted", "must"):
    print(k, len(raw[k]))
json.dump(raw, open("raw.json", "w"))
