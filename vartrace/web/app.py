"""The web interface. Started with: vartrace serve

Uploaded files are written to a temporary directory, read, and deleted as soon
as the scores are computed. Nothing is kept on the server.
"""

import os
import tempfile

import numpy as np
import pandas as pd
import streamlit as st

from vartrace import __version__, cli, io, model as models, report
from vartrace.web import charts

ABSENT = "absent from the genome"
PRESENT = "in the genome"


@st.cache_data(show_spinner=False)
def compute(files, use_ensembl):
    """Files are (name, bytes). Cached so moving the slider does not rerun it."""
    tmp = tempfile.mkdtemp(prefix="vartrace_")
    paths = []
    try:
        for name, data in files:
            path = os.path.join(tmp, os.path.basename(name))
            with open(path, "wb") as fh:
                fh.write(data)
            paths.append(path)

        variants, n_cells, skipped = io.read_all(paths)
        if not variants:
            return None
        messages = []
        model, why, scores, annotation = cli.score(
            variants, n_cells, models.MODELS_DIR, use_ensembl,
            log=messages.append)

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
        }
    finally:
        for path in paths:
            os.remove(path)
        os.rmdir(tmp)


st.set_page_config(page_title="VarTrace", layout="wide")

with st.sidebar:
    st.header("VarTrace")
    st.caption("The probability that a variant seen in single-cell RNA "
               "sequencing is not present in the genome of that cell.")
    uploaded = st.file_uploader("One VCF or ANNOVAR table per cell",
                                type=["vcf", "tsv", "txt"],
                                accept_multiple_files=True)
    threshold = st.slider("Threshold", 0.0, 1.0, 0.5, 0.001)
    st.caption("0.867 is the high confidence mode, where precision is about "
               "0.99")
    use_ensembl = st.checkbox("Look up dbSNP and gnomAD at Ensembl if the "
                              "input has no annotation", value=True)
    st.caption(f"vartrace {__version__}. Uploaded files are deleted as soon as "
               "the scores are computed.")

if not uploaded:
    st.info("Upload one file per cell to begin. A plain VCF is enough; if it "
            "carries no database annotation, VarTrace can fetch dbSNP and "
            "gnomAD itself.")
    st.stop()

with st.spinner("Reading, annotating and scoring"):
    result = compute(tuple((f.name, f.getvalue()) for f in uploaded), use_ensembl)

if result is None:
    st.error("No usable variants in those files.")
    st.stop()

table = result["table"].copy()
table["call"] = np.where(table["score"] >= threshold, ABSENT, PRESENT)
for message in result["messages"]:
    st.sidebar.info(message)

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
