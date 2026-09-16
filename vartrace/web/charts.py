"""The plots on the page, kept apart from it so they can be checked alone."""

import altair as alt
import numpy as np
import pandas as pd


def chromosome_order(name):
    body = name[3:] if name.startswith("chr") else name
    return (0, int(body)) if body.isdigit() else (1, body)


def spectrum(table):
    """The twelve substitution types, as a share of each class."""
    counts = table.groupby(["substitution", "call"]).size().reset_index(name="n")
    totals = counts.groupby("call")["n"].transform("sum")
    counts["percent"] = 100 * counts["n"] / totals
    return alt.Chart(counts).mark_bar().encode(
        x=alt.X("substitution:N", title=None),
        y=alt.Y("percent:Q", title="percent of class"),
        color=alt.Color("call:N", title=None),
        xOffset="call:N",
        tooltip=["substitution", "call", "n",
                 alt.Tooltip("percent:Q", format=".1f")],
    ).properties(height=260)


def vaf_against_score(table, threshold):
    """Allele fraction against score, with the threshold drawn across it."""
    points = alt.Chart(table).mark_circle(size=28, opacity=0.5).encode(
        x=alt.X("vaf:Q", title="allele fraction", scale=alt.Scale(domain=[0, 1])),
        y=alt.Y("score:Q", title="score", scale=alt.Scale(domain=[0, 1])),
        color=alt.Color("call:N", title=None),
        tooltip=["cell", "chrom", "pos", "substitution", "vaf", "score"],
    )
    line = alt.Chart(pd.DataFrame({"y": [threshold]})).mark_rule(
        strokeDash=[4, 4]).encode(y="y:Q")
    return (points + line).properties(height=320)


def score_histogram(table):
    counts, edges = np.histogram(table["score"], bins=20, range=(0, 1))
    data = pd.DataFrame({"score": np.round(edges[:-1] + 0.025, 3),
                         "variants": counts})
    return alt.Chart(data).mark_bar().encode(
        x=alt.X("score:Q", title="score"),
        y=alt.Y("variants:Q", title="variants"),
    ).properties(height=260)
