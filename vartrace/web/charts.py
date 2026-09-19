"""The plots on the page, kept apart from it so they can be checked alone."""

import altair as alt
import numpy as np
import pandas as pd

# The two calls keep the same colour in every plot, so the eye does not have to
# relearn them from one chart to the next. They live here rather than in the
# page because the plots are what give them meaning.
ABSENT = "absent from the genome"
PRESENT = "in the genome"
GREEN = "#4E9A63"
BLUE = "#6E93B8"


def by_call():
    return alt.Color("call:N", title=None,
                     scale=alt.Scale(domain=[ABSENT, PRESENT],
                                     range=[GREEN, BLUE]))


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
        color=by_call(),
        xOffset="call:N",
        tooltip=["substitution", "call", "n",
                 alt.Tooltip("percent:Q", format=".1f")],
    ).properties(height=260)


def vaf_against_score(table, threshold):
    """Allele fraction against score, with the threshold drawn across it."""
    points = alt.Chart(table).mark_circle(size=28, opacity=0.5).encode(
        x=alt.X("vaf:Q", title="allele fraction", scale=alt.Scale(domain=[0, 1])),
        y=alt.Y("score:Q", title="score", scale=alt.Scale(domain=[0, 1])),
        color=by_call(),
        tooltip=["cell", "chrom", "pos", "substitution", "vaf", "score"],
    )
    line = alt.Chart(pd.DataFrame({"y": [threshold]})).mark_rule(
        strokeDash=[4, 4]).encode(y="y:Q")
    return (points + line).properties(height=320)


def score_histogram(table):
    counts, edges = np.histogram(table["score"], bins=20, range=(0, 1))
    data = pd.DataFrame({"score": np.round(edges[:-1] + 0.025, 3),
                         "variants": counts})
    return alt.Chart(data).mark_bar(color=GREEN).encode(
        x=alt.X("score:Q", title="score"),
        y=alt.Y("variants:Q", title="variants"),
    ).properties(height=260)
