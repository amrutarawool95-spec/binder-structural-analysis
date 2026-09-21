"""Streamlit UI for the curated PDB structural analysis."""

from __future__ import annotations

import tempfile
import urllib.error
from pathlib import Path

import pandas as pd
import streamlit as st

from main import STRUCTURES, create_plot, summarize_target


st.set_page_config(
    page_title="Binder Structural Analysis",
    page_icon="🧬",
    layout="wide",
)


@st.cache_data(ttl=3600, show_spinner=False)
def load_results() -> tuple[pd.DataFrame, list[str]]:
    """Download and analyze the curated dataset, caching results for one hour."""

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="binder-streamlit-") as temp:
        work_dir = Path(temp)
        for target in STRUCTURES:
            try:
                rows.append(summarize_target(target, work_dir))
            except (OSError, ValueError, urllib.error.URLError) as error:
                errors.append(f"{target.pdb_id} chain {target.chain}: {error}")

    if not rows:
        raise RuntimeError(
            "No structures could be analyzed. Check the app's network access "
            "to files.rcsb.org and try refreshing the data."
        )
    return pd.DataFrame(rows), errors


def create_plot_bytes(data: pd.DataFrame) -> bytes:
    """Render the shared analysis plot without writing into the repository."""

    with tempfile.TemporaryDirectory(prefix="binder-plot-") as temp:
        plot_path = Path(temp) / "confidence_comparison.png"
        create_plot(data.to_dict("records"), plot_path)
        return plot_path.read_bytes()


st.title("Computational Post-Analysis of PDB Structures")
st.write(
    "Compare chain-level B-factors from curated de novo miniprotein and "
    "natural-antibody structures without running GPU-heavy models."
)

with st.sidebar:
    st.header("Controls")
    if st.button("Refresh data"):
        st.cache_data.clear()
        st.rerun()
    st.caption("Results are cached for one hour to avoid repeated PDB downloads.")

    st.header("Scientific scope")
    st.info(
        "These entries are experimental X-ray structures. Their B-factors are "
        "crystallographic temperature factors, not pLDDT."
    )

try:
    with st.spinner("Downloading and analyzing PDB structures…"):
        results, errors = load_results()
except Exception as error:
    st.error("The analysis could not load.")
    st.exception(error)
    st.stop()

if errors:
    st.warning("Some structures could not be analyzed:")
    for error in errors:
        st.write(f"- {error}")

summary = (
    results.groupby("group", sort=False)["mean_b_factor"]
    .agg(["count", "mean", "median"])
    .reset_index()
)

metric_columns = st.columns(len(summary))
for column, row in zip(metric_columns, summary.to_dict("records")):
    column.metric(
        label=f"{row['group']} mean B-factor",
        value=f"{row['mean']:.2f}",
        help="Source-specific experimental B-factor units; not pLDDT.",
    )

st.subheader("Comparison")
st.image(
    create_plot_bytes(results),
    caption="Mean atom B-factor with one point per analyzed PDB chain.",
    width="stretch",
)

st.subheader("Structure-level results")
display_columns = [
    "pdb_id",
    "group",
    "chain",
    "description",
    "atom_count",
    "mean_b_factor",
    "median_b_factor",
    "stdev_b_factor",
]
st.dataframe(
    results[display_columns].round({"mean_b_factor": 2, "median_b_factor": 2, "stdev_b_factor": 2}),
    width="stretch",
    hide_index=True,
)
st.download_button(
    "Download results CSV",
    data=results.to_csv(index=False),
    file_name="structural_metrics.csv",
    mime="text/csv",
)

with st.expander("Dataset and interpretation"):
    st.dataframe(
        pd.DataFrame(
            [
                {
                    "PDB": target.pdb_id,
                    "Group": target.group,
                    "Chain": target.chain,
                    "Structure": target.description,
                }
                for target in STRUCTURES
            ]
        ),
        width="stretch",
        hide_index=True,
    )
    st.write(
        "A B-factor can represent pLDDT only when the source prediction "
        "pipeline explicitly documents that convention. Do not compare these "
        "experimental temperature factors as AI confidence scores."
    )