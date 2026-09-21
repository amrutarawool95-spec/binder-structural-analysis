"""Streamlit UI for experimental B-factor and AI pLDDT analysis."""

from __future__ import annotations

import tempfile
import urllib.error
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from main import (
    STRUCTURES,
    create_plot,
    download_pdb,
    extract_b_factors,
    summarize_target,
)


LOCAL_RESULTS = Path(__file__).with_name("structural_metrics.csv")
EXPERIMENTAL_MODE = "Experimental structures — B-factor"
AI_MODE = "AI-predicted structures — pLDDT"


@st.cache_data
def load_bundled_results() -> tuple[pd.DataFrame, list[str]]:
    """Load the checked-in snapshot so the app can render without network I/O."""

    if not LOCAL_RESULTS.exists():
        raise RuntimeError(
            "The bundled structural_metrics.csv file is missing from the repository."
        )
    return pd.read_csv(LOCAL_RESULTS), []


@st.cache_data(ttl=3600, show_spinner=False)
def load_remote_results() -> tuple[pd.DataFrame, list[str]]:
    """Download and analyze the curated dataset, caching live results for one hour."""

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


def load_results(
    refresh_from_rcsb: bool,
) -> tuple[pd.DataFrame, list[str], str]:
    """Prefer the bundled snapshot and make live refresh opt-in and resilient."""

    if not refresh_from_rcsb:
        results, errors = load_bundled_results()
        return results, errors, "Showing the bundled verified snapshot."

    try:
        results, errors = load_remote_results()
        return results, errors, "Showing fresh results downloaded from RCSB."
    except Exception as error:
        results, errors = load_bundled_results()
        errors.append(
            f"Live refresh failed ({error}); showing the bundled snapshot instead."
        )
        return results, errors, "Live refresh failed safely; using the bundled snapshot."


def create_plot_bytes(data: pd.DataFrame) -> bytes:
    """Render the experimental analysis plot without writing to the repository."""

    with tempfile.TemporaryDirectory(prefix="binder-plot-") as temp:
        plot_path = Path(temp) / "confidence_comparison.png"
        create_plot(data.to_dict("records"), plot_path)
        return plot_path.read_bytes()


def analyze_ai_source(
    pdb_id: str,
    chain: str,
    uploaded_file: object | None,
) -> dict[str, object]:
    """Read one AI-predicted PDB source and interpret B-factors as pLDDT."""

    chain = chain.strip()
    if len(chain) != 1:
        raise ValueError("Enter exactly one chain identifier, such as A or B.")

    with tempfile.TemporaryDirectory(prefix="binder-ai-") as temp:
        work_dir = Path(temp)
        if uploaded_file is not None:
            filename = getattr(uploaded_file, "name", "uploaded_structure.pdb")
            pdb_path = work_dir / Path(filename).name
            pdb_path.write_bytes(uploaded_file.getvalue())
            source = f"Uploaded file: {filename}"
        else:
            pdb_id = pdb_id.strip().upper()
            if len(pdb_id) != 4 or not pdb_id.isalnum():
                raise ValueError("Enter a four-character PDB ID or upload a PDB file.")
            pdb_path = work_dir / f"{pdb_id}.pdb"
            download_pdb(pdb_id, pdb_path)
            source = f"RCSB PDB: {pdb_id}"

        values = extract_b_factors(pdb_path, chain)

    if not values:
        raise ValueError(
            f"No ATOM B-factors found for chain {chain}. "
            "Check the chain ID and PDB format."
        )
    outside_plddt_range = [value for value in values if not 0 <= value <= 100]
    if outside_plddt_range:
        raise ValueError(
            "This chain contains B-factor values outside the pLDDT range of "
            "0–100. It may be an experimental structure, or its prediction "
            "pipeline may use a different confidence convention."
        )

    series = pd.Series(values, dtype="float64")
    return {
        "source": source,
        "chain": chain,
        "atom_count": len(values),
        "values": values,
        "mean": float(series.mean()),
        "median": float(series.median()),
        "high_confidence_percent": float((series >= 90).mean() * 100),
        "low_confidence_percent": float((series < 50).mean() * 100),
    }


def create_plddt_plot(values: list[float]) -> plt.Figure:
    """Create a pLDDT distribution with standard confidence bands."""

    figure, axis = plt.subplots(figsize=(9, 5))
    axis.hist(
        values,
        bins=[0, 50, 70, 90, 100],
        color="#4C78A8",
        edgecolor="white",
        linewidth=1,
    )
    axis.axvline(50, color="#E45756", linestyle="--", label="50: low-confidence cutoff")
    axis.axvline(70, color="#F2CF5B", linestyle="--", label="70: good-confidence cutoff")
    axis.axvline(90, color="#54A24B", linestyle="--", label="90: high-confidence cutoff")
    axis.set_xlim(0, 100)
    axis.set_xlabel("pLDDT confidence score")
    axis.set_ylabel("Atom count")
    axis.set_title("AI structure pLDDT distribution", weight="bold")
    axis.legend(loc="upper left", fontsize="small")
    figure.tight_layout()
    return figure


def render_experimental_mode() -> None:
    """Render the crystallographic B-factor workflow."""

    st.header("Experimental structures → B-factor")
    st.write(
        "For experimental crystallography, the B-factor is a temperature/disorder "
        "measure. Within comparable experiments, lower values generally indicate "
        "less atomic mobility; they are not pLDDT confidence scores."
    )

    with st.sidebar:
        refresh_from_rcsb = st.button("Refresh experimental data from RCSB")
        if refresh_from_rcsb:
            st.cache_data.clear()

    try:
        with st.spinner(
            "Refreshing from RCSB…"
            if refresh_from_rcsb
            else "Loading the bundled analysis…"
        ):
            results, errors, data_source = load_results(refresh_from_rcsb)
    except Exception as error:
        st.error("The experimental analysis could not load.")
        st.exception(error)
        st.stop()

    st.success(data_source)
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
            help="Lower is generally less mobile within comparable experimental data.",
        )

    st.subheader("B-factor comparison")
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
        results[display_columns].round(
            {"mean_b_factor": 2, "median_b_factor": 2, "stdev_b_factor": 2}
        ),
        width="stretch",
        hide_index=True,
    )
    st.download_button(
        "Download B-factor results CSV",
        data=results.to_csv(index=False),
        file_name="experimental_b_factor_metrics.csv",
        mime="text/csv",
    )

    with st.expander("Experimental dataset"):
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


def render_ai_mode() -> None:
    """Render the AI confidence workflow for uploaded or public PDB structures."""

    st.header("AI-predicted structures → pLDDT")
    st.write(
        "For RFdiffusion, ProteinMPNN/AlphaFold workflows, some pipelines store "
        "pLDDT in the PDB B-factor column. Higher pLDDT means higher predicted "
        "confidence."
    )
    st.warning(
        "Use this mode only when the structure producer confirms that the "
        "B-factor column contains pLDDT. The parser cannot infer provenance "
        "from a PDB file alone."
    )

    with st.form("ai_structure_form"):
        source_type = st.radio(
            "Structure source",
            ["RCSB PDB ID", "Upload a PDB file"],
            horizontal=True,
        )
        pdb_id = ""
        uploaded_file = None
        if source_type == "RCSB PDB ID":
            pdb_id = st.text_input("PDB ID", placeholder="Example: 7JZL")
        else:
            uploaded_file = st.file_uploader(
                "Upload an AI-predicted PDB file",
                type=["pdb", "ent", "txt"],
            )
        chain = st.text_input("Chain to analyze", value="A", max_chars=1)
        submitted = st.form_submit_button("Analyze pLDDT")

    if submitted:
        try:
            with st.spinner("Reading the AI structure…"):
                st.session_state["ai_result"] = analyze_ai_source(
                    pdb_id, chain, uploaded_file
                )
            st.session_state.pop("ai_error", None)
        except Exception as error:
            st.session_state.pop("ai_result", None)
            st.session_state["ai_error"] = str(error)

    if "ai_error" in st.session_state:
        st.error(st.session_state["ai_error"])
    result = st.session_state.get("ai_result")
    if not result:
        st.info(
            "Enter an RCSB PDB ID or upload a predicted PDB file, choose its "
            "chain, and select Analyze pLDDT."
        )
        return

    st.success(f"Analyzed {result['source']} chain {result['chain']}.")
    metric_columns = st.columns(4)
    metric_columns[0].metric("Mean pLDDT", f"{result['mean']:.2f}")
    metric_columns[1].metric("Median pLDDT", f"{result['median']:.2f}")
    metric_columns[2].metric(
        "High confidence (≥90)",
        f"{result['high_confidence_percent']:.1f}%",
    )
    metric_columns[3].metric(
        "Low confidence (<50)",
        f"{result['low_confidence_percent']:.1f}%",
    )

    figure = create_plddt_plot(result["values"])
    st.pyplot(figure)
    plt.close(figure)
    st.caption(
        "pLDDT interpretation: <50 very low, 50–70 low, 70–90 good, "
        "and ≥90 high confidence."
    )


def render_app() -> None:
    """Render the full Streamlit page."""

    st.set_page_config(
        page_title="Binder Structural Analysis",
        page_icon="🧬",
        layout="wide",
    )
    st.title("Computational Post-Analysis of Protein Structures")
    st.write(
        "Use the analysis mode that matches the provenance of your structure. "
        "The two metrics have opposite interpretations and are never combined "
        "into one score."
    )

    mode = st.sidebar.radio(
        "Analysis mode",
        [EXPERIMENTAL_MODE, AI_MODE],
        help="Choose B-factor for experimental structures or pLDDT for AI predictions.",
    )
    if mode == EXPERIMENTAL_MODE:
        render_experimental_mode()
    else:
        render_ai_mode()


if __name__ == "__main__":
    render_app()