"""Compare structural metrics from curated public PDB structures.

This script intentionally does not run a structure-prediction model. It downloads
small public PDB files, extracts the B-factor column for a configured chain, and
summarizes the result. For predicted structures, B-factors may contain pLDDT;
for experimental structures, they are usually crystallographic temperature
factors and must not be interpreted as pLDDT.
"""

from __future__ import annotations

import argparse
import csv
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


OUTPUT_PLOT = Path("confidence_comparison.png")
OUTPUT_CSV = Path("structural_metrics.csv")


@dataclass(frozen=True)
class StructureTarget:
    """One structure and the chain that represents the analyzed binder."""

    pdb_id: str
    group: str
    chain: str
    description: str
    experimental_method: str
    metric_interpretation: str


# These entries were checked against the RCSB structure records. The analyzed
# chains are deliberately explicit because chain B is not a universal binder
# chain across PDB entries.
STRUCTURES = (
    StructureTarget(
        "6V67",
        "De novo miniprotein",
        "A",
        "De novo PD-1-binding miniprotein GR918.2 (apo)",
        "X-ray diffraction",
        "Experimental B-factor; not pLDDT",
    ),
    StructureTarget(
        "7S5B",
        "De novo miniprotein",
        "A",
        "De novo binder to human interleukin-7 receptor (unbound)",
        "X-ray diffraction",
        "Experimental B-factor; not pLDDT",
    ),
    StructureTarget(
        "8VEI",
        "De novo miniprotein",
        "A",
        "De novo designed colic-acid binder CHD_r1",
        "X-ray diffraction",
        "Experimental B-factor; not pLDDT",
    ),
    StructureTarget(
        "1MLC",
        "Natural antibody",
        "E",
        "Fab D44.1 bound to chicken egg-white lysozyme",
        "X-ray diffraction",
        "Experimental B-factor; not pLDDT",
    ),
    StructureTarget(
        "4JHW",
        "Natural antibody",
        "H",
        "RSV fusion glycoprotein bound to antibody D25",
        "X-ray diffraction",
        "Experimental B-factor; not pLDDT",
    ),
    StructureTarget(
        "1ADQ",
        "Natural antibody",
        "H",
        "Human rheumatoid-factor Fab bound to IgG Fc",
        "X-ray diffraction",
        "Experimental B-factor; not pLDDT",
    ),
)


def download_pdb(pdb_id: str, destination: Path) -> None:
    """Download one PDB file from RCSB with a bounded network timeout."""

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "binder-structural-analysis/1.0"},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        destination.write_bytes(response.read())


def extract_b_factors(pdb_path: Path, target_chain: str) -> list[float]:
    """Extract atom B-factors for one chain from fixed-width ATOM records."""

    b_factors: list[float] = []
    with pdb_path.open("r", encoding="ascii", errors="replace") as pdb_file:
        for line in pdb_file:
            if not line.startswith("ATOM") or len(line) < 66:
                continue
            if line[21].strip() != target_chain:
                continue
            try:
                b_factors.append(float(line[60:66].strip()))
            except ValueError:
                # Malformed individual records should not discard the whole
                # structure, but the skipped value is visible in atom_count.
                continue
    return b_factors


def summarize_target(target: StructureTarget, work_dir: Path) -> dict[str, object]:
    """Download and summarize a configured PDB chain."""

    pdb_path = work_dir / f"{target.pdb_id}.pdb"
    download_pdb(target.pdb_id, pdb_path)
    values = extract_b_factors(pdb_path, target.chain)
    if not values:
        raise ValueError(
            f"No ATOM B-factors found for {target.pdb_id} chain {target.chain}"
        )

    return {
        "pdb_id": target.pdb_id,
        "group": target.group,
        "chain": target.chain,
        "description": target.description,
        "experimental_method": target.experimental_method,
        "metric_interpretation": target.metric_interpretation,
        "atom_count": len(values),
        "mean_b_factor": mean(values),
        "median_b_factor": median(values),
        "stdev_b_factor": pstdev(values),
        "minimum_b_factor": min(values),
        "maximum_b_factor": max(values),
    }


def write_csv(rows: Iterable[dict[str, object]], output_path: Path) -> None:
    """Write machine-readable results for downstream analysis."""

    rows = list(rows)
    if not rows:
        return
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def create_plot(rows: list[dict[str, object]], output_path: Path) -> None:
    """Create a compact comparison plot from the extracted metrics."""

    data = pd.DataFrame(rows)
    sns.set_theme(style="whitegrid")
    figure, axis = plt.subplots(figsize=(9, 6))
    sns.boxplot(
        data=data,
        x="group",
        y="mean_b_factor",
        hue="group",
        palette="Set2",
        legend=False,
        ax=axis,
    )
    sns.stripplot(
        data=data,
        x="group",
        y="mean_b_factor",
        color="black",
        size=7,
        jitter=0.12,
        ax=axis,
    )
    axis.set_title("Mean PDB B-factor by structure group", fontsize=14, weight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("Mean atom B-factor (source-specific units)")
    axis.tick_params(axis="x", rotation=10)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Write CSV results but skip PNG generation.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=OUTPUT_CSV,
        help=f"CSV output path (default: {OUTPUT_CSV}).",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=OUTPUT_PLOT,
        help=f"PNG output path (default: {OUTPUT_PLOT}).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows: list[dict[str, object]] = []

    print("Downloading and analyzing curated PDB chains...")
    with tempfile.TemporaryDirectory(prefix="binder-structural-analysis-") as temp:
        work_dir = Path(temp)
        for target in STRUCTURES:
            try:
                row = summarize_target(target, work_dir)
                rows.append(row)
                print(
                    f"  {target.pdb_id} chain {target.chain}: "
                    f"mean B-factor {row['mean_b_factor']:.2f} "
                    f"({row['atom_count']} atoms)"
                )
            except (OSError, ValueError, urllib.error.URLError) as error:
                print(f"  ERROR {target.pdb_id} chain {target.chain}: {error}")

    if not rows:
        raise RuntimeError("No structures could be analyzed.")

    write_csv(rows, args.csv)
    if not args.no_plot:
        create_plot(rows, args.plot)

    print(f"\nSaved tabular results to {args.csv}")
    if not args.no_plot:
        print(f"Saved plot to {args.plot}")

    summary = (
        pd.DataFrame(rows)
        .groupby("group", sort=False)["mean_b_factor"]
        .agg(["count", "mean", "median"])
    )
    print("\nGroup summary:")
    print(summary.to_string(float_format=lambda value: f"{value:.2f}"))
    print(
        "\nInterpretation note: these entries are experimental X-ray structures. "
        "Their B-factors are crystallographic temperature factors, not pLDDT. "
        "Do not claim that one group has higher AI confidence from this plot."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())