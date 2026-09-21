# Computational Post-Analysis of De Novo Miniprotein Structures

This project downloads a small, curated set of public Protein Data Bank (PDB)
structures and compares the mean B-factor of a specific chain. It demonstrates
PDB retrieval, fixed-width structural parsing, summary statistics, and
visualization without running GPU-heavy models such as RFdiffusion or
ProteinMPNN.

## Scientific scope

The original project idea described the PDB B-factor as a universal pLDDT
confidence field. That is only true for some predicted-structure files where
the producer intentionally stores pLDDT in that column. The structures in this
curated dataset are experimental X-ray structures, so their B-factors are
**crystallographic temperature factors**. They reflect atomic disorder and
refinement behavior, not AlphaFold confidence.

The analysis therefore labels the result as **mean B-factor**, not pLDDT, and
prints a warning against treating the two groups as directly comparable AI
confidence scores. The same parser can be reused for predicted PDB files after
their provenance and B-factor convention have been verified.

## Dataset

| Group | PDB | Chain | Structure |
| --- | --- | --- | --- |
| De novo miniprotein | `6V67` | `A` | De novo PD-1-binding miniprotein GR918.2 |
| De novo miniprotein | `7S5B` | `A` | De novo binder to human interleukin-7 receptor |
| De novo miniprotein | `8VEI` | `A` | De novo designed colic-acid binder CHD_r1 |
| Natural antibody | `1MLC` | `E` | Fab D44.1 bound to lysozyme |
| Natural antibody | `4JHW` | `H` | RSV fusion glycoprotein bound to antibody D25 |
| Natural antibody | `1ADQ` | `H` | Rheumatoid-factor Fab bound to IgG Fc |

PDB IDs and chain assignments are explicit in `main.py`. Chain `B` is not a
universal binder chain: chain identifiers vary between structures.

## Run

```bash
python -m pip install -r requirements.txt
python main.py
```

## Streamlit deployment

The repository includes `app.py` as the Streamlit entrypoint. To deploy on
Streamlit Community Cloud:

1. Select the repository
   `amrutarawool95-spec/binder-structural-analysis`.
2. Set the branch to `main`.
3. Set **Main file path** to `app.py`.
4. Deploy.

The app loads the checked-in `structural_metrics.csv` snapshot first, so the
page renders immediately even when Streamlit Cloud has slow or restricted
outbound network access. It has two separate analysis modes:

- **Experimental structures → B-factor**: lower generally means less atomic
  mobility within comparable crystallographic experiments.
- **AI-predicted structures → pLDDT**: higher means greater predicted
  confidence, when the prediction pipeline stored pLDDT in the B-factor field.

The **Refresh experimental data from RCSB** button is optional; it downloads the
experimental PDB files, caches live results for one hour, and falls back to the
bundled snapshot if a refresh fails. The AI mode accepts either an uploaded PDB
file or an RCSB PDB ID and chain. All Python dependencies, including
`streamlit`, are listed in `requirements.txt`.

The repository also keeps `main.py` compatible with Streamlit Cloud if an
existing deployment still points to that file: it automatically hands off to
the same UI when run by Streamlit, while remaining a command-line script when
run with Python.

To run it locally:

```bash
streamlit run app.py
```

The script creates:

- `structural_metrics.csv` — one row per analyzed PDB chain
- `confidence_comparison.png` — a boxplot with individual structure points

To generate only the CSV:

```bash
python main.py --no-plot
```

## Parsing workflow

1. Download each PDB file from the RCSB file service.
2. Read `ATOM` records using the PDB fixed-width format.
3. Filter to the configured chain.
4. Extract columns 61–66, the PDB B-factor field.
5. Calculate mean, median, standard deviation, minimum, and maximum.
6. Save a CSV and visualization.

The downloader uses a temporary directory and removes the files automatically
when the run finishes.

## Responsible extension

For a valid pLDDT comparison, use predicted PDB files whose documentation
confirms that pLDDT was written to the B-factor column. Keep experimental
temperature factors and predicted pLDDT in separate columns or analyses. A
strong follow-up would compare experimentally solved structures against their
corresponding prediction confidence, rather than comparing pLDDT to antibody
crystallographic B-factors.

## Sources

- [RCSB PDB](https://www.rcsb.org/)
- [RCSB PDB file download service](https://files.rcsb.org/)
- [PDB format specification](https://www.wwpdb.org/documentation/file-format)