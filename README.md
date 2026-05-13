# Mutation Analysis Pipeline

A Python pipeline for processing and analyzing mutation data from [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) outputs.

## Overview

This pipeline supports the following:

- **Data cleaning and filtering**: Removes low-coverage mutations, ancestor mutations, and nonstandard characters
- **Mutation classification**: Categorizes mutations as nonsynonymous, synonymous, intergenic, nonsense, noncoding, or pseudogene
- **Frequency-based filtering**: Generates filtered datasets at user-defined frequency thresholds
- **Statistical summaries**: Calculates mutation class proportions, average frequencies, and (optionally) average coverage per strain
- **Mutation analysis**: Identifies parallel mutations at the site and gene level, and unique mutations per strain
- **Neutral-model simulations**: Estimates expected dN/dS ratios and parallel-mutation counts under a null (neutral) model
- **Plotting**: Generates bubble plots, mutation spectrum charts, parallel-mutation heatmaps, and allele frequency distributions

## Requirements

- Python 3.9+
- pandas
- numpy
- matplotlib
- openpyxl
- gdtools*

\* *gdtools* is included in the *breseq* conda package.

## Installation

```bash
# Clone the repository
git clone https://github.com/Alebraco/mutation_analysis
cd mutation_analysis

# Create and activate the conda environment
conda create -n mutation_analysis -c bioconda pandas numpy matplotlib openpyxl breseq
conda activate mutation_analysis

# Optional: Install as pip package to use `mutanalysis` instead of `python main.py`
pip install .
```

After `pip install .`, the pipeline can be run as `mutanalysis` instead of `python main.py`.

## Usage

The pipeline has three subcommands:

| Subcommand | Purpose |
|------------|---------|
| `process` | Clean breseq outputs and run downstream analysis |
| `simulate-dnds` | Expected dN/dS under a neutral null model |
| `simulate-parallel` | Expected parallel mutation counts under a neutral null model |

---

### `process` — Clean and analyze mutation data

The `process` subcommand has two input modes depending on what data is available.

---

#### Mode 1 — Raw breseq sample directories

Use this mode if you have the original breseq output folders and a reference genome file.

**Expected directory structure:**

```
samples_dir/
├── sample_A/
│   └── data/
│       ├── output.gd
│       └── summary.json
├── sample_B/
│   └── data/
│       ├── output.gd
│       └── summary.json
└── ...
```

> **Note:** The reference file can be a GenBank (`.gbk`), FASTA (`.fasta`), or GFF (`.gff`) file.

```bash
mutanalysis process <ancestor> \
  --samples-dir <samples_dir> \
  --reference <reference.gbk> \
  [--gdtools <path/to/gdtools>] \
  [--output <output_dir>] \
  [--threshold 0.25 0.5 0.75 1.0] \
  [--plot bubble spectrum parallel allele]
```

**Example:**
```bash
mutanalysis process KZ_19 \
  --samples-dir samples \
  --reference reference.gbk \
  --output results \
  --threshold 0.25 0.5 1.0 \
  --plot spectrum parallel
```
> **Note:** The `process` keyword can be omitted for backwards compatibility:`mutanalysis KZ_19 --samples-dir ...` still works.

> **Note:** A single or multiple frequency threshold values may be used (e.g. `--threshold 0.75`, `--threshold 0.5 0.8`)

---

#### Mode 2 — Pre-processed mutation table

Use this mode if a mutation table is already available (generated from a previous run). Accepted formats: `.xlsx`, `.xls`, `.csv`, `.tsv`.

**Place your input file anywhere accessible and pass its path directly:**

```bash
mutanalysis process <ancestor> <input_file> \
  [--header-row <int>] \
  [--output <output_dir>] \
  [--threshold 0.25 0.5 0.75 1.0] \
  [--plot bubble spectrum parallel allele]
```

**Example:**
```bash
mutanalysis process KZ_19 sample_breseq_output.xlsx \
  --output results \
  --threshold 0.5 1.0 \
  --plot bubble
```

---

#### `process` arguments

| Argument | Mode | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `ancestor` | Both | Yes | — | Ancestor strain column name (e.g., `KZ_19`) |
| `--output` | Both | No | `output_files/` | Output directory for results |
| `--threshold` | Both | No | — | Space-separated frequency thresholds for filtering |
| `--samples-dir` | 1 | Yes* | — | Directory containing breseq sample folders |
| `--reference` | 1 | Yes* | — | Reference genome file used in the original breseq run |
| `--gdtools` | 1 | No | `gdtools` | Path to gdtools executable (not needed if using **conda**) |
| `input_file` | 2 | Yes* | — | Path to pre-processed mutation table |
| `--header-row` | 2 | No | `0` | Header row index in the input file |
| `--plot` | Both | No | — | Plot types: `bubble`, `spectrum`, `parallel`, `allele` |

\* Required within that mode.

---

### `simulate-dnds` — Expected dN/dS under neutrality

Runs a neutral-model simulation to estimate what fraction of mutations would be synonymous, nonsynonymous, nonsense, RNA genes, or intergenic purely by chance. Deviations from these expected values indicate selection.

Outputs **expected** values only (under neutral model), written to `<output>/expected/`. Compare against the observed mutation proportions in `mutation_summary.csv` from a prior `process` run.

**Requires `cleaned_data.csv` from a prior `process` run.**

```bash
mutanalysis simulate-dnds \
  --input <output_dir>/cleaned_data.csv \
  --reference <reference.gbk> \
  --ancestor <ancestor> \
  [--companion <companion.gff>] \
  [--output <output_dir>] \
  [--replicates 1000] \
  [--seed <int>]
```

**Example:**
```bash
mutanalysis simulate-dnds \
  --input results/cleaned_data.csv \
  --reference reference.gbk \
  --ancestor KZ_19 \
  --output results \
  --replicates 2000
```

#### `simulate-dnds` arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--input` | Yes | — | Path to `cleaned_data.csv` from `process` |
| `--reference` | Yes | — | Reference genome (`.gbk`, `.gbff`, `.fasta`, `.fna`, `.fa`, `.gff`, `.gff3`) |
| `--ancestor` | Yes | — | Ancestor column name (must match the `process` run) |
| `--companion` | No | — | Companion GFF/FASTA if `--reference` lacks sequence or annotations |
| `--output` | No | `output_files/` | Output directory |
| `--replicates` | No | `1000` | Number of replicates |
| `--seed` | No | — | Seed for reproducibility |

---

### `simulate-parallel` — Expected parallel mutations under neutrality

Estimates how many sites and genes would mutate independently in multiple strains at random (neutral expectation).

Outputs **expected** values only (under neutral model), written to `<output>/expected/`. Compare against the observed counts in `site_parallel_mutations.csv` / `gene_parallel_mutations.csv` from a prior `process` run.

**Requires `cleaned_data.csv` from a prior `process` run.**

```bash
mutanalysis simulate-parallel \
  --input <output_dir>/cleaned_data.csv \
  --reference <reference.gbk> \
  --ancestor <ancestor> \
  [--companion <companion.gff>] \
  [--output <output_dir>] \
  [--replicates 10000] \
  [--seed <int>]
```

**Example:**
```bash
mutanalysis simulate-parallel \
  --input results/cleaned_data.csv \
  --reference reference.gbk \
  --ancestor KZ_19 \
  --output results \
  --replicates 10000 \
  --seed 42
```

#### `simulate-parallel` arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--input` | Yes | — | Path to `cleaned_data.csv` from `process` |
| `--reference` | Yes | — | Reference genome (`.gbk`, `.gbff`, `.fasta`, `.fna`, `.fa`, `.gff`, `.gff3`) |
| `--ancestor` | Yes | — | Ancestor column name (must match the `process` run) |
| `--companion` | No | — | Companion GFF/FASTA if `--reference` lacks sequence or annotations |
| `--output` | No | `output_files/` | Output directory |
| `--replicates` | No | `10000` | Number of neutral-model replicates |
| `--seed` | No | — | Integer seed for reproducibility |

---

## Output Files

All outputs are written to the directory specified by `--output` (default: `output_files/`).

### 1. Cleaned Data

| File | Description |
|------|-------------|
| `cleaned_data.csv` | Filtered dataset with ancestor mutations and low-coverage rows removed |
| `low_coverage_rows.csv` | Rows flagged with low coverage (`?` values), if any exist |

### 2. Frequency-Filtered Tables

Generated for each threshold passed to `--threshold`:

| File | Description |
|------|-------------|
| `frequency_(threshold).csv` | Mutations that meet the frequency threshold in at least one strain |

### 3. Statistical Summary

| File | Description |
|------|-------------|
| `mutation_summary.csv` | Per-strain statistics: mutation type proportions, total count, average frequency, and average coverage* |

\* Average Coverage only available in Mode 1.

**Example output:**

| Line | Nonsynonymous | Synonymous | Intergenic | NonSense | Noncoding | Pseudogene | Total Mutations | Average Frequency | Average Coverage |
|------|---------------|------------|------------|----------|-----------|------------|-----------------|-------------------|-----------------|
| Strain_A | 0.65 | 0.20 | 0.10 | 0.03 | 0.02 | 0.00 | 120 | 0.875 | 142.3 |
| Strain_B | 0.58 | 0.25 | 0.12 | 0.04 | 0.01 | 0.00 | 98 | 0.820 | 138.7 |

### 4. Plots

Generated when `--plot` is used:

| File | Description |
|------|-------------|
| `bubble_plot.png` | Bubble plot of mutation frequency vs. nonsynonymous proportion (sized by total mutations) |
| `mutation_spectrum.png` | Stacked bar chart of mutation type proportions per group/timepoint |
| `parallel_mutation_heatmap.png` | Heatmap of genes mutated in parallel across strains |
| `allele_distribution.png` | Allele frequency distribution by group/timepoint |

### 5. Mutation Analysis

| File | Description |
|------|-------------|
| `site_parallel_mutations.csv` | Mutations shared at the exact same position across multiple strains |
| `gene_parallel_mutations.csv` | Mutations in the same gene across multiple strains (any position) |
| `unique_mutations.csv` | Mutations present in only one strain |

### 6. Simulation Outputs

Written to `<output_dir>/expected/` by `simulate-dnds` and `simulate-parallel`. These commands write only to the `expected/` subdirectory. The `--output` directory will not obtain new files. The matching **observed** values come from `process` (see §5 above).

| File | Subcommand | Description |
|------|------------|-------------|
| `expdNdS_<stem>.csv` | `simulate-dnds` | Expected mutation category percentages per strain |
| `avg_expdNdS_<stem>.csv` | `simulate-dnds` | Average expected percentages across all strains |
| `parallel_sites_<stem>.csv` | `simulate-parallel` | Expected number of parallel sites per strain pair |
| `parallel_genes_<stem>.csv` | `simulate-parallel` | Expected number of parallel genes per strain pair |

`<stem>` is derived from the input CSV filename (e.g. `cleaned_data` for `cleaned_data.csv`).

---

## Mutation Classification

Mutations are classified from the `annotation` column in the following priority order:

| Type | Criteria |
|------|----------|
| **Nonsense** | Annotation contains a stop codon (`*`) |
| **Noncoding** | Annotation contains "noncoding" |
| **Intergenic** | Annotation contains "intergenic" |
| **Pseudogene** | Annotation contains "pseudogene" |
| **Nonsynonymous (indel)** | Annotation contains "coding" |
| **Synonymous** | Contains `→` and codon change produces the same amino acid |
| **Nonsynonymous** | Contains `→` and codon change produces a different amino acid |
| **Unknown** | Annotation does not match any of the above |

---

## Project Structure

```
mutation_analysis/
├── main.py              # Orchestrating script
├── pyproject.toml       # Package config
└── modules/
    ├── analysis.py             # Parallel and unique mutation detection
    ├── data_cleaner.py         # Removes nonstandard characters
    ├── data_loader.py          # Loads, standardizes, and filters input data
    ├── gdtools_runner.py       # Launches gdtools COMPARE (Mode 1)
    ├── mutation_classifier.py  # Classifies mutations by type
    ├── plotting.py             # Bubble plot, spectrum, heatmap, and allele distribution charts
    ├── statistics.py           # Summary statistics and frequency filtering
    ├── utils.py                # Shared data structures
    └── simulation/
        ├── dnds.py             # dN/dS simulation
        ├── parallel.py         # Parallel mutation simulation
        ├── mutation_model.py   # Substitution probability estimation and replicate drawing
        ├── reference_loader.py # Parses GenBank, FASTA, and GFF reference files
        └── sequence_utils.py   # Helpers
```

---

## Troubleshooting

**Issue:** `AttributeError` or `KeyError` for ancestor column
- The ancestor name must match the column header exactly (case-insensitive, spaces treated as underscores)
- Use `--header-row` if your file has extra rows before the column headers

**Issue:** No mutations in output
- Check that the ancestor column name is correct. All rows matching the ancestor are excluded by design

**Issue:** `simulate-parallel` / `simulate-dnds` ran with no errors but the output directory looks empty
- Check the `expected/` subfolder inside `--output`. Both simulate commands only write there. Observed parallelism and dN/dS values are produced by `process`, not by the simulate commands.

**Issue:** Mode 1 fails with `gdtools` not found
- Ensure gdtools is installed and on your PATH, or pass the full path with `--gdtools /path/to/gdtools`
- The `--reference` argument must point to the genome file breseq was originally run against (`.gff`, `.fasta`, or `.gbk`), not a sample directory
