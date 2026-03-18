# Mutation Analysis Pipeline

A Python pipeline for processing and analyzing mutation data from [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) outputs.

## Overview

This pipeline supports two input modes and produces a set of mutation tables:

- **Data cleaning and filtering**: Removes low-coverage mutations, ancestor mutations, and nonstandard characters
- **Mutation classification**: Categorizes mutations as nonsynonymous, synonymous, intergenic, nonsense, noncoding, or pseudogene
- **Frequency-based filtering**: Generates filtered datasets at user-defined frequency thresholds
- **Statistical summaries**: Calculates mutation class proportions, average frequencies, and (optionally) average coverage per strain
- **Mutation analysis**: Identifies parallel mutations at the site and gene level, and unique mutations per strain

## Requirements

- Python 3
- pandas
- numpy
- openpyxl
- gdtools*

\* gdtools is included in the breseq conda package.

## Installation

```bash
# Clone the repository
git clone https://github.com/Alebraco/mutation_analysis
cd mutation_analysis

# Install dependencies
conda create -n mutation_analysis -c bioconda pandas numpy openpyxl breseq
conda activate mutation_analysis
```

## Usage

The pipeline has two input modes depending on what data is available.

---

### Mode 1 — Raw breseq sample directories

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
python main.py <ancestor> \
  --samples-dir <samples_dir> \
  --reference <reference.gbk> \
  [--gdtools <path/to/gdtools>] \
  [--output <output_dir>] \
  [--threshold 0.25 0.5 0.75 1.0]
```

**Example:**
```bash
python main.py KZ_19 \
  --samples-dir samples \
  --reference reference.gbk \
  --output results \
  --threshold 0.25 0.5 1.0
```
> **Note:** A single or multiple frequency threshold values may be used (e.g. `--threshold 0.75`; `--threshold 0.5 0.8`)

---

### Mode 2 — Pre-processed mutation table

Use this mode if a mutation table is already available (generated from a previous run). Accepted formats: `.xlsx`, `.xls`, `.csv`, `.tsv`.

**Place your input file anywhere accessible and pass its path directly:**

```bash
python main.py <ancestor> <input_file> \
  [--header-row <int>] \
  [--output <output_dir>] \
  [--threshold 0.25 0.5 0.75 1.0]
```

**Example:**
```bash
python main.py KZ_19 sample_breseq_output.xlsx \
  --output results \
  --threshold 0.5 1.0
```

---

### All Arguments

| Argument | Mode | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `ancestor` | Both | Yes | — | Ancestor strain column name (e.g., `KZ_19`) |
| `--output` | Both | No | `output_files/` | Output directory for results |
| `--threshold` | Both | No | — | Space-separated frequency thresholds for filtering |
| `--samples-dir` | 1 | Yes* | — | Directory containing breseq sample folders |
| `--reference` | 1 | Yes* | — | Reference genome file used in the original breseq run |
| `--gdtools` | 1 | No | `gdtools` | Path to gdtools executable (not needed if using **conda**)|
| `input-file` | 2 | Yes* | — | Path to pre-processed mutation table |
| `--header-row` | 2 | No | `0` | Header row index in the input file |

\* Required within that mode.

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

### 4. Mutation Analysis

| File | Description |
|------|-------------|
| `site_parallel_mutations.csv` | Mutations shared at the exact same genomic position across multiple strains |
| `gene_parallel_mutations.csv` | Mutations in the same gene across multiple strains (any position) |
| `unique_mutations.csv` | Mutations present in exactly one strain |

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
└── modules/
    ├── analysis.py          # Parallel and unique mutation detection
    ├── data_cleaner.py      # Removes nonstandard characters
    ├── data_loader.py       # Loads, standardizes, and filters input data
    ├── gdtools_runner.py    # Launches gdtools COMPARE (Mode 1)
    ├── mutation_classifier.py  # Classifies mutations by type
    ├── statistics.py        # Summary statistics and frequency filtering
    └── utils.py             # Shared data structures (e.g. codon table)
```

---

## Troubleshooting

**Issue:** `AttributeError` or `KeyError` for ancestor column
- The ancestor name must match the column header exactly (case-insensitive, spaces treated as underscores)
- Use `--header-row` if your file has extra rows before the column headers

**Issue:** No mutations in output
- Check that the ancestor column name is correct. All rows matching the ancestor are excluded by design

**Issue:** Mode 1 fails with `gdtools` not found
- Ensure gdtools is installed and on your PATH, or pass the full path with `--gdtools /path/to/gdtools`
- The `--reference` argument must point to the genome file breseq was originally run against (`.gff`, `.fasta`, or `.gbk`), not a sample directory
