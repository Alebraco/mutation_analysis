# Mutation Analysis Pipeline

A Python pipeline for processing and analyzing mutation data from [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) outputs.

## Overview

This pipeline processes breseq outputs in Excel format and performs statistical and comparative analyses:

- **Data cleaning and filtering**: Removes low-coverage mutations, ancestor mutations, and deletions
- **Mutation classification**: Categorizes mutations as nonsynonymous, synonymous, intergenic, nonsense, noncoding, or pseudogene
- **Frequency-based filtering**: Generates filtered datasets at multiple frequency thresholds (25%, 50%, 75%, 100%)
- **Statistical summaries**: Calculates mutation type proportions and average frequencies per strain
- **Parallel mutation detection**: Identifies mutations shared across strains at both site and gene levels

## Requirements

- Python 3
- pandas
- numpy

## Installation

```bash
# Clone the repository
git clone https://github.com/Alebraco/mutation_analysis
cd mutation_analysis

# Install dependencies
pip install pandas numpy
```

## Usage

### Basic Usage

```bash
python main.py <input_file> <ancestor_name>
```

**Example:**
```bash
python main.py breseq_output.xlsx KZ_19
```

### Command-Line Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `input_file` | Yes | - | Path to input Excel file containing breseq mutation data |
| `ancestor` | Yes | - | Name of the ancestor strain column (e.g., "KZ_19") |
| `--header-row` | No | 1 | Row index containing column headers |
| `--output-dir` | No | `.` | Output directory for results |

## Input Format

The input Excel file should contain a standard breseq output.

## Output Files

The pipeline generates the following outputs in the specified output directory:

### 1. Data Files

| File | Description |
|------|-------------|
| `cleaned_data.xlsx` | Filtered dataset with ancestor mutations, low coverage, and deletions removed |
| `low_coverage_rows.xlsx` | Rows flagged with low coverage (contains `?` in data) |

### 2. Frequency-Filtered Tables

| File | Description |
|------|-------------|
| `frequency_0.25.xlsx` | Mutations present at ≥25% frequency in at least one strain |
| `frequency_0.5.xlsx` | Mutations present at ≥50% frequency in at least one strain |
| `frequency_0.75.xlsx` | Mutations present at ≥75% frequency in at least one strain |
| `frequency_1.0.xlsx` | Mutations fixed (100% frequency) in at least one strain |

### 3. Statistical Summaries

| File | Description |
|------|-------------|
| `mutation_summary.xlsx` | Per-strain statistics including mutation type proportions, total count, and average frequency |

**Example output:**

| Line | Nonsynonymous | Synonymous | Intergenic | NonSense | Noncoding | Total Mutations | Average Frequency |
|------|---------------|------------|------------|----------|-----------|-----------------|-------------------|
| Strain_A | 0.65 | 0.20 | 0.10 | 0.03 | 0.02 | 120 | 0.875 |
| Strain_B | 0.58 | 0.25 | 0.12 | 0.04 | 0.01 | 98 | 0.820 |

### 4. Parallel Mutation Analysis

| File | Description |
|------|-------------|
| `site_parallel_mutations.xlsx` | Shared mutations occurring at the same position |
| `gene_parallel_mutations.xlsx` | Shared mutations occurring at same gene (any site) |

## Mutation Classification

Mutations are classified based on their annotation:

| Type | Criteria |
|------|----------|
| **Nonsense** | Contains stop codon (`*`) |
| **Noncoding** | Annotation contains "noncoding" |
| **Intergenic** | Annotation contains "intergenic" |
| **Pseudogene** | Annotation contains "pseudogene" |
| **Synonymous** | Codon change with same amino acid |
| **Nonsynonymous** | Codon change with different amino acid |

## Module Overview

### Core Modules

- **`main.py`**: Orchestrating script
- **`data_loader.py`**: Loads and filters raw mutation data
- **`data_cleaner.py`**: Cleans non-standard characters from data
- **`mutation_classifier.py`**: Classifies mutations
- **`statistics.py`**: Calculates summary statistics and frequency filtering
- **`analysis.py`**: Performs mutation analysis (shared and unique)
- **`utils.py`**: Utility functions and common data structures

## Example Workflow

```bash
# Run the full analysis pipeline
python main.py breseq_example_output.xlsx KZ_19 --output-dir results/
```

## Troubleshooting

**Issue**: `KeyError` for ancestor column
- Ancestor name must match exactly (case-sensitive) with the column name in the input Excel file

**Issue**: `ValueError` when parsing frequencies
- Frequency values must be numeric (0.0 to 1.0) or empty