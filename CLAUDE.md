# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

- Install (editable): `pip install -e .` — exposes the `mutanalysis` entrypoint (see [pyproject.toml](pyproject.toml): `mutanalysis = "main:main"`).
- Run without install: `python main.py <subcommand> ...`.
- Conda environment (required for Mode 1, which shells out to `gdtools` from breseq):
  `conda create -n mutation_analysis -c bioconda pandas numpy matplotlib openpyxl breseq`.
- No test suite, linter, or formatter is configured in this repo.

## CLI shape

`main.py` defines three subcommands ([main.py:16](main.py:16)):

- `process` — clean breseq output and run downstream analysis.
- `simulate-dnds` — neutral-model expected dN/dS.
- `simulate-parallel` — neutral-model expected parallel-site/gene counts.

Non-obvious behaviors:

- **Back-compat shim** at [main.py:212](main.py:212): if the first argv token isn't a known subcommand or a top-level help flag, `process` is injected. The legacy flat CLI shown in the README (`mutanalysis KZ_19 --samples-dir ...`) still works.
- **`process` mode validation** in `_run_process` ([main.py:99](main.py:99)): Mode 1 (`--samples-dir` + `--reference`) and Mode 2 (positional `input_file`) are mutually exclusive; at least one is required.
- **Simulations consume `cleaned_data.csv`** produced by `process`, and require `--reference` + `--ancestor` matching the `process` run.

## Architecture

`main.py` is a thin orchestrator. All real work lives in `modules/`.

**Mode 1 pipeline** (raw breseq outputs):
`gdtools_runner.run_gdtools_compare` (invokes `gdtools COMPARE`) → writes `mutation_data.csv` → `data_loader.load_and_filter` cleans/standardizes and drops ancestor + low-coverage rows → `statistics.calculate_basic_stats` writes `mutation_summary.csv` (uses each sample's `summary.json` to add average coverage, Mode 1 only) → `statistics.frequency_filter` (one CSV per `--threshold`) → `analysis.mutation_analysis` writes parallel-site, parallel-gene, and unique-mutation tables.

**Mode 2** skips `gdtools_runner` and feeds the user-supplied table directly into `data_loader.load_and_filter`.

**Classification**: `mutation_classifier` applies the priority-ordered rules from the README (nonsense → noncoding → intergenic → pseudogene → nonsynonymous/synonymous via codon table in `utils.py`).

**Plotting** is imported lazily ([main.py:143](main.py:143)) only when `--plot` is set, so matplotlib isn't required otherwise.

**`modules/simulation/`** (not covered in README): neutral-model simulations. `run_dnds_simulation` and `run_parallel_simulation` are re-exported from `modules/simulation/__init__.py`; internals are split across `dnds.py`, `parallel.py`, `mutation_model.py`, `reference_loader.py`, `sequence_utils.py`.

## Conventions

- **Ancestor name normalization**: both `process` and `simulate-*` normalize the ancestor with `re.sub(r'\s+', '_', s.lower().strip())` (see [main.py:174](main.py:174) and `data_loader`). Pass the same string across subcommands so the column lookup succeeds.
- **Reference files**: accepted extensions are `.gbk/.gb/.gbff`, `.fasta/.fna/.fa`, `.gff/.gff3`. Use `--gff` on `simulate-*` when sequence and annotations live in separate files.
- **Output**: defaults to `output_files/`, created if missing. Simulation outputs are prefixed with the input CSV's basename (`file_stem`).

## Pointers

See [README.md](README.md) for the full argument table, output-file catalog, and mutation classification rules — don't duplicate them here.
