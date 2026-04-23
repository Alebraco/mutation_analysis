#!/usr/bin/env python3
"""
Mutation Analysis Pipeline.

Subcommands:
    process           Process breseq outputs (or an existing mutation table) and run downstream analysis.
    simulate-dnds     Expected dN/dS under a neutral null model.
    simulate-parallel Expected parallel-site and parallel-gene counts under a neutral null model.
"""

import argparse
import os
import re
import sys

SUBCOMMANDS = ('process', 'simulate-dnds', 'simulate-parallel')


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='mutanalysis',
        description='Analyze mutation data from breseq output files and run neutral-model simulations.',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    _add_process_parser(subparsers)
    _add_simulate_dnds_parser(subparsers)
    _add_simulate_parallel_parser(subparsers)

    return parser


def _add_process_parser(subparsers) -> None:
    p = subparsers.add_parser(
        'process',
        help='Process breseq outputs (or a pre-built mutation table) and run analysis.',
    )
    p.add_argument('ancestor', help='Ancestor name (e.g., "KZ_19")')
    p.add_argument('--output', default='output_files',
                   help='Output directory (default: output_files/)')
    p.add_argument('--threshold', nargs='+', type=float, default=[],
                   help='Mutation filtering based on frequency threshold (default: none)')
    p.add_argument('--plot', nargs='+', choices=['bubble', 'spectrum'],
                   help='Generate plots: bubble (summary bubble plot), spectrum (stacked mutation types)')

    mode1 = p.add_argument_group('Mode 1: raw breseq output')
    mode1.add_argument('--samples-dir', metavar='DIR',
                       help='Directory containing breseq sample folders (data/output.gd)')
    mode1.add_argument('--reference',
                       help='Reference genome file for gdtools COMPARE')
    mode1.add_argument('--gdtools', default='gdtools',
                       help='Path to gdtools executable (default: gdtools)')

    mode2 = p.add_argument_group('Mode 2: pre-processed mutation table')
    mode2.add_argument('input_file', nargs='?',
                       help='Input file path (.xlsx, .tsv, .csv)')
    mode2.add_argument('--header-row', type=int, default=0,
                       help='Header row index for input file (default: 0)')


def _add_simulate_dnds_parser(subparsers) -> None:
    p = subparsers.add_parser(
        'simulate-dnds',
        help='Expected dN/dS under neutrality (null model).',
    )
    _add_simulation_common_args(p, default_replicates=1000)


def _add_simulate_parallel_parser(subparsers) -> None:
    p = subparsers.add_parser(
        'simulate-parallel',
        help='Expected parallel-mutation counts under neutrality (null model).',
    )
    _add_simulation_common_args(p, default_replicates=10000)


def _add_simulation_common_args(p: argparse.ArgumentParser, default_replicates: int) -> None:
    p.add_argument('--input', required=True,
                   help='Path to cleaned_data.csv produced by `mutanalysis process`.')
    p.add_argument('--reference', required=True,
                   help='Reference genome: .gbk/.gb/.gbff, .fasta/.fna/.fa, or .gff/.gff3.')
    p.add_argument('--gff', default=None,
                   help='Companion GFF/FASTA file if the --reference does not include both sequence and annotations.')
    p.add_argument('--ancestor', required=True,
                   help='Ancestor column name (matches the argument given to `mutanalysis process`).')
    p.add_argument('--output', default='output_files',
                   help='Output directory (default: output_files/)')
    p.add_argument('--replicates', type=int, default=default_replicates,
                   help=f'Number of replicates (default: {default_replicates}).')
    p.add_argument('--seed', type=int, default=None,
                   help='Random seed for reproducibility (default: nondeterministic).')


def _run_process(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    from modules.data_loader import load_and_filter
    from modules.statistics import calculate_basic_stats, frequency_filter
    from modules.analysis import mutation_analysis

    use_mode1 = bool(args.samples_dir or args.reference)
    use_mode2 = bool(args.input_file)

    if use_mode1 and use_mode2:
        parser.error('Cannot use --samples-dir/--reference together with input_file.')
    if not use_mode1 and not use_mode2:
        parser.error('Provide either input_file (processed file) or both --samples-dir and --reference (raw breseq outputs).')
    if use_mode1 and not (args.samples_dir and args.reference):
        parser.error('Both --samples-dir and --reference are required.')

    os.makedirs(args.output, exist_ok=True)

    if use_mode1:
        from modules.gdtools_runner import find_gd_files, find_summary_jsons, run_gdtools_compare
        gd_files = find_gd_files(args.samples_dir)
        if not gd_files:
            sys.exit(f'No output.gd files found in {args.samples_dir}')
        json_files = find_summary_jsons(args.samples_dir)
        compare_path = os.path.join(args.output, 'mutation_data.csv')
        run_gdtools_compare(args.gdtools, args.reference, gd_files, compare_path)
        input_file = compare_path
        header_row = 0
    else:
        input_file = args.input_file
        header_row = args.header_row
        json_files = None

    df_clean, question_df, ancestor = load_and_filter(input_file, args.ancestor, header_row)

    clean_file = os.path.join(args.output, 'cleaned_data.csv')
    df_clean.to_csv(clean_file, index=False)
    print(f'Cleaned data contains {len(df_clean)} rows.')
    print(f'  Saved: {clean_file}')

    if question_df is not None:
        question_file = os.path.join(args.output, 'low_coverage_rows.csv')
        question_df.to_csv(question_file, index=False)
        print(f'  Saved: {question_file}')

    calculate_basic_stats(df_clean, ancestor, args.output, json_files=json_files)

    if args.plot:
        summary_file = os.path.join(args.output, 'mutation_summary.csv')
        if 'bubble' in args.plot:
            from modules.plotting import plot_mutation_bubble
            plot_mutation_bubble(summary_file,
                                 output_file=os.path.join(args.output, 'bubble_plot.png'),
                                 show=False)
        if 'spectrum' in args.plot:
            from modules.plotting import plot_mutation_spectrum
            plot_mutation_spectrum(summary_file,
                                   output_file=os.path.join(args.output, 'mutation_spectrum.png'),
                                   show=False)
        print("Plots saved to output directory.")

    for threshold in args.threshold:
        df_filtered = frequency_filter(df_clean, threshold, ancestor)
        filtered_file = os.path.join(args.output, f'frequency_{threshold}.csv')
        df_filtered.to_csv(filtered_file, index=False)

    mutation_analysis(df_clean, ancestor, args.output)
    print("Analysis complete.")


def _run_simulation(args: argparse.Namespace, parser: argparse.ArgumentParser, kind: str) -> None:
    import pandas as pd

    if args.replicates < 1:
        parser.error('--replicates must be >= 1.')
    if args.replicates > 50_000:
        print(f'[warn] --replicates={args.replicates}: runtime may be very long.', file=sys.stderr)

    os.makedirs(args.output, exist_ok=True)
    df_clean = pd.read_csv(args.input)

    ancestor = re.sub(r'\s+', '_', str(args.ancestor).lower().strip())
    file_stem = os.path.splitext(os.path.basename(args.input))[0]

    if kind == 'simulate-dnds':
        from modules.simulation import run_dnds_simulation
        run_dnds_simulation(
            df_clean=df_clean,
            ancestor=ancestor,
            reference_path=args.reference,
            output_dir=args.output,
            file_stem=file_stem,
            gff_path=args.gff,
            n_replicates=args.replicates,
            seed=args.seed,
        )
    elif kind == 'simulate-parallel':
        from modules.simulation import run_parallel_simulation
        run_parallel_simulation(
            df_clean=df_clean,
            ancestor=ancestor,
            reference_path=args.reference,
            output_dir=args.output,
            file_stem=file_stem,
            gff_path=args.gff,
            n_replicates=args.replicates,
            seed=args.seed,
        )
    else:
        parser.error(f'Unknown simulation kind: {kind}')


_TOP_LEVEL_PASSTHROUGH = {'-h', '--help'}


def main():
    argv = sys.argv[1:]
    # Back-compat shim: if the first token isn't a known subcommand and isn't a top-level
    # help/version flag, assume the caller is using the old flat CLI and inject `process`.
    if argv and argv[0] not in SUBCOMMANDS and argv[0] not in _TOP_LEVEL_PASSTHROUGH:
        argv = ['process'] + argv

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == 'process':
        _run_process(args, parser)
    elif args.command in ('simulate-dnds', 'simulate-parallel'):
        _run_simulation(args, parser, args.command)


if __name__ == '__main__':
    main()
