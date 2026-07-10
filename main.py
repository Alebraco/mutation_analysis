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

def build_parser():
    parser = argparse.ArgumentParser(
        prog='mutanalysis',
        description='Analyze mutation data from breseq output files and run neutral-model simulations.',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    add_process_parser(subparsers)
    add_simulate_dnds_parser(subparsers)
    add_simulate_parallel_parser(subparsers)

    return parser


def add_process_parser(subparsers):
    p = subparsers.add_parser(
        'process',
        help='Process breseq outputs (or a pre-built mutation table) and run analysis.',
    )
    p.add_argument('ancestor', help='Ancestor name (e.g., "KZ_19")')
    p.add_argument('--output', default='output_files',
                   help='Output directory (default: output_files/)')
    p.add_argument('--reference',
                   help='Reference genome (.gbk/.gb/.gbff, .fasta/.fna/.fa, or .gff/.gff3). '
                        'Required for Mode 1, optional in Mode 2 (enables --pseudoclones).')
    p.add_argument('--threshold', nargs='+', type=float, default=[],
                   help='Mutation filtering based on frequency threshold (default: none)')
    p.add_argument('--pseudoclones', action='store_true',
                   help='Generate one pseudoclone genome per sample and frequency bin '
                        '(clonal deconvolution / contamination check). Requires --reference.')
    p.add_argument('--pseudoclone-bin-width', type=float, default=0.1,
                   help='Bin width (frequency) for pseudoclones (default: 0.1).')
    p.add_argument('--pseudoclone-min-freq', type=float, default=0.05,
                   help='Minimum frequency for a SNP to enter a pseudoclone '
                        '(default: 0.05).')
    p.add_argument('--extract-sequences', action='store_true',
                   help='Write nucleotide/protein FASTA files for mutated genes. '
                        'Requires --reference.')
    p.add_argument('--plot', nargs='+',
                   choices=['bubble', 'spectrum', 'parallel', 'allele', 'trajectory', 'genome'],
                   help='Generate plots: bubble (summary bubble plot), spectrum (stacked mutation types), '
                        'parallel (gene-level parallel-mutation heatmap), '
                        'allele (allele frequency distribution by group/day), '
                        'trajectory (total mutations over time per group), '
                        'genome (Kosterlitz zoomed genome plot of shared mutations; requires --reference)')
    p.add_argument('--genome-min-strains', type=float, default=0.25,
                   help='For --plot genome: minimum strains a genomic region must have mutations in'
                        'to be shown and focus on parallel mutations. A fraction in (0,1) is '
                        'a proportion of the strains. An integer is an absolute strain count.')
    p.add_argument('--specificity-permutations', type=int, default=10000,
                   help='Permutations for the treatment-specificity Dice test '
                        '(default: 10000; 0 skips the permutation p-value)')
    p.add_argument('--seed', type=int, default=None,
                   help='Random seed for the specificity permutation test (default: none)')

    mode1 = p.add_argument_group('Mode 1: raw breseq output')
    mode1.add_argument('--samples-dir', metavar='DIR',
                       help='Directory containing breseq sample folders (data/output.gd)')
    mode1.add_argument('--gdtools', default='gdtools',
                       help='Path to gdtools executable (default: gdtools)')

    mode2 = p.add_argument_group('Mode 2: pre-processed mutation table')
    mode2.add_argument('input_file', nargs='?',
                       help='Input file path (.xlsx, .tsv, .csv)')
    mode2.add_argument('--header-row', type=int, default=0,
                       help='Header row index for input file (default: 0)')


def add_simulate_dnds_parser(subparsers):
    p = subparsers.add_parser(
        'simulate-dnds',
        help='Expected dN/dS under neutrality (null model).',
    )
    add_simulation_common_args(p, default_replicates=1000)


def add_simulate_parallel_parser(subparsers):
    p = subparsers.add_parser(
        'simulate-parallel',
        help='Expected parallel-mutation counts under neutrality (null model).',
    )
    add_simulation_common_args(p, default_replicates=10000)


def add_simulation_common_args(parser, default_replicates):
    parser.add_argument('--input', required=True,
                   help='Path to cleaned_data.csv produced by `mutanalysis process`.')
    parser.add_argument('--reference', required=True,
                   help='Reference genome: .gbk/.gb/.gbff, .fasta/.fna/.fa, or .gff/.gff3.')
    parser.add_argument('--companion', default=None,
                   help='Companion GFF/FASTA file if the --reference does not include both sequence and annotations.')
    parser.add_argument('--ancestor', required=True,
                   help='Ancestor column name (matches the argument given to `mutanalysis process`).')
    parser.add_argument('--output', default='output_files',
                   help='Output directory (default: output_files/)')
    parser.add_argument('--replicates', type=int, default=default_replicates,
                   help=f'Number of replicates (default: {default_replicates}).')
    parser.add_argument('--seed', type=int, default=None,
                   help='Random seed for reproducibility (default: none).')


def run_process(args, parser):
    from modules.data_loader import load_and_filter
    from modules.statistics import calculate_basic_stats, frequency_filter
    from modules.analysis import mutation_analysis

    use_mode1 = bool(args.samples_dir)
    use_mode2 = bool(args.input_file)

    if use_mode1 and use_mode2:
        parser.error('Cannot use --samples-dir together with input_file.')
    if not use_mode1 and not use_mode2:
        parser.error('Provide either input_file (processed file) or --samples-dir with --reference (raw breseq outputs).')
    if use_mode1 and not args.reference:
        parser.error('--reference is required together with --samples-dir.')
    if args.pseudoclones and not args.reference:
        parser.error('--pseudoclones requires --reference.')
    if args.extract_sequences and not args.reference:
        parser.error('--extract-sequences requires --reference.')

    os.makedirs(args.output, exist_ok=True)
    stats_dir = os.path.join(args.output, 'statistical_tests')

    if use_mode1:
        import pandas as pd
        from modules.gdtools_runner import (
            find_gd_files, find_summary_jsons, run_gdtools_compare, run_gdtools_count,
        )

        gd_files = find_gd_files(args.samples_dir)
        if not gd_files:
            sys.exit(f'No output.gd files found in {args.samples_dir}')
        json_files = find_summary_jsons(args.samples_dir)
        compare_path = os.path.join(args.output, 'mutation_data.csv')
        run_gdtools_compare(args.gdtools, args.reference, gd_files, compare_path)
        input_file = compare_path
        header_row = 0

        from modules.dnds_test import compute_dnds_test
        os.makedirs(stats_dir, exist_ok=True)
        count_path = os.path.join(stats_dir, 'base_substitution_counts.csv')
        run_gdtools_count(args.gdtools, args.reference, gd_files, count_path)
        dnds_df = compute_dnds_test(pd.read_csv(count_path))
        dnds_file = os.path.join(stats_dir, 'dnds_test.csv')
        dnds_df.to_csv(dnds_file, index=False)
        print(f'  Saved: {dnds_file}')
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

    for threshold in args.threshold:
        df_filtered = frequency_filter(df_clean, threshold, ancestor)
        filtered_file = os.path.join(args.output, f'frequency_{threshold}.csv')
        df_filtered.to_csv(filtered_file, index=False)

    mutation_analysis(df_clean, ancestor, args.output)

    from modules.specificity import run_specificity_analysis
    run_specificity_analysis(
        df_clean, ancestor, stats_dir,
        permutations=args.specificity_permutations, seed=args.seed,
    )

    if args.pseudoclones:
        from modules.pseudoclones import run_pseudoclone_analysis
        run_pseudoclone_analysis(
            df_clean, ancestor, args.reference,
            os.path.join(args.output, 'pseudoclones'),
            bin_width=args.pseudoclone_bin_width, min_freq=args.pseudoclone_min_freq,
        )

    if args.extract_sequences:
        from modules.sequence_extraction import run_sequence_extraction
        run_sequence_extraction(
            df_clean, args.reference, os.path.join(args.output, 'gene_sequences'),
        )

    if args.plot:
        plots_dir = os.path.join(args.output, 'plots')
        os.makedirs(plots_dir, exist_ok=True)
        summary_file = os.path.join(args.output, 'mutation_summary.csv')
        if 'bubble' in args.plot:
            from modules.plotting import plot_mutation_bubble
            plot_mutation_bubble(summary_file,
                                 output_file=os.path.join(plots_dir, 'bubble_plot.png'),
                                 show=False)
        if 'spectrum' in args.plot:
            from modules.plotting import plot_mutation_spectrum
            plot_mutation_spectrum(summary_file,
                                   output_file=os.path.join(plots_dir, 'mutation_spectrum.png'),
                                   show=False)
        if 'parallel' in args.plot:
            from modules.plotting import plot_parallel_mutation_heatmap
            gene_parallel_file = os.path.join(args.output, 'gene_parallel_mutations.csv')
            if os.path.exists(gene_parallel_file):
                plot_parallel_mutation_heatmap(
                    gene_parallel_file,
                    output_file=os.path.join(plots_dir, 'parallel_mutation_heatmap.png'),
                    show=False,
                )
            else:
                print(f'Warning: {gene_parallel_file} not found, cannot generate parallel heatmap.', file=sys.stderr)
        if 'allele' in args.plot:
            from modules.plotting import plot_allele_distribution
            plot_allele_distribution(clean_file,
                                     output_file=os.path.join(plots_dir, 'allele_distribution.png'),
                                     show=False)
        if 'trajectory' in args.plot:
            from modules.plotting import plot_time_trajectory
            plot_time_trajectory(summary_file,
                                 output_file=os.path.join(plots_dir, 'time_trajectory.png'),
                                 show=False)
        if 'genome' in args.plot:
            from modules.plotting import plot_zoomed_genome
            if args.reference:
                plot_zoomed_genome(
                    clean_file, args.reference, ancestor=ancestor,
                    min_strains=args.genome_min_strains,
                    output_file=os.path.join(plots_dir, 'zoomed_genome.png'),
                    show=False,
                )
            else:
                print('Warning: --reference is required for the zoomed genome plot; skipping.', file=sys.stderr)
        print(f"Plots saved to {plots_dir}")

    print("Analysis complete.")


def run_simulation(args, parser, kind):
    import pandas as pd

    if args.replicates < 1:
        parser.error('--replicates must be >= 1.')
    if args.replicates > 50_000:
        print(f'Warning: --replicates={args.replicates}: runtime may be very long.', file=sys.stderr)

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
            companion_path=args.companion,
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
            companion_path=args.companion,
            n_replicates=args.replicates,
            seed=args.seed,
        )
    else:
        parser.error(f'Unknown simulation type: {kind}')


HELP_COMMANDS = {'-h', '--help'}


def main():
    argv = sys.argv[1:]
    # Backwards compatibility
    if argv and argv[0] not in SUBCOMMANDS and argv[0] not in HELP_COMMANDS:
        argv = ['process'] + argv

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == 'process':
        run_process(args, parser)
    elif args.command in ('simulate-dnds', 'simulate-parallel'):
        run_simulation(args, parser, args.command)


if __name__ == '__main__':
    main()
