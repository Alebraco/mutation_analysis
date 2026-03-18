#!/usr/bin/env python3
"""
Mutation Analysis Pipeline
Process mutation data from breseq output files and generate statistics

Two input modes:
    Mode 1 — from raw breseq sample directories:
        python main.py <ancestor_name> --samples-dir <dir> --reference <ref.gbk>
    Mode 2 — from a pre-existing file:
        python main.py <ancestor_name> <input_file>
"""

import argparse
import os
import sys
from modules.data_loader import load_and_filter
from modules.statistics import calculate_basic_stats, frequency_filter
from modules.analysis import mutation_analysis

def main():
    parser = argparse.ArgumentParser(description='Analyze mutation data from breseq output files.')

    parser.add_argument('ancestor', help='Ancestor name (e.g., "KZ_19")')
    parser.add_argument('--output', default='output_files', 
                       help='Output directory (default: output_files/)')
    parser.add_argument('--threshold', nargs='+', type=float, default=[],
                        help='Mutation filtering based on frequency threshold (default: None)')

    # Mode 1: raw breseq sample directory
    mode1 = parser.add_argument_group('Mode 1: raw breseq output')
    mode1.add_argument('--samples-dir', metavar='DIR',
                        help='Directory containing breseq sample folders (data/output.gd)')
    mode1.add_argument('--reference', 
                       help='Reference genome file for gdtools COMPARE')
    mode1.add_argument('--gdtools', default='gdtools',
                       help='Path to gdtools executable (default: gdtools)')
    
    # Mode 2: pre-processed mutation table
    mode2 = parser.add_argument_group('Mode 2: pre-processed mutation table')
    mode2.add_argument('input_file', nargs='?',
                       help='Input file path (.xlsx, .tsv, .csv)')
    mode2.add_argument('--header-row', type=int, default=0,
                          help='Header row index for input file (default: 0)')

    args = parser.parse_args()

    use_mode1 = bool(args.samples_dir or args.reference)
    use_mode2 = bool(args.input_file)

    if use_mode1 and use_mode2:
        parser.error('Cannot use --samples-dir/--reference together with --input-file.')
    if not use_mode1 and not use_mode2:
        parser.error('Provide either --input-file (processed file) or both --samples-dir and --reference (raw breseq outputs).')
    if use_mode1 and not (args.samples_dir and args.reference):
        parser.error('Both --samples-dir and --reference are required.')

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    if use_mode1:
        from modules.gdtools_runner import find_gd_files, find_summary_jsons, run_gdtools_compare
        gd_files = find_gd_files(args.samples_dir)
        if not gd_files:
            sys.exit(f'No output.gd files found in {args.samples_dir}')
        json_files = find_summary_jsons(args.samples_dir)
        tsv_path = os.path.join(args.output, 'mutation_data.tsv')
        run_gdtools_compare(args.gdtools, args.reference, gd_files, tsv_path)
        input_file = tsv_path
        header_row = 0
    else:
        input_file = args.input_file
        header_row = args.header_row
        json_files = None

    # Load and filter data
    df_clean, question_df = load_and_filter(input_file, args.ancestor, header_row)

    # Save cleaned data
    clean_file = os.path.join(args.output, 'cleaned_data.csv')
    df_clean.to_csv(clean_file, index=False)
    print(f'Saved cleaned data: {clean_file}')
    print(f'Cleaned data contains {len(df_clean)} rows.')

    # Save low coverage rows if any
    if question_df is not None:
        question_file = os.path.join(args.output, 'low_coverage_rows.csv')
        question_df.to_csv(question_file, index=False)
        print(f'Saved low coverage rows: {question_file}')

    # Calculate basic statistics
    calculate_basic_stats(df_clean, args.ancestor, args.output, json_files=json_files)

    # Filter mutations based on frequency threshold
    for threshold in args.threshold:
        df_filtered = frequency_filter(df_clean, threshold, args.ancestor)
        filtered_file = os.path.join(args.output, f'frequency_{threshold}.csv')
        df_filtered.to_csv(filtered_file, index=False)

    # Run mutations analysis
    mutation_analysis(df_clean, args.ancestor, args.output)

    print("Analysis complete.")

if __name__ == '__main__':
    main()