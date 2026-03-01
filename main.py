#!/usr/bin/env python3
"""
Mutation Analysis Pipeline
Process mutation data from Excel files and generate statistics
"""

import argparse
import os
import pandas as pd
from data_loader import load_and_filter
from statistics import calculate_basic_stats, frequency_filter
from analysis import mutation_analysis

def main():
    parser = argparse.ArgumentParser(description='Analyze mutation data from breseq output files in Excel format.')
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('ancestor', help='Ancestor name (e.g., "KZ_19")')
    parser.add_argument('--header-row', type=int, default=1, 
                       help='Header row index (default: 1)')
    parser.add_argument('--output-dir', default='.', 
                       help='Output directory (default: current directory)')
    
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load and filter data
    df_clean, question_df = load_and_filter(
        args.input_file, 
        args.ancestor, 
        args.header_row
    )

    # Save cleaned data
    clean_file = os.path.join(args.output_dir, 'cleaned_data.xlsx')
    df_clean.to_excel(clean_file, index=False)
    print(f'Saved cleaned data: {clean_file}')

    # Save low coverage rows if any
    if question_df is not None:
        question_file = os.path.join(args.output_dir, 'low_coverage_rows.xlsx')
        question_df.to_excel(question_file, index=False)
        print(f'Saved low coverage rows: {question_file}')

    # Calculate basic statistics
    calculate_basic_stats(df_clean, args.ancestor, args.output_dir)

    # Filter mutations based on frequency threshold
    frequency_threshold = [0.25, 0.5, 0.75, 1.0]
    for threshold in frequency_threshold:
        df_filtered = frequency_filter(df_clean, threshold, args.ancestor)
        filtered_file = os.path.join(args.output_dir, f'frequency_{threshold}.xlsx')
        df_filtered.to_excel(filtered_file, index=False)

    # Run mutations analysis
    mutation_analysis(df_clean, args.ancestor, args.output_dir)

    print("Analysis complete.")

if __name__ == '__main__':
    main()