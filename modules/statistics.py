import os
import json
import pandas as pd
import numpy as np
from .utils import get_strain_columns

def load_coverage_averages(json_files):
    '''
    Read each sample's summary.json and extract average coverage depth.
    Extracts: summary["references"]["reference"]["seq_id"]["coverage_average"]
    Returns {sample_name: average} dict. Returns "NA" per sample on error.
    '''
    
    coverage = {}
    for sample, json_path in json_files.items():
        try:
            with open(json_path) as f:
                summary = json.load(f)
            
            reference = summary.get('references', {}).get('reference', {})
            cov_total = 0.0
            seq_number = 0

            for seq_id, seq_data in reference.items():
                avg = seq_data.get('coverage_average')
                if avg is not None:
                    cov_total += float(avg)
                    seq_number += 1
                else:
                    print(f'No coverage_average found in {json_path} for sample {sample}')

            if seq_number > 0:
                coverage[sample] = round(cov_total / seq_number, 4)
            else:
                coverage[sample] = pd.NA

        except Exception as e:
            print(f'Error loading coverage from {json_path} for sample {sample}: {e}')
            coverage[sample] = pd.NA
    return coverage


def calculate_basic_stats(df, ancestor, output_dir='.', json_files=None):
    '''
    Calculate Basic Statistics per line:
    Mutation Classification, Number of Mutations, Average Frequency, Average Coverage
    '''

    os.makedirs(output_dir, exist_ok=True)
    strain_cols = get_strain_columns(df, ancestor)

    results = []

    for strain in strain_cols:
        # Get only rows with mutations in this strain
        strain_df = df[df[strain].notna()]

        type_counts = strain_df['mutation_type'].value_counts()
        total = len(strain_df)
        avg_frequency = pd.to_numeric(strain_df[strain], errors='coerce').mean()

        # Calculate mutation proportions
        proportions = (type_counts / total) if total > 0 else type_counts
        proportions = proportions.round(6)

        results.append({
            'Line': strain,
            'Nonsynonymous': proportions.get('nonsynonymous', 0),
            'Synonymous': proportions.get('synonymous', 0),
            'Intergenic': proportions.get('intergenic', 0),
            'Nonsense': proportions.get('nonsense', 0),
            'Noncoding': proportions.get('noncoding', 0),
            'Pseudogene': proportions.get('pseudogene', 0),
            'Unknown': proportions.get('unknown', 0),
            'Total Mutations': total,
            'Average Frequency': avg_frequency,
        })

    summary_df = pd.DataFrame(results)

    if json_files is not None:
        coverage = load_coverage_averages(json_files)
        coverage = {k.lower(): v for k, v in coverage.items()}
        summary_df['Average Coverage'] = summary_df['Line'].map(coverage)
    summary_file = os.path.join(output_dir, 'mutation_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    print(f'  Saved: {summary_file}')
    return summary_df

def frequency_filter(df, min_frequency, ancestor):
    '''
    Filter mutations based on minimum frequency threshold
    Keeps entire row if at least one strain meets the threshold
    '''
    strain_cols = get_strain_columns(df, ancestor)
    values = df[strain_cols].apply(pd.to_numeric, errors='coerce')

    keep_rows = (values >= min_frequency).any(axis=1)
    df_filtered = df[keep_rows]
    print(f'{len(df_filtered)} mutations meet frequency threshold: {min_frequency}')

    return df_filtered