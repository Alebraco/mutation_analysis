import os
import pandas as pd
import numpy as np
from utils import get_strain_columns

def calculate_basic_stats(df, ancestor, output_dir='.'):
    '''
    Calculate Basic Statistics per line:
    Mutation Classification, Number of Mutations, and Average Frequency
    '''
    os.makedirs(output_dir, exist_ok=True)
    strain_cols = get_strain_columns(df, ancestor)

    results = []

    for strain in strain_cols:
        # Get only rows with mutations in this strain
        strain_df = df[df[strain].notna()]

        type_counts = strain_df['mutation_type'].value_counts()
        total = len(strain_df)
        avg_frequency = pd.to_numeric(strain_df[strain]).mean()

        # Calculate mutation proportions
        proportions = (type_counts / total)
        proportions = proportions.round(6)

        results.append({
            'Line': strain,
            'Nonsynonymous': proportions.get('nonsynonymous', 0),
            'Synonymous': proportions.get('synonymous', 0),
            'Intergenic': proportions.get('intergenic', 0),
            'NonSense': proportions.get('nonsense', 0),
            'Noncoding': proportions.get('noncoding', 0),
            'Total Mutations': total,
            'Average Frequency': avg_frequency
        })

    summary_df = pd.DataFrame(results)
    summary_file = os.path.join(output_dir, 'mutation_summary.xlsx')
    summary_df.to_excel(summary_file, index=False)
    print(f'Saved mutation summary: {summary_file}')
    return summary_df