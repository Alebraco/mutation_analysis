import pandas as pd
import numpy as np
import re
import os 

from .utils import get_strain_columns
from .data_cleaner import clean_text
from .mutation_classifier import classify_mutation

def load_and_filter(input_file, ancestor, header_row=0):
    '''
    Load data and filter invalid rows.
    Accepts both .xlsx and tabular files.
    '''

    print(f'Reading input file: {input_file}')
    ext = os.path.splitext(input_file)[1].lower()
    if ext in ('.tsv', '.txt'):
        df = pd.read_csv(input_file, sep='\t', header=header_row)
    elif ext == '.csv':
        df = pd.read_csv(input_file, header=header_row)
    elif ext in ('.xlsx', '.xls'):
        df = pd.read_excel(input_file, header=header_row)

    # Standardize column names
    df.columns = (
                df.columns.astype(str)
                    .str.lower()                                # Make columns lower-case
                    .str.strip()                                # Remove trailing spaces
                    .str.replace(r'\s+', '_', regex=True)       # Convert single/multiple spaces into underscore
                    .str.replace('seqid', 'seq_id')             # Always use underscore for seq_id

    )

    ancestor = re.sub(r'\s+', '_', str(ancestor)    # Substitute spaces for underscores to match column
                      .lower().strip())             # Standardize ancestor name

    # Replace #ERROR! with NA
    df = df.replace('#ERROR!', pd.NA)

    if ancestor.endswith('.gd'):
        ancestor = ancestor.replace('.gd', '')

    # Get strain columns
    strain_cols = get_strain_columns(df, ancestor)
    
    # Remove triangles from strain columns only, preserving NaN values
    for col in strain_cols:
        df[col] = df[col].where(df[col].isna(), df[col].astype(str).str.replace('Δ', ''))

    valid_rows = []
    question_rows = []
    excluded_ancestor = 0
    excluded_low_coverage = 0

    for index, row in df.iterrows():
        row_string = str(row.values)

        # Exclude ancestor mutations
        if pd.notna(df.loc[index, ancestor]):
            excluded_ancestor += 1
            continue

        if '?' in row_string:
            question_rows.append(row)
            strain_values = [str(row[col]) for col in strain_cols]
            other_values = [val for val in strain_values if val not in ('?', 'nan', 'NA', 'None')]

            if not other_values:
                excluded_low_coverage += 1
                continue  # Skip rows with only '?' values

        # Append valid rows only
        valid_rows.append(row)

    if excluded_ancestor:
        print(f'Excluded {excluded_ancestor} ancestor mutation rows')
    if excluded_low_coverage:
        print(f'Excluded {excluded_low_coverage} low coverage rows (only "?" values)')

    # Save valid rows
    df_clean = pd.DataFrame(valid_rows)
    # Save low coverage rows
    question_df = pd.DataFrame(question_rows) if question_rows else None

    print('Cleaning nonstandard characters.')
    for col in df_clean.columns:
        df_clean[col] = df_clean[col].apply(clean_text)

    print('Classifying mutations.')
    df_clean['mutation_type'] = df_clean['annotation'].apply(classify_mutation)
    
    unknowns = df_clean[df_clean['mutation_type'] == 'unknown']
    if not unknowns.empty:
        print(f'{len(unknowns)} mutations classified as unknown. Consider reviewing their annotations for potential issues.')
        print(unknowns.value_counts().to_string())

    return df_clean, question_df