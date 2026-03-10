import pandas as pd
import numpy as np
import re
from utils import get_strain_columns
from data_cleaner import clean_text
from mutation_classifier import classify_mutation

def load_and_filter(input_file, ancestor, header_row=0):
    '''
    Load data and filter invalid rows
    '''
    print(f'Reading input file: {input_file}')
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
    
    # Remove triangles from strain columns only
    for col in strain_cols:
        df[col] = df[col].astype(str).str.replace('Δ', '')

    valid_rows = []
    question_rows = []

    for index, row in df.iterrows():
        row_string = str(row.values)

        # Exclude ancestor mutations
        if pd.notna(df.loc[index, ancestor]):
            continue

        # Exclude rows with low coverage
        if '?' in row_string:
            question_rows.append(row)
            continue

        # Append valid rows only
        valid_rows.append(row)

    print(f'Excluded {len(df) - len(valid_rows)} rows')

    # Save valid rows
    df_clean = pd.DataFrame(valid_rows)
    # Save low coverage rows
    question_df = pd.DataFrame(question_rows) if question_rows else None

    print('Cleaning nonstandard characters.')
    for col in df_clean.columns:
        df_clean[col] = df_clean[col].apply(clean_text)

    print('Classifying mutations.')
    df_clean['mutation_type'] = df_clean['annotation'].apply(classify_mutation)

    return df_clean, question_df