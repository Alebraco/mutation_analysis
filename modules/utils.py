#!/usr/bin/env python3

import re

import pandas as pd

# Codon table and descriptive columns for mutation analysis
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

descriptive_cols = ['seq_id', 'position', 'mutation', 'annotation',
                   'gene', 'description', 'mutation_type']

def get_strain_columns(df, ancestor):
    '''
    Get list of strain columns (non-descriptive columns only)
    '''
    # Start with descriptive columns and add ancestor column
    excluded_cols = descriptive_cols + [ancestor]

    # Strain columns are all columns that are not in descriptive_cols
    strain_cols = [col for col in df.columns if col not in excluded_cols]

    return strain_cols


def parse_line_label(line):
    '''
    Parse labels like:
    - sm-d120-me1-p
    - d120-me1
    - d180-me-6-p   (separator between treatment and replicate)
    Returns: (day, group, replicate)

    The `group` token is the treatment identifier used by the specificity
    analysis, so both plotting and specificity share this one parser.
    '''
    if pd.isna(line):
        return None, None, None

    match = re.search(r"([dD]\d+)[-_]([A-Za-z]+)[-_]?(\d+)", str(line))
    if not match:
        return None, None, None

    day = match.group(1)
    group = match.group(2)
    replicate = match.group(3)
    return day, group, replicate