import pandas as pd
from .utils import codon_table

def classify_mutation(annotation):
    '''
    Classify mutation type based on annotation.
    '''
    if pd.isna(annotation):
        return 'unknown'
    annotation = str(annotation).lower()

    # Nonsense mutations
    if '*' in annotation:
        return 'nonsense'

    # Noncoding mutations
    elif 'noncoding' in annotation:
        return 'noncoding'
    elif 'intergenic' in annotation:
        return 'intergenic'
    elif 'pseudogene' in annotation:
        return 'pseudogene'

    # Coding mutations
    # Indels
    elif 'coding' in annotation:
        return 'nonsynonymous'

    # Point mutations
    elif '→' in annotation:
        cods = annotation.split('→')

        oldcodon = cods[0][-3:].upper()
        newcodon = cods[1][:3].upper()

        old_aa = codon_table[oldcodon]
        new_aa = codon_table[newcodon]

        if old_aa == new_aa:
            return 'synonymous'
        elif new_aa == '*':
            return 'nonsense'
        else:
            return 'nonsynonymous'

    else:
        print(f'Warning: Unrecognized annotation format: {annotation}')
        return 'unknown'


def classify_from_gdtools_json(snp_type, gene_position):
    '''
    Classify mutation type from gdtools JSON output fields (Mode 1 only).
    '''
    valid_types = {'nonsense', 'noncoding', 'intergenic', 'pseudogene', 'synonymous', 'nonsynonymous'}
    if snp_type in valid_types:
        return snp_type

    text = str(gene_position).lower() if pd.notna(gene_position) else ''

    if 'noncoding' in text:
        return 'noncoding'
    elif 'intergenic' in text:
        return 'intergenic'
    elif 'pseudogene' in text:
        return 'pseudogene'
    elif 'coding' in text:
        return 'nonsynonymous'
    else:
        return 'unknown'
