import pandas as pd
from utils import codon_table

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
    else:
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