from ..utils import codon_table

COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N',
}


def rev_comp(seq: str) -> str:
    seq = seq.upper()
    out = []
    for base in reversed(seq):
        comp = COMPLEMENT.get(base)
        if comp is None:
            raise ValueError(f'Unexpected nucleotide symbol: {base!r}')
        out.append(comp)
    return ''.join(out)


def translate(seq: str) -> str:
    seq = seq.upper()
    aas = []
    i = 0
    while i <= len(seq) - 3:
        codon = seq[i:i + 3]
        aa = codon_table.get(codon, 'X')
        aas.append(aa)
        if aa == '*':
            break
        i += 3
    return ''.join(aas)
