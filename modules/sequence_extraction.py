'''
Extract nucleotide and protein sequences for mutated genes.

Reuses the same pattern as modules/simulation/dnds.py
position -> GFF feature -> CDS -> translate.
'''

import os

import pandas as pd

from .simulation.reference_loader import load_reference
from .simulation.sequence_utils import rev_comp, translate


def run_sequence_extraction(df, reference_path, output_dir, companion_path=None):
    seq, gff = load_reference(reference_path, companion_path)

    genes = {}
    for _, row in df.iterrows():
        contig = str(row['seq_id'])
        pos = pd.to_numeric(row['position'], errors='coerce')
        if pd.isna(pos):
            continue
        for feature in gff.get(contig, {}).get(int(pos), []):
            name, start, end, strand, gene_type = feature
            entry = genes.setdefault(name, {'feature': (contig, start, end, strand, gene_type), 'hits': []})
            entry['hits'].append((int(pos), str(row.get('gene', '')), str(row.get('description', ''))))

    if not genes:
        print('Sequence extraction skipped: no mutations overlapped an annotated gene.')
        return

    os.makedirs(output_dir, exist_ok=True)
    nt_path = os.path.join(output_dir, 'gene_sequences_nt.fasta')
    aa_path = os.path.join(output_dir, 'gene_sequences_aa.fasta')
    manifest = []

    with open(nt_path, 'w') as nt_h, open(aa_path, 'w') as aa_h:
        for name, entry in genes.items():
            contig, start, end, strand, gene_type = entry['feature']
            nt = seq[contig][start - 1:end]
            if strand == '-':
                nt = rev_comp(nt)
            nt_h.write(f'>{name}\n{nt}\n')

            aa = translate(nt) if gene_type == 'CDS' and len(nt) % 3 == 0 else ''
            if aa:
                aa_h.write(f'>{name}\n{aa}\n')

            positions = sorted({h[0] for h in entry['hits']})
            manifest.append({
                'feature': name, 'contig': contig, 'start': start, 'end': end,
                'strand': strand, 'gene_type': gene_type,
                'gene_name': next((h[1] for h in entry['hits'] if h[1]), ''),
                'description': next((h[2] for h in entry['hits'] if h[2]), ''),
                'protein_length': len(aa),
                'n_mutations': len(entry['hits']),
                'positions': ';'.join(str(p) for p in positions),
            })

    pd.DataFrame(manifest).to_csv(os.path.join(output_dir, 'gene_sequences_manifest.csv'), index=False)
    print(f'Sequence extraction: wrote {len(genes)} genes to {nt_path} / {aa_path}')
