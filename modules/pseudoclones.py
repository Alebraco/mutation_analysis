'''
Pseudoclone generation for clonal deconvolution / contamination checks.
For each sample, SNPs are grouped into frequency bins (deciles by default).
Each (sample, bin) is written as one genome, the reference with the SNPs applied.

SNP-only: other mutations are ignored.
Adapted from Dr. Louis-Marie Bobay's script. 
'''

import os
import re

import pandas as pd

from .utils import get_strain_columns
from .simulation.reference_loader import load_sequences

SNP = re.compile(r'^([ACGT])→([ACGT])$')


def _resolve(contigs, seq_id):
    '''Table seq_id -> reference contig: by name, else by trailing index (contig_1 -> 1st).'''
    if seq_id in contigs:
        return seq_id
    m = re.search(r'(\d+)$', seq_id)
    return contigs[int(m.group(1)) - 1] if m and int(m.group(1)) <= len(contigs) else None


def run_pseudoclone_analysis(df, ancestor, reference_path, output_dir,
                             companion_path=None, bin_width=0.1, min_freq=0.05):
    strains = get_strain_columns(df, ancestor)
    seqs = load_sequences(reference_path, companion_path)
    contigs = list(seqs)
    n_bins = max(1, round(1.0 / bin_width))

    # binned[strain][bin] -> list of (contig, pos, ref, alt)
    binned = {s: {} for s in strains}
    for _, row in df.iterrows():
        m = SNP.match(str(row['mutation']).strip())
        pos = pd.to_numeric(row['position'], errors='coerce')
        if not m or pd.isna(pos):
            continue
        rec = (_resolve(contigs, str(row['seq_id'])), int(pos), m.group(1), m.group(2))
        for s in strains:
            freq = pd.to_numeric(row[s], errors='coerce')
            if pd.notna(freq) and freq >= min_freq:
                # +1e-9 avoids float-floor errors, e.g. 0.3/0.1 == 2.9999996
                b = min(int(freq / bin_width + 1e-9), n_bins - 1)
                binned[s].setdefault(b, []).append(rec)

    if not any(binned.values()):
        print('Pseudoclone generation skipped: no single-base substitutions with usable frequencies.')
        return

    os.makedirs(output_dir, exist_ok=True)
    fasta = os.path.join(output_dir, 'pseudoclones.fa')
    manifest, dropped = [], 0
    with open(fasta, 'w') as h:
        for s in strains:
            for d in sorted(binned[s]):
                genome = {c: bytearray(v, 'ascii') for c, v in seqs.items()}
                applied = skipped = 0
                for contig, pos, ref, alt in binned[s][d]:
                    ba = genome.get(contig)
                    if ba is not None and 1 <= pos <= len(ba) and ba[pos - 1] == ord(ref):
                        ba[pos - 1] = ord(alt)
                        applied += 1
                    else:
                        skipped += 1
                dropped += skipped
                if not applied:
                    continue
                lo, hi = round(d * bin_width, 4), round(min((d + 1) * bin_width, 1.0), 4)
                name = f'{s}__f{int(lo * 100)}-{int(hi * 100)}'
                h.write(f'>{name}\n' + ''.join(g.decode() for g in genome.values()) + '\n')
                manifest.append({'pseudoclone': name, 'sample': s, 'freq_low': lo,
                                 'freq_high': hi, 'snps_applied': applied,
                                 'snps_skipped': skipped})

    pd.DataFrame(manifest).to_csv(os.path.join(output_dir, 'pseudoclone_manifest.csv'), index=False)
    print(f'Pseudoclones: wrote {len(manifest)} genomes to {fasta}')
    if dropped:
        print(f'  Note: {dropped} SNP placement(s) skipped (reference-base mismatch or unmapped contig).')
