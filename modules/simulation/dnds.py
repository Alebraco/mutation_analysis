import csv
import os
import random
from typing import Optional

import numpy as np
import pandas as pd

from .mutation_model import (
    build_contig_pointers,
    draw_replicate,
    estimate_mutation_model,
)
from .reference_loader import load_reference
from .sequence_utils import rev_comp, translate

_CATEGORIES = ('intergenic', 'RNA', 'S', 'NS', 'STOP')


def run_dnds_simulation(
    df_clean: pd.DataFrame,
    ancestor: str,
    reference_path: str,
    output_dir: str,
    file_stem: str,
    gff_path: Optional[str] = None,
    n_replicates: int = 1000,
    freq_threshold: float = 0.05,
    seed: Optional[int] = None,
) -> None:
    '''
    Simulate `n_replicates` neutral-mutation replicates and write the expected
    per-category fractions to <output_dir>/expected/expdNdS_<stem>.csv and
    avg_expdNdS_<stem>.csv.
    '''
    seq, gff = load_reference(reference_path, gff_path=gff_path)
    alpha, outcome, load = estimate_mutation_model(df_clean, ancestor, freq_threshold)
    if not load:
        raise RuntimeError(
            f'No strain has mutations at freq >= {freq_threshold}; nothing to simulate.'
        )
    pointer, total_length = build_contig_pointers(seq)
    rng = random.Random(seed)

    # typing[st][category] -> list of per-replicate counts
    typing: dict[str, dict[str, list[float]]] = {
        st: {cat: [] for cat in _CATEGORIES} for st in load
    }

    for rep_idx in range(1, n_replicates + 1):
        if rep_idx % 100 == 0:
            print(f'[simulate-dnds] replicate {rep_idx}/{n_replicates}')

        for st in load:
            picks = draw_replicate(seq, pointer, total_length, alpha, outcome, load[st], rng)
            counts = _classify_replicate(seq, gff, picks)
            for cat in _CATEGORIES:
                typing[st][cat].append(float(counts[cat]))

    os.makedirs(os.path.join(output_dir, 'expected'), exist_ok=True)
    per_strain_path = os.path.join(output_dir, 'expected', f'expdNdS_{file_stem}.csv')
    avg_path = os.path.join(output_dir, 'expected', f'avg_expdNdS_{file_stem}.csv')

    averages: dict[str, dict[str, float]] = {
        st: {cat: float(np.mean(typing[st][cat])) for cat in _CATEGORIES}
        for st in typing
    }

    with open(per_strain_path, 'w', newline='') as h:
        writer = csv.writer(h)
        writer.writerow(['strain', 'category', 'expected_percent'])
        for st, cats in averages.items():
            total = sum(cats.values())
            if total == 0:
                continue
            coeff = 100.0 / total
            for cat, mean in cats.items():
                writer.writerow([st, cat, coeff * mean])
    print(f'  Saved: {per_strain_path}')

    super_average: dict[str, list[float]] = {cat: [] for cat in _CATEGORIES}
    for cats in averages.values():
        for cat in _CATEGORIES:
            super_average[cat].append(cats[cat])

    total = sum(float(np.mean(super_average[cat])) for cat in _CATEGORIES)
    with open(avg_path, 'w', newline='') as h:
        writer = csv.writer(h)
        writer.writerow(['category', 'expected_percent'])
        if total > 0:
            coeff = 100.0 / total
            for cat in _CATEGORIES:
                writer.writerow([cat, coeff * float(np.mean(super_average[cat]))])
    print(f'  Saved: {avg_path}')


def _classify_replicate(
    seq: dict[str, str],
    gff: dict[str, dict[int, list]],
    picks: list[tuple[str, int, str, str]],
) -> dict[str, int]:
    '''Apply the simulated mutations and tally the categories for a single replicate.'''
    counts = {cat: 0 for cat in _CATEGORIES}

    mutated = {c: list(s) for c, s in seq.items()}
    for contig, pos, n1, n2 in picks:
        idx = pos - 1
        if mutated[contig][idx] != n1:
            raise RuntimeError(
                f'Ancestral nucleotide mismatch at {contig}:{pos} '
                f'(expected {n1!r}, found {mutated[contig][idx]!r}).'
            )
        mutated[contig][idx] = n2
    mutated_seqs = {c: ''.join(chars) for c, chars in mutated.items()}

    for contig, pos, _n1, _n2 in picks:
        features = gff.get(contig, {}).get(pos)
        if not features:
            counts['intergenic'] += 1
            continue

        for feature in features:
            _name, start, end, strand, gene_type = feature
            if gene_type != 'CDS':
                counts['RNA'] += 1
                break

            cds_new = mutated_seqs[contig][start - 1:end]
            cds_old = seq[contig][start - 1:end]
            if strand == '-':
                cds_new = rev_comp(cds_new)
                cds_old = rev_comp(cds_old)
            if len(cds_old) % 3 != 0:
                # Skip CDS whose length isn't a multiple of 3 (matches source behavior).
                continue
            new_prot = translate(cds_new)
            old_prot = translate(cds_old)
            if len(new_prot) < len(old_prot):
                counts['STOP'] += 1
            elif new_prot == old_prot:
                counts['S'] += 1
            else:
                counts['NS'] += 1
            break
    return counts
