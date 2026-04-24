import csv
import os
import random
from itertools import combinations

import pandas as pd

from .mutation_model import (
    build_contig_pointers,
    draw_replicate,
    estimate_mutation_model,
)
from .reference_loader import load_reference


def run_parallel_simulation(
    df_clean: pd.DataFrame,
    ancestor: str,
    reference_path: str,
    output_dir: str,
    file_stem: str,
    companion_path: str | None = None,
    n_replicates: int = 10000,
    freq_threshold: float = 0.05,
    seed: int | None = None,
) -> None:
    '''
    Simulate `n_replicates` neutral mutation replicates and write the expected
    number of parallel sites and parallel genes per strain pair.
    '''

    seq, gff = load_reference(reference_path, companion_path=companion_path)
    alpha, outcome, load = estimate_mutation_model(df_clean, ancestor, freq_threshold)
    if len(load) < 2:
        raise RuntimeError(
            'Need at least two strains with mutations to simulate parallel mutations.'
        )
    pointer, total_length = build_contig_pointers(seq)
    rng = random.Random(seed)

    strains = sorted(load.keys())
    pair_list = list(combinations(strains, 2))

    site_totals: dict[tuple[str, str], int] = {p: 0 for p in pair_list}
    gene_totals: dict[tuple[str, str], int] = {p: 0 for p in pair_list}

    for rep_idx in range(1, n_replicates + 1):
        if rep_idx % 100 == 0:
            print(f'[simulate-parallel] replicate {rep_idx}/{n_replicates}')

        replicate_picks: dict[str, list[tuple[str, int, str, str]]] = {}
        for st in strains:
            replicate_picks[st] = draw_replicate(
                seq, pointer, total_length, alpha, outcome, load[st], rng
            )

        by_site: dict[str, set[tuple[str, int]]] = {}
        by_gene: dict[str, set[tuple[str, str]]] = {}
        for st, picks in replicate_picks.items():
            site_set: set[tuple[str, int]] = set()
            gene_set: set[tuple[str, str]] = set()
            for contig, pos, n1, n2 in picks:
                site_set.add((contig, pos))
                features = gff.get(contig, {}).get(pos)
                if features:
                    for feature in features:
                        gene_set.add((contig, feature[0]))
            by_site[st] = site_set
            by_gene[st] = gene_set

        for st1, st2 in pair_list:
            site_totals[(st1, st2)] += len(by_site[st1] & by_site[st2])
            gene_totals[(st1, st2)] += len(by_gene[st1] & by_gene[st2])

    os.makedirs(os.path.join(output_dir, 'expected'), exist_ok=True)
    sites_path = os.path.join(output_dir, 'expected', f'parallel_sites_{file_stem}.csv')
    genes_path = os.path.join(output_dir, 'expected', f'parallel_genes_{file_stem}.csv')

    with open(sites_path, 'w', newline='') as h:
        writer = csv.writer(h)
        writer.writerow(['strain_1', 'strain_2', 'expected_parallel_sites'])
        for (st1, st2), total in site_totals.items():
            writer.writerow([st1, st2, total / n_replicates])
    print(f'  Saved: {sites_path}')

    with open(genes_path, 'w', newline='') as h:
        writer = csv.writer(h)
        writer.writerow(['strain_1', 'strain_2', 'expected_parallel_genes'])
        for (st1, st2), total in gene_totals.items():
            writer.writerow([st1, st2, total / n_replicates])
    print(f'  Saved: {genes_path}')
