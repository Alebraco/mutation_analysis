'''
Treatment-specificity analysis: per-gene Fisher's exact test, Dice 
and rank tests. Methods adapted from barricklab/breseq-ext-specificity
(Maddamsetti & Deatherage). Treatments come from the
`group` label in each strain name (utils.parse_line_label)
'''

import os
import random
from collections import Counter
from itertools import combinations

import pandas as pd
from scipy.stats import mannwhitneyu, kruskal, fisher_exact

from .utils import get_strain_columns, parse_line_label


def _assign_treatments(strain_cols):
    '''Map strains to treatments and drop strains whose names lack a group token.'''
    strain_to_treatment, treatment_to_strains = {}, {}
    for strain in strain_cols:
        _, group, _ = parse_line_label(strain)
        if group is None:
            continue
        strain_to_treatment[strain] = group
        treatment_to_strains.setdefault(group, []).append(strain)
    return strain_to_treatment, treatment_to_strains


def _valid_gene(gene):
    '''Single-gene name, or None for intergenic/multi-gene (include '/').'''
    if pd.isna(gene):
        return None
    name = str(gene).strip()
    if not name or '/' in name or name in ('-', '–', '.'):
        return None
    return name


def _gene_hit_sets(df, strain_cols):
    '''Per strain, the set of single genes it carries a mutation in.'''
    gene_sets = {strain: set() for strain in strain_cols}
    genes = df['gene'].tolist()
    for strain in strain_cols:
        for gene, present in zip(genes, df[strain].notna().tolist()):
            valid = _valid_gene(gene) if present else None
            if valid is not None:
                gene_sets[strain].add(valid)
    return gene_sets


def gene_specificity_table(df, gene_sets, strain_to_treatment, treatments):
    '''Two-sided Fisher per gene, scored on genomes hit.'''
    total_genomes = len(strain_to_treatment)
    genomes_per_gene = {}
    for strain, genes in gene_sets.items():
        if strain not in strain_to_treatment:
            continue
        for gene in genes:
            genomes_per_gene.setdefault(gene, []).append(strain)

    descriptions = df.drop_duplicates('gene').set_index('gene')['description'].to_dict()

    rows = []
    for gene, hits in genomes_per_gene.items():
        if len(hits) < 2:
            continue
        counts = Counter(strain_to_treatment[s] for s in hits)
        top, top_count = counts.most_common(1)[0]
        n_top = len(treatments[top])
        s_c1, f_c1 = top_count, n_top - top_count
        s_c2 = len(hits) - top_count
        f_c2 = (total_genomes - n_top) - s_c2
        row = {'gene': gene, 'description': descriptions.get(gene, ''),
               'genomes_hit': len(hits), 'top_treatment': top,
               'p_value': fisher_exact([[s_c1, f_c1], [s_c2, f_c2]], alternative='two-sided')[1]}
        for t in treatments:
            row[f'hits_{t}'] = counts.get(t, 0)
        rows.append(row)

    if not rows:
        return None
    return pd.DataFrame(rows).sort_values(
        ['p_value', 'genomes_hit'], ascending=[True, False]).reset_index(drop=True)


def _dice(x, y):
    '''Dice similarity of two sets.'''
    denom = len(x) + len(y)
    return 2.0 * len(x & y) / denom if denom else None


def _within_between(strains, gene_sets, treatment_of):
    within, between = [], []
    for a, b in combinations(strains, 2):
        d = _dice(gene_sets[a], gene_sets[b])
        if d is None:
            continue
        (within if treatment_of[a] == treatment_of[b] else between).append(d)
    return within, between


def dice_permutation_test(gene_sets, strain_to_treatment, n_permutations, rng):
    '''Within- and between-treatment Dice with a permutation p-value.'''
    strains = list(strain_to_treatment)
    labels = [strain_to_treatment[s] for s in strains]
    within, between = _within_between(strains, gene_sets, dict(zip(strains, labels)))
    if not within or not between:
        return None

    obs = sum(within) / len(within) - sum(between) / len(between)
    ge = 0
    for _ in range(n_permutations):
        shuffled = labels[:]
        rng.shuffle(shuffled)
        w, b = _within_between(strains, gene_sets, dict(zip(strains, shuffled)))
        if w and b and (sum(w) / len(w) - sum(b) / len(b)) >= obs:
            ge += 1

    return {'within_mean': sum(within) / len(within),
            'between_mean': sum(between) / len(between), 'difference': obs,
            'n_pairs_within': len(within), 'n_pairs_between': len(between),
            'permutations': n_permutations,
            'p_value': (1 + ge) / (n_permutations + 1) if n_permutations else float('nan')}


def rank_tests(burden, treatment_to_strains):
    '''Kruskal-Wallis + pairwise Mann-Whitney on mutation burden.'''
    treatments = sorted(treatment_to_strains)
    vals = {t: [burden[s] for s in treatment_to_strains[t]] for t in treatments}

    try:
        h, p = kruskal(*[vals[t] for t in treatments])
    except ValueError:
        h, p = float('nan'), float('nan')
    results = [{'test': 'kruskal-wallis', 'groups': ' + '.join(treatments),
                'statistic': h, 'p_value': p}]

    for a, b in combinations(treatments, 2):
        try:
            u, p = mannwhitneyu(vals[a], vals[b], alternative='two-sided')
        except ValueError:
            u, p = float('nan'), float('nan')
        results.append({'test': 'mann-whitney', 'groups': f'{a} vs {b}',
                        'statistic': u, 'p_value': p})
    return results


def run_specificity_analysis(df, ancestor, output_dir='.', permutations=10000, seed=None):
    '''Run the specificity analysis.'''
    strain_cols = get_strain_columns(df, ancestor)
    strain_to_treatment, treatment_to_strains = _assign_treatments(strain_cols)

    if len(treatment_to_strains) < 2:
        print(f'Specificity analysis skipped: found {len(treatment_to_strains)} '
              f'treatment in strain names (need 2 or more, e.g. "d120-me1").')
        return None

    os.makedirs(output_dir, exist_ok=True)

    print(f'Treatment specificity: {len(treatment_to_strains)} treatments '
          f'({", ".join(f"{t}={len(s)}" for t, s in sorted(treatment_to_strains.items()))}).')

    gene_sets = _gene_hit_sets(df, strain_cols)

    gene_table = gene_specificity_table(df, gene_sets, strain_to_treatment, treatment_to_strains)
    if gene_table is not None:
        gene_file = os.path.join(output_dir, 'gene_specificity.csv')
        gene_table.to_csv(gene_file, index=False)
        print(f'  Saved: {gene_file}')

    rng = random.Random(seed)
    dice_result = dice_permutation_test(gene_sets, strain_to_treatment, permutations, rng)
    if dice_result is not None:
        dice_file = os.path.join(output_dir, 'treatment_dice.csv')
        pd.DataFrame([dice_result]).to_csv(dice_file, index=False)
        print(f'  Saved: {dice_file} (within={dice_result["within_mean"]:.3f}, '
              f'between={dice_result["between_mean"]:.3f}, p={dice_result["p_value"]:.4f})')

    burden = {s: int(df[s].notna().sum()) for s in strain_cols}
    rank_result = rank_tests(burden, treatment_to_strains)
    rank_file = os.path.join(output_dir, 'treatment_rank_tests.csv')
    pd.DataFrame(rank_result).to_csv(rank_file, index=False)
    print(f'  Saved: {rank_file}')

    return gene_table, dice_result, rank_result
