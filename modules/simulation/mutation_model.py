import random

import pandas as pd

from ..utils import get_strain_columns

Alpha = dict[str, float]
Outcome = dict[str, dict[str, float]]
Load = dict[str, int]

BASES = ('A', 'C', 'G', 'T')


def parse_point_mutation(raw: object) -> tuple[str, str] | None:
    '''Return (N1, N2) for a single-base substitution, else None.'''
    if not isinstance(raw, str):
        return None
    text = raw.strip()
    if '→' in text:
        left, right = text.split('→', 1)
    elif '->' in text:
        left, right = text.split('->', 1)
    elif len(text) == 3 and text[1] == '_': # e.g. "A_C"
        left, right = text[0], text[2]
    else:
        return None
    
    left = left.strip().upper()
    right = right.strip().upper()

    if len(left) != 1 or len(right) != 1:
        return None
    if left not in BASES or right not in BASES or left == right:
        return None
    return left, right


def force_frequency(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def estimate_mutation_model(
    df_clean: pd.DataFrame,
    ancestor: str,
    freq_threshold: float = 0.05,
) -> tuple[Alpha, Outcome, Load]:
    '''
    Derive the neutral-model inputs from the cleaned mutation table.

    Returns:
        alpha   — base composition of source in point mutations.
        outcome — substitution probabilities for each source base.
        load    — count of mutations (rows) meeting threshold frequency per strain.
    '''
    alpha: dict[str, int] = {b: 0 for b in BASES}
    outcome_counts: dict[str, dict[str, int]] = {
        b1: {b2: 0 for b2 in BASES if b2 != b1} for b1 in BASES
    }

    strain_cols = get_strain_columns(df_clean, ancestor)
    load: Load = {st: 0 for st in strain_cols}

    for _, row in df_clean.iterrows():
        parsed = parse_point_mutation(row.get('mutation'))
        if parsed is None:
            continue
        n1, n2 = parsed
        alpha[n1] += 1
        outcome_counts[n1][n2] += 1
        for st in strain_cols:
            freq = force_frequency(row.get(st))
            if freq is not None and freq >= freq_threshold:
                load[st] += 1

    total_alpha = sum(alpha.values())
    if total_alpha == 0:
        raise RuntimeError(
            'No base substitutions found in the mutation table.'
            'Cannot estimate the neutral model.'
        )
    alpha_norm: Alpha = {b: alpha[b] / total_alpha for b in BASES}

    outcome: Outcome = {}
    for n1 in BASES:
        row_total = sum(outcome_counts[n1].values())
        if row_total == 0:
            outcome[n1] = {n2: 1.0 / 3 for n2 in BASES if n2 != n1}
        else:
            outcome[n1] = {n2: outcome_counts[n1][n2] / row_total for n2 in outcome_counts[n1]}

    load = {st: ct for st, ct in load.items() if ct > 0}
    return alpha_norm, outcome, load


def build_contig_pointers(seq: dict[str, str]
                          ) -> tuple[dict[str, tuple[int, int]], int]:
    '''
    Flatten contigs into a coordinate space.
    Returns (pointer[contig] = (start, end), L).
    '''
    pointer: dict[str, tuple[int, int]] = {}
    offset = 0
    for contig, s in seq.items():
        pointer[contig] = (offset, offset + len(s))
        offset += len(s)
    return pointer, offset


def pick_contig(pointer: dict[str, tuple[int, int]], rando: int
                ) -> tuple[str, int]:
    for contig, (start, end) in pointer.items():
        if start <= rando < end:
            return contig, end - start + 1
    raise RuntimeError(f'Coordinate {rando} fell outside every contig range.')


def draw_replicate(
    seq: dict[str, str],
    pointer: dict[str, tuple[int, int]],
    total_length: int,
    alpha: Alpha,
    outcome: Outcome,
    load_count: int,
    rng: random.Random,
    ) -> list[tuple[str, int, str, str]]:

    '''
    Simulate `load_count` mutations for a single strain.

    Rejection-sampling framework:
        - uniform genomic position;
        - accept with probability alpha[ancestral_base];
        - draw the substituted base from outcome[N1];
        - reject positions already chosen in this replicate.
    '''

    chosen_positions: set[tuple[str, int]] = set()
    picks: list[tuple[str, int, str, str]] = []

    while len(picks) < load_count:
        rando = rng.randrange(total_length)             # Uniform position in genome
        contig, longueur = pick_contig(pointer, rando)  # Map to contig and local position
        local = rng.randrange(longueur - 1)
        pos = local + 1
        n1 = seq[contig][local]                         # Extract ancestral base

        if n1 not in alpha:
            continue
        if rng.random() >= alpha[n1]:                   # Accept/reject based on alpha
            continue
        key = (contig, pos)
        if key in chosen_positions:
            continue

        mutant = rng.random()
        combined = 0.0
        chosen: str | None = None
        for n2, p in outcome[n1].items():               # Draw mutant base from outcome distribution
            combined += p
            if mutant <= combined:
                chosen = n2
                break
        if chosen is None:
            chosen = next(iter(outcome[n1]))
        chosen_positions.add(key)
        picks.append((contig, pos, n1, chosen))
    return picks
