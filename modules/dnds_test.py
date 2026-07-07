import pandas as pd
from scipy.stats import fisher_exact


def _build_row(sample, obs_syn, obs_nonsyn, obs_nonsense, pos_syn, pos_nonsyn):
    result = {
        'sample': sample,
        'observed_synonymous': obs_syn,
        'observed_nonsynonymous': obs_nonsyn,
        'observed_nonsense': obs_nonsense,
        'observed_dNdS': obs_nonsyn / obs_syn if obs_syn > 0 else float('nan'),
        'expected_dNdS': pos_nonsyn / pos_syn if pos_syn > 0 else float('nan'),
    }
    if obs_syn + obs_nonsyn > 0:
        result['p_value_elevated'] = fisher_exact(
            [[obs_nonsyn, obs_syn], [pos_nonsyn, pos_syn]], alternative='greater')[1]
    else:
        result['p_value_elevated'] = float('nan')
    return result


def compute_dnds_test(count_df):
    '''
    For each sample and pooled across all samples, test if the observed
    genome-wide dN/dS ratio is elevated compared to neutral expectation. 
    Neutral expectation: genome-wide ratio of nonsynonymous to synonymous
    possible sites (from gdtools COUNT --base-substitution-statistics)

    `p_value_elevated` is a one-sided Fisher's exact test comparing each sample's
    observed nonsynonymous and synonymous counts against the possible-site ratio.
    '''
    rows = []
    pooled_syn = pooled_nonsyn = pooled_nonsense = 0
    pos_syn = pos_nonsyn = None

    for _, row in count_df.iterrows():
        obs_syn = int(row['OBSERVED.SYNONYMOUS.TOTAL'])
        obs_nonsyn = int(row['OBSERVED.NONSYNONYMOUS.TOTAL'])
        obs_nonsense = int(row['OBSERVED.NONSENSE.TOTAL'])

        pos_syn = int(row['POSSIBLE.SYNONYMOUS.TOTAL'])
        pos_nonsyn = int(row['POSSIBLE.NONSYNONYMOUS.TOTAL'])

        pooled_syn += obs_syn
        pooled_nonsyn += obs_nonsyn
        pooled_nonsense += obs_nonsense

        rows.append(_build_row(row['sample'], obs_syn, obs_nonsyn, obs_nonsense, pos_syn, pos_nonsyn))

    if pos_syn is not None:
        rows.append(_build_row('POOLED', pooled_syn, pooled_nonsyn, pooled_nonsense, pos_syn, pos_nonsyn))

    return pd.DataFrame(rows)
