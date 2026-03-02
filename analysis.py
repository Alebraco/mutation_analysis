import os
import pandas as pd
from utils import get_strain_columns

def mutation_analysis(df, ancestor, output_dir='.'):
    '''
    Identify parallel mutations (same site OR same gene)
    Indicates how many and which lines share the mutation
    '''

    os.makedirs(output_dir, exist_ok=True)
    strain_cols = get_strain_columns(df, ancestor)

    # Copy relevant columns for parallel mutation analysis
    analysis_df = df[['position', 'mutation', 'annotation', 'gene', 'description']].copy()

    # Create a True/False matrix of mutation presence for each strain
    mutation_presence = df[strain_cols].notna()

    # 1. SITE-LEVEL PARALLEL MUTATIONS (SAME POSITION)

    # Identify rows where mutation is shared at the site level
    same_site_count = mutation_presence.sum(axis=1) 
    same_site_shared = same_site_count > 1

    if same_site_shared.any():
        site_df = analysis_df[same_site_shared].copy()
        site_presence = mutation_presence[same_site_shared]
        site_df['sample_count'] = site_presence.sum(axis=1)

        # Add strain columns indicating presence or absence (1 or 0)
        for strain in strain_cols:
            site_df[strain] = site_presence[strain].astype(int)

        # Add strain name list column
        def get_strain_list(row):
            strains = [strain for strain in strain_cols if row[strain] == 1]
            return ','.join(strains)
        
        site_df['strains'] = site_df.apply(get_strain_list, axis=1)

        # Save site-level parallel mutations to file
        site_file = os.path.join(output_dir, 'site_parallel_mutations.xlsx')
        site_df.to_excel(site_file, index=False)
        print(f'Saved site-level parallel mutations: {site_file}')
        print(f'Total site-level parallel mutations: {len(site_df)}')
    else:
        print('No site-level parallel mutations found.')
        site_df = None

    # 2. GENE-LEVEL PARALLEL MUTATIONS (SAME GENE)
    mutation_presence['gene'] = df['gene']

    gene_strain_matrix = mutation_presence.groupby('gene')[strain_cols].any()
    gene_count = gene_strain_matrix.sum(axis=1)
    gene_shared = gene_count > 1

    if gene_shared.any():
        gene_list = []

        for gene in gene_strain_matrix[gene_shared].index:
            if pd.isna(gene):
                continue
            
            gene_row = df[df['gene'] == gene].iloc[0]
            description = gene_row['description']

            shared_gene_strains = [strain for strain in strain_cols
                                   if gene_strain_matrix.loc[gene, strain]]
            
            gene_info = {
                'gene': gene,
                'description': description,
                'shared_strains': ', '.join(shared_gene_strains),
                'strain_count': len(shared_gene_strains)
            }

            for strain in strain_cols:
                gene_info[strain] = 1 if strain in shared_gene_strains else 0

            gene_list.append(gene_info)
        
        gene_df = pd.DataFrame(gene_list).sort_values(by='strain_count', ascending=False)
        
        gene_file = os.path.join(output_dir, 'gene_parallel_mutations.xlsx')
        gene_df.to_excel(gene_file, index=False)
        print(f'Saved gene-level parallel mutations: {gene_file}')
        print(f'Total gene-level parallel mutations: {len(gene_df)}')
    else:
        print('No gene-level parallel mutations found.')
        gene_df = None

    # 3. ISOLATE MUTATIONS (UNIQUE TO ONE STRAIN)
    isolate_mutations = same_site_count == 1

    if isolate_mutations.any():
        isolate_df = analysis_df[isolate_mutations].copy()
        isolate_presence = mutation_presence[isolate_mutations]

        for strain in strain_cols:
            isolate_df[strain] = isolate_presence[strain].astype(int)
        
        isolate_df.insert(0, 'mutation_number', range(1, len(isolate_df) + 1))

        isolate_file = os.path.join(output_dir, 'unique_mutations.xlsx')
        isolate_df.to_excel(isolate_file, index=False)
        print(f'Saved unique mutations: {isolate_file}')
        print(f'Total unique mutations: {len(isolate_df)}')
    else:
        print('No unique mutations found.')
        isolate_df = None

    return site_df, gene_df, isolate_df