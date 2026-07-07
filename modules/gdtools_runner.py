import json
import os
import shutil
import subprocess
import tempfile
import pandas as pd

from .mutation_classifier import classify_from_gdtools_json


def find_gd_files(samples_dir):
    """
    Walk samples_dir looking for sample_name/data/output.gd (breseq default structure).
    Returns {sample_name: path} dict sorted by sample name.
    """
    gd_files = {}
    for entry in os.scandir(samples_dir):
        if not entry.is_dir():
            continue
        gd_path = os.path.join(entry.path, 'data', 'output.gd')
        if os.path.isfile(gd_path):
            gd_files[entry.name] = gd_path
    return dict(sorted(gd_files.items()))


def find_summary_jsons(samples_dir):
    """
    Walk samples_dir looking for sample_name/data/summary.json (breseq default structure).
    Returns {sample_name: path} dict sorted by sample name.
    """
    json_files = {}
    for entry in os.scandir(samples_dir):
        if not entry.is_dir():
            continue
        json_path = os.path.join(entry.path, 'data', 'summary.json')
        if os.path.isfile(json_path):
            json_files[entry.name] = json_path
    return dict(sorted(json_files.items()))


def _build_annotation(entry):
    """
    Reconstruct annotation string from JSON fields.
    """
    snp_type = entry.get('snp_type', '')
    if snp_type in ('synonymous', 'nonsynonymous', 'nonsense'):
        aa_ref = entry.get('aa_ref_seq', '')
        aa_pos = entry.get('aa_position', '')
        aa_new = entry.get('aa_new_seq', '')
        codon_ref = entry.get('codon_ref_seq', '')
        codon_new = entry.get('codon_new_seq', '')
        return f'{aa_ref}{aa_pos}{aa_new} ({codon_ref}→{codon_new})'
    return entry.get('gene_position', '')


def _build_mutation(entry):
    """
    Reconstruct mutation string.
    SNPs need to match the "N1->N2" format for simulate-dnds/simulate-parallel
    """
    mtype = entry.get('type', '')
    if mtype == 'SNP':
        return f"{entry.get('ref_seq', '')}→{entry.get('new_seq', '')}"
    elif mtype == 'DEL':
        return f"Δ{entry.get('size', '')} bp"
    elif mtype == 'INS':
        new_seq = entry.get('new_seq', '')
        return f'+{len(new_seq)} bp' if new_seq else f"+{entry.get('size', '')} bp"
    else:
        return mtype


def run_gdtools_compare(gdtools_path, reference, gd_files, output_path):
    """
    Run gdtools COMPARE -f JSON on the provided gd_files, then flatten the JSON
    entries into the same wide table shape the rest of the pipeline expects:
    one row per mutation, one column per sample (frequency), and a
    mutation_type precomputed from gdtools' classification fields.

    Because gdtools names output columns after input filenames, passing output.gd
    directly would produce a column named 'output' for every sample.
    Each file is copied to a temporary directory as sample_name.gd.
    The temp directory is then automatically removed.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        renamed_gd_paths = []
        for sample_name, gd_path in gd_files.items():
            dest = os.path.join(tmpdir, f'{sample_name}.gd')
            shutil.copy2(gd_path, dest)
            renamed_gd_paths.append(dest)

        json_path = os.path.join(tmpdir, 'compare.json')
        cmd = [
            gdtools_path, 'COMPARE',
            '-f', 'JSON',
            '-o', json_path,
            '-r', reference,
        ] + renamed_gd_paths

        print(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)

        with open(json_path) as f:
            entries = json.load(f)['entries']

        rows = []
        for entry in entries:
            row = {
                'seq_id': entry.get('seq_id', ''),
                'position': entry.get('position', ''),
                'mutation': _build_mutation(entry),
                'annotation': _build_annotation(entry),
                'gene': entry.get('gene_name', ''),
                'description': entry.get('gene_product', ''),
                'mutation_type': classify_from_gdtools_json(
                    entry.get('snp_type', ''), entry.get('gene_position', '')
                ),
            }
            for sample_name in gd_files:
                freq = entry.get(f'frequency_{sample_name}', '0')
                row[sample_name] = '' if freq == '0' else freq
            rows.append(row)

        columns = ['seq_id', 'position', 'mutation', 'annotation', 'gene', 'description', 'mutation_type'] + list(gd_files)
        df = pd.DataFrame(rows, columns=columns)
        df.to_csv(output_path, index=False)

        print(f'Saved gdtools output: {output_path}\n')


def run_gdtools_count(gdtools_path, reference, gd_files, output_path):
    """
    Run gdtools COUNT -b -p on .gd files by category and substitution type
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        renamed_gd_paths = []
        for sample_name, gd_path in gd_files.items():
            dest = os.path.join(tmpdir, f'{sample_name}.gd')
            shutil.copy2(gd_path, dest)
            renamed_gd_paths.append(dest)

        cmd = [
            gdtools_path, 'COUNT',
            '-r', reference,
            '-o', output_path,
            '-b',
            '-p',
        ] + renamed_gd_paths

        print(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)

        print(f'Saved gdtools base-substitution counts: {output_path}\n')
