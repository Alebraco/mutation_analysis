import os
import shutil
import subprocess
import tempfile
import pandas as pd


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


def run_gdtools_compare(gdtools_path, reference, gd_files, output_path):
    """
    Run gdtools COMPARE -f TABLE on the provided gd_files.

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

        cmd = [
            gdtools_path, 'COMPARE',
            '-f', 'TABLE',
            '-o', output_path,
            '-r', reference,
        ] + renamed_gd_paths

        print(f'Running: {" ".join(cmd)}')
        subprocess.run(cmd, check=True)

        # Normalize absent-mutation marker.
        df = pd.read_csv(output_path, dtype=str, keep_default_na=False)
        sample_cols = [c for c in gd_files if c in df.columns]
        for col in sample_cols:
            df.loc[df[col] == '0', col] = ''
        df = df.drop(columns='type', errors='ignore')
        df.to_csv(output_path, index=False)

        print(f'Saved gdtools output: {output_path}')
