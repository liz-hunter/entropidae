import os
import pandas as pd
import re
from glob import glob

def process_and_sample(input_dir, output_dir, sample_size=1_000_000):
    os.makedirs(output_dir, exist_ok=True)
    taxid_data = {}

    filename_pattern = re.compile(r'_(\d+)_filtered\.tsv$')

    for filepath in glob(os.path.join(input_dir, '*_filtered.tsv')):
        filename = os.path.basename(filepath)
        match = filename_pattern.search(filename)
        if not match:
            print(f"Skipping file {filename}: couldn't extract taxid.")
            continue

        taxid = match.group(1)

        try:
            df = pd.read_csv(filepath, sep='\t', low_memory=False)
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue

        df = df[[col for col in df.columns if not col.startswith('virid')]]

        if taxid not in taxid_data:
            taxid_data[taxid] = [df]
        else:
            taxid_data[taxid].append(df)

    for taxid, dfs in taxid_data.items():
        combined = pd.concat(dfs, ignore_index=True)

        if len(combined) > sample_size:
            combined = combined.sample(n=sample_size, random_state=42)
        else:
            print(f"Warning: Only {len(combined)} rows for taxid {taxid}, saving all.")

        output_file = os.path.join(output_dir, f"{taxid}_1M_subsample.tsv")
        combined.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python subsample.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    process_and_sample(input_dir, output_dir)
