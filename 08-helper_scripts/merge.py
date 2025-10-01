import pandas as pd
import glob
import os

# Directory containing summary files
data_dir = "./summaries1"  # Change to your actual directory
output_dir = "./summaries"
os.makedirs(output_dir, exist_ok=True)

# Get all *_virid_summary.txt files
virid1_6_files = sorted(glob.glob(os.path.join(data_dir, '*_virid_summary.txt')))
log_file = os.path.join(output_dir, 'merge_log.txt')

with open(log_file, 'w') as log:
    for file1 in virid1_6_files:
        prefix = file1.replace('_virid_summary.txt', '')
        file2 = f'{prefix}_virid6-10_summary.txt'

        if not os.path.exists(file2):
            log.write(f"[MISSING] {file2} not found for {file1}\n")
            continue

        try:
            df1 = pd.read_csv(file1, sep='\t')
            df2 = pd.read_csv(file2, sep='\t')
        except Exception as e:
            log.write(f"[ERROR] Could not read files: {file1} or {file2} — {e}\n")
            continue

        # Merge on SequenceID
        merged = pd.merge(df1, df2, on='SequenceID', suffixes=('_1', '_2'))

        # Check for mismatches in true_taxid
        taxid_mismatch = (merged['true_taxid_1'] != merged['true_taxid_2']).any()
        virid6_mismatch = (merged['virid6_1'] != merged['virid6_2']).any()

        if taxid_mismatch or virid6_mismatch:
            log.write(f"[MISMATCH] {prefix} — true_taxid: {taxid_mismatch}, virid6: {virid6_mismatch}\n")
            continue

        # Drop redundant columns and rename
        merged = merged.drop(columns=['true_taxid_2', 'virid6_2'])
        merged = merged.rename(columns={'true_taxid_1': 'true_taxid', 'virid6_1': 'virid6'})

        # Order columns
        ordered_cols = ['SequenceID', 'true_taxid'] + [f'virid{i}' for i in range(1, 11)]
        merged = merged[ordered_cols]

        # Output filename
        output_file = os.path.join(output_dir, f'{os.path.basename(prefix)}_virid1-10_summary.txt')
        merged.to_csv(output_file, sep='\t', index=False)
        log.write(f"[MERGED] {output_file}\n")
