import pandas as pd
import os
import sys

index = int(sys.argv[1])
files_of_interest = pd.read_csv("/hpc/scratch/Elizabeth.Hunter/entropy/keys/files_of_interest_miseq.tsv", sep="\t")
row = files_of_interest.iloc[index]

filename = row['filename']
taxid = row['taxid']

input_path = f"/hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy/{filename}_taxdist_per_read.tsv"
output_path = f"/hpc/scratch/Elizabeth.Hunter/entropy/analysis/filtered_entropy/{filename}_{taxid}_filtered.tsv"

try:
    df = pd.read_csv(input_path, sep="\t")
    df_filtered = df[df['true_taxid'] == taxid]
    df_filtered.to_csv(output_path, sep="\t", index=False)
    print(f"✅ {filename} done: {len(df_filtered)} rows")
except Exception as e:
    print(f"❌ Failed on {filename}: {e}")
