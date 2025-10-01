import pandas as pd
import sys
import os
from pathlib import Path

# Get SLURM array task ID
task_id = int(sys.argv[1])

# Load list of files
with open("filtered_file_list.txt") as f:
    file_list = [line.strip() for line in f.readlines()]

# Pick the correct file
input_file = file_list[task_id]

# Extract filename and taxid from input
basename = os.path.basename(input_file).replace("_filtered.tsv", "")
taxid = int(basename.split("_")[-1])  # Assumes file ends in _<taxid>_filtered.tsv

# Load taxid-to-name key
name_key = pd.read_csv("/hpc/scratch/Elizabeth.Hunter/entropy/keys/taxidofi_to_name.tsv", sep="\t")
tax_name = name_key.loc[name_key["taxid"] == taxid, "name"].values[0]

# Read filtered data
df = pd.read_csv(input_file, sep="\t")

# Drop everything but dist columns and taxid
cols_to_keep = [col for col in df.columns if col.startswith("dist")] + ["true_taxid"]
df = df[cols_to_keep]

# Melt to long
df_long = df.melt(id_vars=["true_taxid"], var_name="database", value_name="entropy")
df_long["database"] = df_long["database"].str.replace("dist", "virid")
df_long["entropy_bin"] = df_long["entropy"].round().astype("Int64")

# Count
counts = df_long.groupby(["database", "entropy_bin"]).size().reset_index(name="count")
counts["taxid"] = taxid
counts["name"] = tax_name

# Write output
Path("analysis/entropy_counts").mkdir(parents=True, exist_ok=True)
out_path = f"/hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy_counts/{basename}_entropy_counts.tsv"
counts.to_csv(out_path, sep="\t", index=False)
print(f"✅ Processed {basename} → {out_path} ({len(counts)} rows)")
