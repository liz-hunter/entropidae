import os
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Define paths
TAXDUMP_DIR = Path("$MYPATH/entropy/accession2taxid_masters/tax")  # contains nodes.dmp
SUMMARY_DIR = Path("$MYPATH/entropy/kraken2_analysis/summaries")  # *_virid_summary.txt files
OUTPUT_DIR = Path("$MYPATH/entropy/kraken2_analysis/entropy")  # where to write scored outputs

# Load NCBI Taxonomy
def load_taxonomy(nodes_file):
    parent = {}
    rank = {}
    with open(nodes_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t|\t")
            tax_id = parts[0].strip()
            parent_id = parts[1].strip()
            rank_name = parts[2].strip()
            parent[tax_id] = parent_id
            rank[tax_id] = rank_name
    return parent, rank

# Traverse tax tree to root
def get_lineage(taxid, parent_map):
    lineage = []
    while taxid in parent_map and taxid != parent_map[taxid]:
        lineage.append(taxid)
        taxid = parent_map[taxid]
    lineage.append(taxid)  # root
    return lineage

# Calculate distance between two taxids
def tax_distance(t1, t2, parent_map):
    if t1 == t2:
        return 0
    lineage1 = get_lineage(t1, parent_map)
    lineage2 = get_lineage(t2, parent_map)
    set1 = set(lineage1)
    for i, anc2 in enumerate(lineage2):
        if anc2 in set1:
            j = lineage1.index(anc2)
            return i + j  # steps up + steps down
    return None  # No common ancestor (shouldn't happen)

# Analyze a single summary file
def analyze_summary_file(file_path, parent_map):
    df = pd.read_csv(file_path, sep="\t", dtype=str)
    base = file_path.stem.replace("_virid_summary", "")
    tax_columns = [col for col in df.columns if col.startswith("virid")]
    n_reads = len(df)

    # Per-read scores
    per_read_rows = []
    per_sample_summary = []

    for idx, row in df.iterrows():
        sid = row["SequenceID"]
        true_tax = row["true_taxid"]
        dist_vals = []
        unclassified = 0

        for col in tax_columns:
            pred_tax = row[col]
            if pred_tax in ("0", None, "", pd.NA):
                dist_vals.append(None)
                unclassified += 1
            elif true_tax in parent_map and pred_tax in parent_map:
                dist = tax_distance(true_tax, pred_tax, parent_map)
                dist_vals.append(dist)
            else:
                dist_vals.append(None)
                unclassified += 1

        dist_only = [d for d in dist_vals if d is not None]
        mean_dist = sum(dist_only) / len(dist_only) if dist_only else None
        range_dist = max(dist_only) - min(dist_only) if dist_only else None

        row_data = {
            "SequenceID": sid,
            "true_taxid": true_tax,
            **{f"{col}": row[col] for col in tax_columns},
            **{f"dist{col[-1]}": dist_vals[i] for i, col in enumerate(tax_columns)},
            "mean_dist": mean_dist,
            "range_dist": range_dist,
            "unclassified_count": unclassified
        }
        per_read_rows.append(row_data)

    # Save per-read output
    per_read_df = pd.DataFrame(per_read_rows)
    out_path = OUTPUT_DIR / f"{base}_taxdist_per_read.tsv"
    per_read_df.to_csv(out_path, sep="\t", index=False)

    # Per-sample summary
    for i, col in enumerate(tax_columns):
        dist_col = f"dist{col[-1]}"
        col_vals = per_read_df[dist_col].dropna().astype(float)
        row_summary = {
            "sample_virid": f"{base}_{col}",
            "mean": col_vals.mean(),
            "median": col_vals.median(),
            "q25": col_vals.quantile(0.25),
            "q75": col_vals.quantile(0.75),
            "sum": col_vals.sum(),
            "normalized_sum": col_vals.sum() / n_reads if n_reads else None,
            "unclassified_reads": (per_read_df[dist_col].isna()).sum()
        }
        per_sample_summary.append(row_summary)

    return per_sample_summary

# Main script
def run_taxonomic_distance_scoring():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    parent_map, rank_map = load_taxonomy(TAXDUMP_DIR / "nodes.dmp")

    all_summaries = []
    for summary_file in sorted(SUMMARY_DIR.glob("*_virid_summary.txt")):
        print(f"Processing: {summary_file.name}")
        file_summary = analyze_summary_file(summary_file, parent_map)
        all_summaries.extend(file_summary)

    # Save global summary
    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv(OUTPUT_DIR / "taxdist_summary_all_samples.tsv", sep="\t", index=False)

run_taxonomic_distance_scoring()
