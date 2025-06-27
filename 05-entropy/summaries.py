import pandas as pd
import re
from pathlib import Path

# Configure directory paths 
KEY_DIR = Path("$MYPATH/entropy/keys/final_keys")           # Directory with *_key_final.txt files
OUT_DIR = Path("$MYPATH/entropy/kraken2_analysis/out")           # Directory with *.out files
SUMMARY_DIR = Path("$MYPATH/entropy/kraken2_analysis/summaries")       # Output directory for merged files

# Clean up sequence IDs
def clean_sequence_id(seq):
    seq = seq.lstrip('@')
    seq = re.sub(r'/\d+$', '', seq)
    return seq

# Parse .out file and extract SequenceID + NCBI taxid
def parse_out_file(file_path):
    cols = ["CoU", "SequenceID", "taxid_raw", "read_length", "lca2kmer"]
    df = pd.read_csv(file_path, sep="\t", header=None, names=cols, dtype=str)
    df[['taxa', 'taxid']] = df['taxid_raw'].str.extract(r'^(.*?) \(taxid (\d+)\)$')
    return df[['SequenceID', 'taxid']]

# Merge each set of experiments into single dataset with validation
def merge_dataset(base_filename):
    print(f"\n[START] Processing: {base_filename}")

    key_file = KEY_DIR / f"{base_filename}_key_final.txt"
    virid_files = [OUT_DIR / f"{base_filename}_virid{i}.out" for i in range(1, 7)]

    # Check for the presence of all expected files
    missing = []
    if not key_file.exists():
        missing.append("key_final.txt")
    for i, vf in enumerate(virid_files, 1):
        if not vf.exists():
            missing.append(f"virid{i}.out")

    if missing:
        print(f"[FAIL] Skipping {base_filename} — missing: {', '.join(missing)}")
        return

    try:
        key_df = pd.read_csv(key_file, sep="\t", header=None, names=["SequenceID", "assembly", "true_taxid"], dtype=str)
        key_df["SequenceID"] = key_df["SequenceID"].apply(clean_sequence_id)
        result_df = key_df[["SequenceID", "true_taxid"]].copy()
        key_ids = set(result_df["SequenceID"])

        for i, virid_file in enumerate(virid_files, 1):
            virid_df = parse_out_file(virid_file)

            # Check row count match
            out_ids = set(virid_df["SequenceID"])
            missing_in_out = key_ids - out_ids
            extra_in_out = out_ids - key_ids

            if missing_in_out:
                print(f"[FAIL] virid{i}.out missing {len(missing_in_out)} SequenceIDs from key file")
            if extra_in_out:
                print(f"[FAIL] virid{i}.out has {len(extra_in_out)} unexpected SequenceIDs not in key file")

            if missing_in_out or extra_in_out:
                print(f"[!] Skipping {base_filename} due to row mismatch in virid{i}.out")
                return

            virid_df.rename(columns={"taxid": f"virid{i}"}, inplace=True)
            result_df = result_df.merge(virid_df, on="SequenceID", how="left")
            print(f"[CORRECT] Merged virid{i}.out — all SequenceIDs accounted for")

        SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
        output_file = SUMMARY_DIR / f"{base_filename}_virid_summary.txt"
        result_df.to_csv(output_file, sep="\t", index=False)
        print(f"[SAVED] Saved: {output_file}")

    except Exception as e:
        print(f"[‼] Error while processing {base_filename}: {e}")

# Detect all samples in OUT_DIR
def auto_detect_base_filenames_from_out_dir():
    out_files = sorted(OUT_DIR.glob("*_virid1.out"))
    return [f.stem.replace("_virid1", "") for f in out_files]

# Main script
if __name__ == "__main__":
    base_files = auto_detect_base_filenames_from_out_dir()
    print(f"[ℹ] Found {len(base_files)} datasets with _virid1.out")
    for base in base_files:
        merge_dataset(base)
