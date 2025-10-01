import pandas as pd
import re
from pathlib import Path
import argparse

# === 🔧 CONFIGURE THESE DIRECTORIES ===
KEY_DIR = Path("/hpc/scratch/Elizabeth.Hunter/entropy/keys/final_keys")
OUT_DIR = Path("/hpc/scratch/Elizabeth.Hunter/entropy/kraken2_analysis/out")
SUMMARY_DIR = Path("/hpc/scratch/Elizabeth.Hunter/entropy/kraken2_analysis/summaries2")

# === 🚿 Clean sequence ID
def clean_sequence_id(seq):
    seq = seq.lstrip('@')
    seq = re.sub(r'/\d+$', '', seq)
    return seq

# === 📖 Parse .out file and extract SequenceID + numeric taxid
def parse_out_file(file_path):
    cols = ["CoU", "SequenceID", "taxid_raw", "read_length", "lca2kmer"]
    df = pd.read_csv(file_path, sep="\t", header=None, names=cols, dtype=str)
    df[['taxa', 'taxid']] = df['taxid_raw'].str.extract(r'^(.*?) \(taxid (\d+)\)$')
    return df[['SequenceID', 'taxid']]

# === 🔍 Merge a single dataset with validation
def merge_dataset(base_filename, start, end, prefix):
    print(f"\n[▶] Processing: {base_filename} ({prefix}{start} to {prefix}{end})")

    key_file = KEY_DIR / f"{base_filename}_key_final.txt"

    # === Check key file
    if not key_file.exists():
        print(f"[✗] Skipping {base_filename} — missing key_final.txt")
        return

    try:
        key_df = pd.read_csv(key_file, sep="\t", header=None, names=["SequenceID", "assembly", "true_taxid"], dtype=str)
        key_df["SequenceID"] = key_df["SequenceID"].apply(clean_sequence_id)
        result_df = key_df[["SequenceID", "true_taxid"]].copy()
        key_ids = set(result_df["SequenceID"])

        # === Process specified prefixN range
        for i in range(start, end + 1):
            out_file = OUT_DIR / f"{base_filename}_{prefix}{i}.out"
            if not out_file.exists():
                print(f"[✗] Skipping {base_filename} — missing {prefix}{i} outfile")
                return

            df = parse_out_file(out_file)

            out_ids = set(df["SequenceID"])
            missing_in_out = key_ids - out_ids
            extra_in_out = out_ids - key_ids

            if missing_in_out:
                print(f"[✗] {prefix}{i}.out missing {len(missing_in_out)} SequenceIDs from key file")
            if extra_in_out:
                print(f"[✗] {prefix}{i}.out has {len(extra_in_out)} unexpected SequenceIDs not in key file")

            if missing_in_out or extra_in_out:
                print(f"[!] Skipping {base_filename} due to row mismatch in {prefix}{i}.out")
                return

            df.rename(columns={"taxid": f"{prefix}{i}"}, inplace=True)
            result_df = result_df.merge(df, on="SequenceID", how="left")
            print(f"[✓] Merged {prefix}{i}.out — all SequenceIDs accounted for")

        # === Write result
        SUMMARY_DIR.mkdir(parents=True, exist_ok=True)
        output_file = SUMMARY_DIR / f"{base_filename}_{prefix}{start}-{end}_summary.txt"
        result_df.to_csv(output_file, sep="\t", index=False)
        print(f"[💾] Saved: {output_file}")

    except Exception as e:
        print(f"[‼] Error while processing {base_filename}: {e}")

# === 🧭 Detect datasets from OUT_DIR
def auto_detect_base_filenames_from_out_dir(start, prefix):
    pattern = f"*_{prefix}{start}.out"
    out_files = sorted(OUT_DIR.glob(pattern))
    return [f.stem.replace(f"_{prefix}{start}", "") for f in out_files]

# === 🚀 MAIN
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge *_prefixN.out files into summary by range.")
    parser.add_argument("--start", type=int, required=True, help="Start index of prefixN files (inclusive)")
    parser.add_argument("--end", type=int, required=True, help="End index of prefixN files (inclusive)")
    parser.add_argument("--prefix", type=str, default="virid", help="Prefix used in file names (default: virid)")
    args = parser.parse_args()

    base_files = auto_detect_base_filenames_from_out_dir(args.start, args.prefix)
    print(f"[ℹ] Found {len(base_files)} datasets with _{args.prefix}{args.start}.out")

    for base in base_files:
        merge_dataset(base, args.start, args.end, args.prefix)
