import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Split columns of a TSV into separate files.")
    parser.add_argument("--input", required=True, help="Path to the input TSV file")
    parser.add_argument("--identifier", required=True, help="Prefix to use for output filenames")

    args = parser.parse_args()

    # Load the TSV file
    df = pd.read_csv(args.input, sep="\t", dtype=str)

    # Ensure no unexpected characters in column names for filenames
    safe_cols = {col: "".join(c if c.isalnum() else "_" for c in col) for col in df.columns}

    # Process each column
    for col in df.columns:
        all_values = df[col].fillna("").astype(str).tolist()
        unique_values = sorted(set(all_values))

        col_safe = safe_cols[col]

        all_file = f"{args.identifier}_{col_safe}.txt"
        unique_file = f"{args.identifier}_{col_safe}_unique.txt"

        with open(all_file, "w") as f:
            f.write("\n".join(all_values) + "\n")

        with open(unique_file, "w") as f:
            f.write("\n".join(unique_values) + "\n")

    print("Files written successfully.")

if __name__ == "__main__":
    main()
