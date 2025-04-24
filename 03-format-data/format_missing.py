import csv

# Input and output files
input_file = "missing_accessions.txt"       # TSV: accession<TAB>taxid
output_file = "missing_accession2taxid.map" # Output in Kraken2 format

# Read the missing accessions and their taxids
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    reader = csv.reader(infile, delimiter="\t")
    for row in reader:
        if len(row) < 2:
            continue  # Skip lines that don't have both accession and taxid
        full_acc = row[0].strip()
        taxid = row[1].strip()
        base_acc = full_acc.split('.')[0]  # Remove version number
        outfile.write(f"{base_acc}\t{full_acc}\t{taxid}\t0\n")

print(f"Formatted {output_file} with missing accessions for Kraken2.")
