import csv

# Input and output files
input_file = "missing_accessions.txt"       # TSV: accession<TAB>taxid
output_file = "missing_accession2taxid.map" # Output in Kraken2 format

success_count = 0
skipped_count = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    reader = csv.reader(infile, delimiter="\t")
    
    for lineno, row in enumerate(reader, start=1):
        if len(row) < 2:
            print(f"[SKIPPED] Line {lineno}: not enough fields → {row}")
            skipped_count += 1
            continue

        full_acc = row[0].strip()
        taxid = row[1].strip()

        if not full_acc or not taxid:
            print(f"[SKIPPED] Line {lineno}: empty field(s) → {row}")
            skipped_count += 1
            continue

        base_acc = full_acc.split('.')[0]
        outfile.write(f"{base_acc}\t{full_acc}\t{taxid}\t0\n")
        success_count += 1

print(f"\n Formatted {success_count} entries to {output_file}")
print(f"\n Skipped {skipped_count} malformed line(s)")
