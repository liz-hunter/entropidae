import csv

# File paths
file1_path = "missing_accessions.txt"
file2_path = "../accession2taxid_masters/missing_accession_key.txt"
output_path = "merged_missing_accessions.txt"

# Read File2 into a dictionary
file2_data = {}
with open(file2_path, newline='') as f2:
    reader = csv.reader(f2, delimiter='\t')
    for row in reader:
        key = row[0]
        file2_data[key] = row[1:]  # skip the first column, which is the key

# Find the number of columns in File2 (excluding the key)
num_file2_cols = max((len(v) for v in file2_data.values()), default=0)

# Merge with File1
with open(file1_path, newline='') as f1, open(output_path, 'w', newline='') as out:
    reader1 = csv.reader(f1, delimiter='\t')
    writer = csv.writer(out, delimiter='\t')
    
    for row1 in reader1:
        key = row1[0]
        row2 = file2_data.get(key, ["NA"] * num_file2_cols)
        writer.writerow(row1 + row2)
