# Genomic Compression and Entropy Workflow

## Part A: Create sample sets for databases
----------
**Objective:** collect publicly available genomes from NCBI and create bootstrap samples (with replacement) to determine how the level of representation of different taxa impacts database identifications and performance

**Expectation:** taxa that are more heavily represented will have more consistent identifications in the bootstrapped samples 

----------

### #1 Grab relevant accessions from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1)

Pull and save lists of accessions from different taxa (Eukaryotes, prokaryotes, viruses, Achaea, and Viridiplantae to start)


### #2 Create bootstrap samples with R (optionally, do some data visualization)

Use samplingwithreplacement.R to output a tsv with the desired number of samples.

Run generate_samples.py to create taxa_sample.txt and taxa_sample_unique.txt on the output from the R script.

```
python gerenate_samples.py --input taxa_samples.tsv --identifier taxa
```

### #3 Download the relevant genomes for each sample 

```
#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 5
#SBATCH --output=%x_%j.out

echo "Starting download"
date

module load ncbi-datasets-cli/16.27

cd /hpc/scratch/Elizabeth.Hunter/entropy/db1
datasets download genome accession --inputfile virid_batch1.txt --filename virid_batch1.zip

echo "END"
date
```

Could definitely use a big array in the future for maximum parallelization - I think I just had the shebang wrong and that's why this wasn't working.


```
#!/usr/bin/env bash

#SBATCH -J bigdata
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --array=1-5
#SBATCH --output=%x_%a_%A.out

echo "START" 
date

module load ncbi-datasets-cli/16.27

cd /hpc/scratch/Elizabeth.Hunter/entropy/db1
FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p virid_batch1.txt)

echo "Current task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Accessing FILE: ${FILE}"

datasets download genome accession --filename ${FILE}.zip ${FILE}

if [ $? -ne 0 ]; then
    echo "${FILE} download failed." >> failed_downloads.txt
else
    echo "${FILE} downloaded successfully."
fi

echo "END"
date
```
Check the downloads using md5 and unzip (can't multithread decompression)

```
md5sum -c md5sum.txt > checks.out

unzip batch1.zip

```

### #4 Count the number of occurrences of each accession in the full sample list, and generate the appropriate number of symlinks

Use symlinks.sh to do this in one shot
```
#!/usr/bin/env bash

accessions_file="$1"
filenames_list="$2"
fna_dir="$3"
output_dir="$4"

# Create output dir if it doesn't exist
mkdir -p "$output_dir"

# Step 1: Count how many times each accession appears
# Output: accession<TAB>count
sort "$accessions_file" | uniq -c | awk '{print $2 "\t" $1}' > accession_counts.tsv

# Step 2: Read all filenames into an array
mapfile -t filenames < "$filenames_list"

# Step 3: Match accessions and create symlinks
while IFS=$'\t' read -r accession count; do
    match=""
    for fname in "${filenames[@]}"; do
        # Extract accession from filename
        acc=$(echo "$fname" | grep -oE '^GC[AF]_[0-9]+\.[0-9]+')
        if [[ "$acc" == "$accession" ]]; then
            match="$fname"
            break
        fi
    done

    if [[ -n "$match" && -f "$fna_dir/$match" ]]; then
        for i in $(seq 1 "$count"); do
            ln -s "$(realpath "$fna_dir/$match")" "$output_dir/${accession}_dup${i}.fna"
        done
    else
        echo "Warning: No file found for accession [$accession] (match: [$match])" >&2
    fi
done < accession_counts.tsv
```
Use with wrapper (this takes while, ~4hrs for virid):
```
bash symlinks.sh accessions.txt filenames.txt fastas/ symlinks/
```

**Bonus:** For if the symlink (or the md5) reports missing files, which can happen if NCBI accessions are suppressed between pulling the accessions and pulling the files. Resample:
```
#grab a random accession from the master virid file to replace it
sample(virid$`Assembly Accession`, size=1)

#find and replace in existing files as needed
's#GCA_047834745\.1#GCA_002916435\.2#g' accessions.txt

#create individual symlink if needed - use absolute paths
ln -s /hpc/scratch/Elizabeth.Hunter/entropy/db1/fastas/GCA_002916435.2_ASM291643v2_genomic.fna /hpc/scratch/Elizabeth.Hunter/entropy/db1/symlinks/GCA_002916435.2_dup2.fna
```

---------

## Part B: Format taxonomy files

### #1 Download taxonomy files

```
#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --output=%x_%j.out

echo "START"
date

cd /hpc/scratch/Elizabeth.Hunter/entropy/
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

echo "END"
date
```

Can't use tar on the computer nodes yet, unpack it on the headnode. 

```
tar -xvzf taxdump.tar.gz
```

Need to create a custom accession2taxid file:
```wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt```
```wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt```

Concatenate these into all_assembly.txt, and then run this python script to format:
```
module load python 3
python accession2taxid.py
```

If there are missing accessions in the output, lookup and add taxids to the file (tsv), and then run format_missing.py
```python format_missing.py```

Concatenate missing_accession2taxid.map with the other file to create a final key.

### #2 Rename fastas using the key

python fix_headers.py $MYPATH/entropy/db1_files/test $MYPATH/entropy/test/taxonomy/accession2taxid.map

 
### #3 Add to library

Using xargs to batch this out. Make sure the CPUS-per-task matches the -P parameter supplied.
```
#!/usr/bin/env bash
#SBATCH -J library
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

module load kraken2/2.1.2

echo "START"
date

cd $MYPATH/entropy/

find $MYPATH/entropy/db1_files/test_sym -name '*.fna' -print0 | xargs -P 20 -0 -I{} -n1 kraken2-build --add-to-library {} --db $MYPATH/entropy/test

echo "END"
date

```
### #3 Build database

taxonomy directory needs to be inside of the db directory

```
#!/usr/bin/env bash
#SBATCH -J build
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time=5-00:00:00
#SBATCH --mem=500G

module load kraken2/2.1.5

echo "START"
date

cd $MYPATH/entropy/
kraken2-build --threads 20 --build --max-db-size 483183820800 --db $MYPATH/entropy/virid_db1

#483183820800 bytes is about 450Gb, need some room for overhead 
#took about 62hrs

echo "END"
date
```
