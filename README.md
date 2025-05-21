# Genomic Compression and Entropy Workflow

## Part A: Sampling
----------
**Objective:** collect publicly available genomes from NCBI and create bootstrap samples (with replacement) to determine how the level of representation of different taxa impacts database identifications and performance

**Expectation:** taxa that are more heavily represented will have more consistent identifications in the bootstrapped samples 

----------

### #1 Grab relevant accessions from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1)

Pull and save lists of accessions from different taxa (Eukaryotes, prokaryotes, viruses, Achaea, and Viridiplantae to start)


### #2 Create bootstrap samples with R (optionally, do some data visualization)

Use sampling.R to output a tsv with the desired number of bootstrapped samples.

Run generate_samples.py to create taxa_sample.txt and taxa_sample_unique.txt on the tsv output from the R script with the desired number of samples.

```
python gerenate_samples.py --input taxa_samples.tsv --identifier taxa
```

### #3 Download the relevant genomes for each sample 

Split the input (unique) accessions into roughly equal sized batches.

```
awk -v n=1000 '{print > sprintf("chunk_%d.txt", int((NR-1)/n)+1)}' accessions_unique.txt
```

For some reason, even with the API key, submitting more than 2 of these at once causes the NCBI API to fail with too many requests. Running 2 at a time, and just setting it up to run when the first one finishes.


```
#!/usr/bin/env bash
#SBATCH -J datasets
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=1-4
#SBATCH --output=%x_%a_%A.out

echo "START"
date

module load ncbi-datasets-cli/16.27

cd $MYPATH/entropy/db_files

sleep $(( (RANDOM % 5 + 1) * 60 ))

export NCBI_API_KEY=myapikey123

FILE=$SLURM_ARRAY_TASK_ID
#FILE=$((SLURM_ARRAY_TASK_ID + 2))

echo "BEGIN DOWNLOAD"
datasets download genome accession --dehydrated --inputfile chunk_${FILE}.txt --filename batch_${FILE}.zip

echo "BEGIN UNZIP"
unzip batch_${FILE}.zip -d batch${FILE}

echo "REHYDRATE"
datasets rehydrate --directory batch${FILE}

echo "CHECK MD5"
md5sum -c batch${FILE}/md5sum.txt > checks${FILE}.out

echo "END"
date
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
ln -s $MYPATH/entropy/db1_files/fastas/GCA_002916435.2_ASM291643v2_genomic.fna $MYPATH/entropy/db1/symlinks/GCA_002916435.2_dup2.fna
```

---------

## Part B: Taxonomy

### #1 Download taxonomy files

Grab the taxonomy files (can reuse this for subsequent builds).
```
#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --output=%x_%j.out

echo "START"
date

cd $MYPATH/entropy/
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

echo "END"
date
```

Unpack it the taxdump.

```
tar -xvzf taxdump.tar.gz
```

Concatenate the assembly summaries into all_assembly.txt, and then run this python script to format:

```
module load python 3
python accession2taxid.py
```

If there are missing accessions in the output, lookup and add taxids to the file (tsv), and then run format_missing.py
```python format_missing.py```

Concatenate missing_accession2taxid.map with the accession2taxid.map to create the final accession2taxid.map file.

### #2 Rename fastas using the key

python fix_headers.py $MYPATH/entropy/db1_files/test $MYPATH/entropy/test/taxonomy/accession2taxid.map


--------------

## Part C: Database

- Make taxonomy directory inside the db directory
- Move the accession2taxid.map into the db directory

### #1 Add to library

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

### #2 Build the database
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

### #3 Run samples on the newly built database

```
#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500G
#SBATCH --array=1-66%10
#SBATCH --output=%x_%A_%a.out

echo "START"
date
module load kraken2/2.1.5

cd $MYPATH/entropy/
FILE=$(head -n $SLURM_ARRAY_TASK_ID list.txt | tail -n 1)

kraken2 --paired --threads 20 --use-names --report-minimizer-data \
--gzip-compressed \
--output $MYPATH/entropy/kraken2_reports/${FILE}_virid1.out \
--report $MYPATH/entropy/kraken2_reports/${FILE}_report.txt \
--db $MYPATH/entropy/virid_db1 \
$MYPATH/metagen/mixtures/${FILE}-1_R1.fastq.gz \
$MYPATH/metagen/mixtures/${FILE}-1_R2.fastq.gz

echo "END"
date
```