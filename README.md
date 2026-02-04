# Genomic Compression and Entropy Workflow

----------

**Objective:** collect publicly available genomes from NCBI and create bootstrap samples (with replacement) to determine how the level of representation of different taxa impacts database identifications and performance

**Expectation:** taxa that are more heavily represented will have more consistent identifications in the bootstrapped samples 

----------

## Part A: Sampling (external)

### #1 Retrieve and Clean Up Accessions

#### Grab relevant accessions from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1)

Retrieve a table of accessions from taxonomic groups of interest (ex. Viridiplantae, Enterobacteriaceae, etc.) from the NCBI genome browser.

The following columns are required:
  - Assembly Name
  - Assembly Accession
  - Assembly Paired Assembly Accession
  - Organism Name
  - Organism Taxonomic ID

#### Run clean_up.R wrapper to remove RefSeq/Genbank duplicates, and organellar genomes

```
Rscript clean_up.R <input_list.tsv> [output_cleaned_list.tsv] \
                    [--remove-organellar] [--quiet-organellar] [--no-report]
```
*Arguments:*
input_list.tsv  
Input genome metadata table (TSV from NCBI genome browser)

output_cleaned_list.tsv (optional)
Output cleaned table (default: <input_basename>_cleaned.tsv, written to the same directory as the input)

*Options:*
--remove-organellar
Remove rows flagged as organellar (mitochondria, plastid, chloroplast, apicoplast)
--quiet-organellar
Do not print the organellar warning/list (still detects and removes if enabled)
--no-report
Suppress cleanup summary and removed-row reports

*Example:*
```
Rscript clean_up.R mydataset_accessions.tsv --remove-organellar
```

#### Run bootstraps.R wrapper to generate bootstrap genome databases

```
Rscript bootstraps.R --in <cleaned.tsv> --n_boot <int> --label <string> \
[--seed <int>] [--outdir DIR] \
[--outputs composite,files,unique_files,accession_counts,taxid_counts] \
[--quiet]
```

*Arguments:*
--in  
Cleaned TSV with columns: accession, name, taxid

--n_boot  
Number of bootstrap samples (with replacement)

--label  
String used for naming outputs (e.g., plants)

*Options:*
--seed  
Integer seed for reproducible sampling

--outdir  
Output directory (default: {label}_boot)

--outputs  
Comma-separated outputs to produce (reduces from default).  
Allowed: composite, files, unique_files, accession_counts, taxid_counts  
Default (if omitted): ALL outputs

--quiet  
Suppress progress messages (still errors on invalid input)

*Example:*
```
Rscript bootstraps.R --in virid_cleaned.tsv --n_boot 10 --label mydataset --seed 1234
```

#### Optional: Run diversity_metrics.R wrapper to compute taxonomic diversity metrics on the bootstraps

```
Rscript diversity_metrics.R --label <string> [--bootdir DIR] [--taxids FILE] [--out FILE]
[--plots] [--plots_outdir DIR] [--plot_metrics m1,m2,...] [--quiet]
```

*Arguments:*
--label  
String used during bootstrapping (must match bootstrap naming)

*Options:*
--bootdir  
Bootstrap directory (default: {label}_boot)

--taxids  
Composite taxid matrix (default: {bootdir}/composite_{label}_taxids.tsv)

--out  
Output diversity summary table (default: {label}_metrics/diversity_{label}_summary.tsv)

--plots  
Generate diversity metric plots (one bar per database, faceted by metric)

--plots_outdir  
Directory for plots (default: {label}_metrics)

--plot_metrics  
Comma-separated metrics to include in plots  
Default: n_unique_taxa,shannon,gini_simpson,inverse_simpson,pielou_evenness

--quiet  
Suppress progress messages

*Example:*
```
Rscript diversity_metrics.R --label mydataset
```


### #2 Retrieve data from NCBI using ncbi-datasets

mydataset_boot/mydataset_db[1-n]_unique.txt contains the unique taxids contained in each bootstrapped sample. 

You can retrieve the data in a linear way, or optionally chunk the accessions into blocks of n rows, to run in a parallel or serial manner.

*Split accession list into chunks:*
```
awk -v n=1000 '{print > sprintf("mydataset_%d.txt", int((NR-1)/n)+1)}' mydataset_db[1-n]_unique.txt
```

*Recommendations for successful retrieval of large datasets:*
- Use an NCBI API key
- Stagger calls by using random sleep intervals 
- Run a maximum of two download jobs in parallel at a time, more will fail with a "too many calls" error from NCBI
- Ensure you have enough memory, depending on your dataset this can be a huge download 

*Basic ncbi-datasets script:*
```
datasets download genome accession --dehydrated --inputfile mydataset.txt --filename mydataset.zip
unzip mydataset.zip -d mydataset
datasets rehydrate --directory mydataset
md5sum -c mydataset/md5sum.txt > checks.out
```

*Example HPC array with checks:*
```
#!/usr/bin/env bash
#SBATCH -J datasets
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --array=1-4%2
#SBATCH --output=%x_%a_%A.out

echo "START"
date

module load ncbi-datasets-cli/16.27 # sub with your own ncbi-datasets install

cd $MYPATH/mydataset_boot/

sleep $(( (RANDOM % 5 + 1) * 60 ))

export NCBI_API_KEY=myapikey123 # sub with your NCBI_API_KEY (or remove/configure globally)

FILE=$SLURM_ARRAY_TASK_ID

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

### #3 Generate symlinks for library build

Use symlinks.sh to count occurrences and generate symlinks 

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
Use with symlink_wrap.sh (this takes while, ~4hrs for virid):
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
ln -s full/path/fastas/GCA_002916435.2_ASM291643v2_genomic.fna full/path/symlinks/GCA_002916435.2_dup2.fna
```

---------

## Part B: Taxonomy

### #1 Download taxonomy files

Grab the master taxonomy files (can reuse this for subsequent builds).
```
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
```

Unpack the taxdump.

```
tar -xvzf taxdump.tar.gz
```

Concatenate the assembly summaries into all_assembly.txt

Run accession2taxid.py to format

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
