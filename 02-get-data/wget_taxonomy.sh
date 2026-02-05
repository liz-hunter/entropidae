#!/usr/bin/env bash
#SBATCH -J tax_grab
#SBATCH -c 1
#SBATCH --output=%x_%j.out

echo "START"
date

cd $MYPATH/entropy/
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt

tar -xzf taxdump.tar.gz -C kraken_tax/taxonomy

echo "END"
date
