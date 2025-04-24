#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --output=%x_%j.out

echo "START"
date

cd /hpc/scratch/Elizabeth.Hunter/entropy/
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C kraken_tax/taxonomy

echo "END"
date
