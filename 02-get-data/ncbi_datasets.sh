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
