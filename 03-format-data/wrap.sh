#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --output=%x_%j.out

cd /hpc/scratch/Elizabeth.Hunter/entropy/db1

echo "START"
date

bash symlinks.sh accessions.txt filenames.txt fastas/ symlinks/

echo "END"
date
