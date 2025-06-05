#!/usr/bin/env bash
#SBATCH -J fixhead
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --output=%x_%j.out

module load python/3.12.5 

echo "START"
date

python fix_headers.py $MYPATH/entropy/db2_files/fastas $MYPATH/entropy/db2_files/accession2taxid.map

echo "END"
date
