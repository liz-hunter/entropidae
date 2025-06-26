#!/usr/bin/env bash
#SBATCH -J fixhead
#SBATCH -c 1
#SBATCH --output=%x_%j.out

module load python/3.12.5 

echo "START"
date

python headers2taxid.py $MYPATH/entropy/db_files/fastas $MYPATH/entropy/db_files/accession2taxid.map

echo "END"
date
