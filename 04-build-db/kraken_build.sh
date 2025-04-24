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

echo "END"
date
