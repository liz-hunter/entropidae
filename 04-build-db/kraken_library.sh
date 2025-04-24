#!/usr/bin/env bash
#SBATCH -J library
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=48:00:00

module load kraken2/2.1.2

echo "START"
date

cd $MYPATH/entropy/

find $MYPATH/entropy/db1_files/fastas -name '*.fna' -print0 | xargs -P 20 -0 -I{} -n1 kraken2-build --add-to-library {} --db $MYPATH/entropy/virid_db1

echo "END"
date
