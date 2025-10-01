#!/usr/bin/env bash
#SBATCH -J inspect
#SBATCH -n 5
#SBATCH --time=30:00:00
#SBATCH --mem=500G

module load kraken2/2.1.5

echo "START"
date

cd $MYPATH/entropy/
kraken2-inspect --db $MYPATH/entropy/virid_db10 > $MYPATH/entropy/kraken2_analysis/inspect/virid_db10_inspect.txt

echo "END"
date
