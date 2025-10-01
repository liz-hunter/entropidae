#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500G
#SBATCH --array=1-11
#SBATCH --output=%x_%A_%a.out

echo "START"
date
module load kraken2/2.1.5

cd $MYPATH/entropy/
FILE=$(head -n $SLURM_ARRAY_TASK_ID sample_lists/real_list.txt | tail -n 1)

kraken2 --paired --threads 20 --use-names --report-minimizer-data \
--gzip-compressed \
--output $MYPATH/entropy/kraken2_analysis/out/${FILE}_virid9.out \
--report $MYPATH/entropy/kraken2_analysis/reports/${FILE}_virid9.txt \
--db $MYPATH/entropy/virid_db9 \
$MYPATH/metagen/euk_meta/trim/${FILE}_F_trimmed_paired.fastq.gz \
$MYPATH/metagen/euk_meta/trim/${FILE}_R_trimmed_paired.fastq.gz

echo "END"
date
