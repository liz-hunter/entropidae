#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500G
#SBATCH --array=1-66%10
#SBATCH --output=%x_%A_%a.out

echo "START"
date
module load kraken2/2.1.5

cd $MYPATH/entropy/
FILE=$(head -n $SLURM_ARRAY_TASK_ID list.txt | tail -n 1)

kraken2 --paired --threads 20 --use-names --report-minimizer-data \
--gzip-compressed \
--output $MYPATH/entropy/kraken2_reports/${FILE}_virid1.out \
--report $MYPATH/entropy/kraken2_reports/${FILE}_report.txt \
--db $MYPATH/entropy/virid_db1 \
$MYPATH/metagen/mixtures/${FILE}-1_R1.fastq.gz \
$MYPATH/metagen/mixtures/${FILE}-1_R2.fastq.gz

echo "END"
date
