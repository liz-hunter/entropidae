#!/usr/bin/env bash
#SBATCH -J kraken2
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=500G
#SBATCH --array=1-126
#SBATCH --output=%x_%A_%a.out

echo "START"
date
module load kraken2/2.1.5

cd $MYPATH/entropy/
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" keys/sample_lists/new_toxin_list.txt)

for i in {1..10}; do
  echo "Running virid_db$i for $FILE"

  kraken2 --paired --threads 20 --use-names --report-minimizer-data \
    --gzip-compressed \
    --output "$MYPATH/entropy/analysis/out/${FILE}_virid$i.out" \
    --report "$MYPATH/entropy/analysis/reports/${FILE}_virid$i.txt" \
    --db "$MYPATH/entropy/virid_db$i" \
    "$MYPATH/metagen/insilicoseq/lod/mixtures/${FILE}-1_R1.fastq.gz" \
    "$MYPATH/metagen/insilicoseq/lod/mixtures/${FILE}-1_R2.fastq.gz"
done

echo "END"
date
