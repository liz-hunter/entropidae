#!/usr/bin/env

#SBATCH -J bigdata
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --array=1-5
#SBATCH --output=%x_%a_%A.out

echo "START" 
date

module load ncbi-datasets-cli/16.27

cd /hpc/scratch/Elizabeth.Hunter/entropy/db1
#FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p virid_batch1.txt)
FILE=$(head -n $SLURM_ARRAY_TASK_ID virid_batch1.txt | tail -n 1)

echo "Current task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Accessing FILE: ${FILE}"

#datasets download genome accession --filename ${FILE}.zip ${FILE}

#if [ $? -ne 0 ]; then
#    echo "${FILE} download failed." >> failed_downloads.txt
#else
#    echo "${FILE} downloaded successfully."
#fi

echo "END"
date
