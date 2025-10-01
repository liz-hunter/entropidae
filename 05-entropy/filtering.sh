#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -c 1
#SBATCH --array=0-137
#SBATCH --output=%x_%a_%A.out

module load micromamba

echo "START"
date

cd $MYPATH/entropy

export MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$($HOME/bin/micromamba shell hook --shell=bash)"

micromamba activate virid_env

python filtering.py ${SLURM_ARRAY_TASK_ID}

echo "END"
date
