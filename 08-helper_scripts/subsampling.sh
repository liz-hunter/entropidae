#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=50G
#SBATCH --output=%x_%j.out

module load micromamba

echo "START"
date

cd $MYPATH/entropy/analysis/entropy_filtered

export MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$($HOME/bin/micromamba shell hook --shell=bash)"

micromamba activate virid_env

python subsampling.py /hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy_filtered /hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy_filtered/subsamples

echo "END"
date
