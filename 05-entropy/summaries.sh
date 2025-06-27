#!/usr/bin/env bash
#SBATCH -J summaries
#SBATCH -c 2
#SBATCH --output=%x_%j.out

module load micromamba

echo "START"
date

cd $MYPATH/entropy

export MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$($HOME/bin/micromamba shell hook --shell=bash)"

micromamba activate virid_env

python summaries.py

echo "END"
date
