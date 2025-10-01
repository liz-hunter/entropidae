#!/usr/bin/env bash
#SBATCH -J entropy
#SBATCH -c 5
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out

module load micromamba

echo "START"
date

cd $MYPATH/entropy

export MAMBA_ROOT_PREFIX=$HOME/micromamba
eval "$($HOME/bin/micromamba shell hook --shell=bash)"

micromamba activate virid_env

python entropy.py

echo "END"
date
