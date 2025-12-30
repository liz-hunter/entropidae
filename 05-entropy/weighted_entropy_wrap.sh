#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --output=%x_%j.out

module load python3
module load micromamba

micromamba activate ~/micromamba/envs/virid_env
micromamba shell reinit --shell bash

echo "START"
date

cd $MYPATH/entropy_onto
python weighted_entropy.py \
  --nodes-dmp $MYPATH/entropy_stan/stan_full/taxonomy/nodes.dmp \
  --input-tsv entropy_test.txt \
  --output-tsv classifications_with_entropy.tsv \
  --alpha-up 0.3 \ #default
  --alpha-down 1.0 #default

echo "END"
date
