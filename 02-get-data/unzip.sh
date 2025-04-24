#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --output=%x_%j.out

echo "START"
date

unzip virid_batch2.zip

echo "END"
date
