#!/usr/bin/env bash
#SBATCH -J check
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --output=%x_%j.out

echo "START"
date

md5sum -c md5sum.txt > checks.out

echo "END"
date
