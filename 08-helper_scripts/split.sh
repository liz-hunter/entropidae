#!/usr/bin/env bash

INPUT_DIR="/hpc/scratch/Elizabeth.Hunter/entropy/analysis/plots"

for file in "$INPUT_DIR"/*.csv; do
  base=$(basename "$file")
  csvcut -c SequenceID,taxid,database,entropy,name,db_score,label,cv_group,cv "$file" > out/"$base"
done
