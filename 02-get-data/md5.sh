#!/usr/bin/env bash
#SBATCH -J bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --output=%x_%j.out

echo "START"
date

while read FILE; do
    DIR="$FILE"
    MD5FILE="$DIR/md5sum.txt"

    if [[ ! -f "$MD5FILE" ]]; then
        echo -e "${FILE}\tMISSING_MD5_FILE"
        continue
    fi

    (cd "$DIR" && md5sum -c md5sum.txt > /dev/null 2>&1)
    if [[ $? -eq 0 ]]; then
        echo -e "${FILE}\tPASS"
    else
        echo -e "${FILE}\tFAIL"
    fi
done < "$1"

echo "END"
date
