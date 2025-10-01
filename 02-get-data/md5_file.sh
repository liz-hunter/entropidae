#!/bin/bash

# Usage: ./check_md5_presence.sh file_list.txt

MISSING_LIST="missing_md5_files.txt"
> "$MISSING_LIST"  # empty or create the file

while read FILE; do
    if [[ -f "$FILE/md5sum.txt" ]]; then
        echo -e "${FILE}\tPRESENT"
    else
        echo -e "${FILE}\tMISSING"
        echo "$FILE" >> "$MISSING_LIST"
    fi
done < "$1"
