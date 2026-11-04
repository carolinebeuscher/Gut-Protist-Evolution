#!/bin/bash
# for macOS
# Usage: bash sum_busco.sh short_summary*.txt

echo -e "inputfile_name\t%_complete\t%_complete_singlecopy\t%_complete_duplicate\t%_fragmented\t%_missing\ttotal_number"

for FILE in "$@"; do
  LINE=$(grep -o "C:[0-9.]*%\[S:[0-9.]*%,D:[0-9.]*%\],F:[0-9.]*%,M:[0-9.]*%,n:[0-9]*" "$FILE")

  if [[ -z "$LINE" ]]; then
    echo "WARNING: Skipping $FILE — BUSCO summary line not found" >&2
    continue
  fi

  C=$(echo "$LINE" | sed 's/C:\([0-9.]*\)%.*/\1/')
  S=$(echo "$LINE" | sed 's/.*S:\([0-9.]*\)%.*/\1/')
  D=$(echo "$LINE" | sed 's/.*D:\([0-9.]*\)%.*/\1/')
  F=$(echo "$LINE" | sed 's/.*F:\([0-9.]*\)%.*/\1/')
  M=$(echo "$LINE" | sed 's/.*M:\([0-9.]*\)%.*/\1/')
  N=$(echo "$LINE" | sed 's/.*n:\([0-9]*\)/\1/')

  echo -e "$(basename "$FILE")\t${C}\t${S}\t${D}\t${F}\t${M}\t${N}"
done
