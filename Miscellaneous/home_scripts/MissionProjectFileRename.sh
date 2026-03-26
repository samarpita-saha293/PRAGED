#!/bin/bash

# Pick either *.fq.gz or *.fastq.gz
shopt -s nullglob
files=(*.fastq.gz *.fq.gz)

if [ ${#files[@]} -eq 0 ]; then
  echo "No FASTQ files found!"
  exit 1
fi

for file in "${files[@]}"; do
  # Strip extension first
  base=$(basename "$file" .fastq.gz)
  base=$(basename "$base" .fq.gz)

  # Split fields by underscores
  IFS="_" read -ra fields <<< "$base"

  sampleID=""
  samplename=""
  pread=""

  for field in "${fields[@]}"; do
    # Keep the sample ID (like 090098F or 050052P)
    if [[ $field =~ ^[0-9]{6}(P[0-9]|P|F|M|R|S|S[0-9]*)$ ]]; then
      sampleID=$field
    # Skip fields like L001, 001, S53, S60
    elif [[ $field =~ ^L[0-9]+$ || $field =~ ^[0-9]+$ || $field =~ ^S[0-9]+$ ]]; then
      continue
    # Read orientation (R1/R2)
    elif [[ $field =~ ^R[12]$ ]]; then
      pread=$field
    # Otherwise, build sample name
    else
      # Capitalize nicely
      field="$(echo "$field" | awk '{print toupper(substr($0,1,1)) tolower(substr($0,2))}')"
      samplename+=$field
    fi
  done

  # New filename
  newname="${sampleID}-${samplename}_${pread}.fastq.gz"
  echo "$file  -->  $newname"

  mv "$file" "$newname"
done
