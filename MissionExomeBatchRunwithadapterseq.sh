#!/bin/bash

set -eu

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <md5sum_filename>"
    exit 1
fi

md5sum_file=$1

echo "Checking for the presence of the md5sum file..."
if [ ! -f "$md5sum_file" ]; then
    echo "File not found: $md5sum_file"
    exit 1
fi

echo "Extracting filenames..."
filenames=$(awk '{ if ($2 ~ /R1\.(fq|fastq)\.gz$/) print $2 }' "$md5sum_file")
echo "Filenames extracted: $filenames"

if [ -z "$filenames" ]; then
    echo "No matching filenames found."
    exit 1
fi

for filename in $filenames; do
  echo "Processing filename: $filename"
  IFS='-_.' read -ra fields <<< "$filename"
  for field in "${fields[@]}"; do
    if [[ $field =~ ^[0-9]{6}(P|P1|P2|P3|P4|P5|P6|M|F)$ ]]; then
      echo "Match found: $field"
      sampleID=$field
      instID=${sampleID:0:2}

      echo "Searching for instdir..."
      #instdir=$(find /home/diag_ashwin/exome/Mission_project/ -maxdepth 1 -type d -name "$instID*" -print -quit)
      instdir=$(find /home/diag_ashwin/exome/Mission_project/ -maxdepth 1 -type d -name "$instID*" ! -name "*Fastqfiles_forAWS" -print -quit)
      echo "instdir found: $instdir"

      echo "Searching for sampledir..."
      sampledir=$(find "$instdir" -maxdepth 2 -type d -name "${sampleID//[[:alpha:]]/}" -print -quit)
      echo "sampledir found: $sampledir"

      echo "Checking for .tsv and .fastq.gz files..."
      if [[ -z $(find $sampledir -type f -name '*.tsv') && -n $(find $sampledir -type f -name "$sampleID*.fastq.gz") ]]; then
         echo "Starting sbatch..."
         (cd $sampledir && sbatch exomescript_withadapterseq-v4.sh &) || { echo "sbatch failed"; exit 1; }
      else
         echo "Required files not found or conditions not met."
      fi
    fi
  done
done

