#!/bin/bash

# Prompt the user to choose between exome and genome
read -p "Enter the type of sequencing (exome/genome): " sequencing_type

# Set the base path based on the user's choice
if [[ "$sequencing_type" == "exome" ]]; then
  base_path="/home/diag_ashwin/exome/Mission_project/"
elif [[ "$sequencing_type" == "genome" ]]; then
  base_path="/home/diag_ashwin/genome/Mission_project/"
else
  echo "Invalid sequencing type. Please enter 'exome' or 'genome'."
  exit 1
fi

# Set the source directory
read -p "Enter the source directory path: " source_directory

# Navigate to the source directory
cd "$source_directory" || { echo "Source directory not found"; exit 1; }

for file in *.fastq.gz; do
  if [ -f "$file" ]; then
    filename=$(basename "$file")
    patientid=$(echo "$filename" | cut -d'-' -f1)
    instid=$(echo "$patientid" | cut -c 1,2)
    stripped_patientid=$(echo "$patientid" | cut -c 1-6)

    if [[ $instid =~ ^(0[1-9]|1[0-6])$ ]]; then
      instdir=$(find "$base_path" -maxdepth 1 -type d -name "$instid*" ! -name "*Fastqfiles_forAWS" -print -quit)

      if [ -n "$instdir" ]; then
        patientdir=$(find "$instdir" -maxdepth 1 -type d -name "$stripped_patientid" -print -quit)

        if [ -z "$patientdir" ]; then
          patientdir="$instdir/$stripped_patientid"
          mkdir -p "$patientdir"
          echo "Created patient directory: $patientdir"
        fi

        cp "$file" "$patientdir"
        echo "$file has been transferred to $patientdir"
      else
        echo "Institution directory not found for $instid, skipping transfer."
      fi
    else
      echo "File renaming is not done for $filename"
    fi
  fi
done

