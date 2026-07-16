#!/bin/bash

# Define the working directory
working_dir="/home/diag_ashwin/exome/Mission_project/IBDC/batch14/"

# Change to the working directory
cd "$working_dir" || exit 1

# Generate R1.fastq.tsv by finding and listing R1 files
for folder in */; do
    find "$folder" -type d -exec sh -c 'cd "{}" && ls *R1*.fastq.gz 2>/dev/null' \;
done > R1.fastq.tsv

# Generate R2.fastq.tsv by finding and listing R2 files
for folder in */; do
    find "$folder" -type d -exec sh -c 'cd "{}" && ls *R2*.fastq.gz 2>/dev/null' \;
done > R2.fastq.tsv

# Concatenate contents of txt files into md5sum.txt
for folder in */; do
    find "$folder" -type d -exec sh -c '(cd "{}" && cat *.txt 2>/dev/null)' \;
done > md5sum.txt

# Remove intermediate files after use
rm -f R1.fastq.tsv R2.fastq.tsv

echo "Script executed successfully in $working_dir. Intermediate files R1.fastq.tsv and R2.fastq.tsv have been removed."
