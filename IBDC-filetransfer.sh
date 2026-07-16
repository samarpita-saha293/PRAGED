#!/bin/bash
#SBATCH --job-name="IBDC_filetransfer"
#SBATCH --output="IBDC_414.%j.%N.out"
#SBATCH --error="IBDC_414.%j.%N.err"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=amit.kumar10nov@gmail.com

source_dir="/home/diag_ashwin/exome/Mission_project"
destination_dir="/home/diag_ashwin/exome/Mission_project/IBDC/batch14"

# Create/clear error log
: > error_log.txt

# Make sure destination base folder exists
mkdir -p "$destination_dir"

echo "Using source_dir: $source_dir"
echo "Using destination_dir: $destination_dir"
echo "Reading input from: $(realpath input.txt)"

while IFS= read -r line || [[ -n "$line" ]]; do
    # Clean line
    line=$(echo "$line" | tr -d '\r' | xargs)

    # Skip blank lines
    [[ -z "$line" ]] && continue

    first_two_digits="${line:0:2}"
    echo "Processing line: $line with prefix $first_two_digits"

    # Find matching top-level folders
    mapfile -t matching_folders < <(find "$source_dir" -maxdepth 1 -type d -name "${first_two_digits}-*" 2>/dev/null)

    echo "Matching folders found: ${matching_folders[*]}"

    if [ ${#matching_folders[@]} -gt 0 ]; then
        for folder in "${matching_folders[@]}"; do
            if [ -d "$folder" ]; then
                echo "Found matching folder: $folder"

                # Corrected find command
                mapfile -t fastq_files < <(
                    find "$folder" -maxdepth 3 -type f \( -name "${line}*.fastq.gz" -o -name "Re-${line}*.fastq.gz" \)
                )

                if [ ${#fastq_files[@]} -gt 0 ]; then
                    echo "Found FASTQ files for $line: ${fastq_files[*]}"

                    destination_subdir="$destination_dir/$line"
                    mkdir -p "$destination_subdir" 2>>error_log.txt
                    echo "Creating destination subdir: $destination_subdir"

                    for file in "${fastq_files[@]}"; do
                        filename=$(basename "$file")
                        new_filename=$(echo "$filename" | sed -E 's/^Re-//; s/-[^_]+_/_/')
                        full_destination_path="$destination_subdir/$new_filename"

                        if [ ! -f "$full_destination_path" ]; then
                            if cp "$file" "$full_destination_path"; then
                                echo "File $file copied as $full_destination_path"
                            else
                                echo "Failed to copy file: $file" | tee -a error_log.txt
                            fi
                        else
                            echo "File already exists: $full_destination_path. Skipping..."
                        fi
                    done

                    echo "Generating md5sums for R1 and R2 files"
                    find "$destination_subdir" -type f -name '*_R1*.fastq.gz' -exec md5sum {} + > "$destination_subdir/${line}_R1_md5sum.txt"
                    find "$destination_subdir" -type f -name '*_R2*.fastq.gz' -exec md5sum {} + > "$destination_subdir/${line}_R2_md5sum.txt"

                else
                    echo "No FASTQ files found for $line in $folder."
                fi
            else
                echo "Matching folder not found: $folder"
            fi
        done
    else
        echo "No matching folders found for line: $line"
    fi
done < input.txt

echo "File copy, rename, and md5sum generation is completed."
