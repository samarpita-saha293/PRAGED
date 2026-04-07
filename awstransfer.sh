#!/bin/bash

set -euo pipefail

# Function to ask the user whether they want to upload exome data or genome data
ask_data_type() {
  echo "Do you want to upload exome data or genome data? (exome/genome)"
  read data_type
  if [[ "$data_type" == "exome" || "$data_type" == "Exome" || "$data_type" == "EXOME" ]]; then
    cd /home/diag_ashwin/exome/Mission_project || { echo "Failed to navigate to /home/diag_ashwin/exome/Mission_project"; exit 1; }
    echo "Navigated to /home/diag_ashwin/exome/Mission_project"
  elif [[ "$data_type" == "genome" || "$data_type" == "Genome" || "$data_type" == "GENOME" ]]; then
    cd /home/diag_ashwin/genome/Mission_project || { echo "Failed to navigate to /home/diag_ashwin/genome/Mission_project"; exit 1; }
    echo "Navigated to /home/diag_ashwin/genome/Mission_project"
  else
    echo "Invalid option. Please enter 'exome' or 'genome'."
    ask_data_type
  fi
}

ask_data_type

# Initialize variables
declare -a uploadStatuses
declare -a uploadPaths

# Define associative array for mapping local folder names to S3 bucket names
declare -A projectMap=(
  ["01-AIIMS-Delhi"]="01-aiims-delhi"
  ["02-BGSBU-Jammu"]="02-bgsbu-jammu"
  ["03-BHU-Varanasi"]="03-bhu-varanasi"
  ["04-CUP-Punjab"]="04-cup-punjab"
  ["05-CDFD-Hyderabad"]="05-cdfd-hyderabad"
  ["06-CHG-Bengaluru"]="06-chg-bengaluru"
  ["07-NIRTH-Jabalpur"]="07-nirth-jabalpur"
  ["08-ILS-Bhubaneshwar"]="08-ils-bhubaneshwar"
  ["09-NIMHANS-Bengaluru"]="09-nimhans-bengaluru"
  ["10-NIRRH-Mumbai"]="10-nirrh-mumbai"
  ["11-NIMS-Hyderabad"]="11-nims-hyderabad"
  ["12-RGCB-Trivandrum"]="12-rgcb-trivandrum"
  ["13-SGPIMS-Lucknow"]="13-sgpgims-lucknow"
  ["14-UOJ-Jammu"]="14-uoj-jammu"
  ["15-UOM-Chennai"]="15-uom-chennai"
  ["16-Others"]="16-others"
)

# Function to upload files to S3 and track status
upload_to_s3() {
  local folder=$1  # Folder name e.g., 010065 or 160041/160041P-Blood
  local baseFolder=${folder%%/*}
  local projectID=${baseFolder:0:2}  # Project ID e.g., 01 or 16
  local projectName=""
  local s3BucketName=""

  # Find the project name and corresponding S3 bucket name
  for key in "${!projectMap[@]}"; do
    if [[ "$key" == *"$projectID"* ]]; then
      projectName="$key"
      s3BucketName="${projectMap[$key]}"
      break
    fi
  done

  if [ -z "$projectName" ] || [ -z "$s3BucketName" ]; then
    echo "Error: Could not map the folder to an S3 bucket."
    return
  fi

  # Adjust base directory depending on the data type
  local baseDir=""
  if [[ "$data_type" =~ ^([Ee]xome)$ ]]; then
    baseDir="/home/diag_ashwin/exome/Mission_project"
  else
    baseDir="/home/diag_ashwin/genome/Mission_project"
  fi

  local localPath="$baseDir/$projectName/$folder"

  # 👇 For genome, create uploads under s3://bucket/genome/<folder>
  local s3Path
  if [[ "$data_type" =~ ^([Gg]enome)$ ]]; then
    s3Path="s3://$s3BucketName/genome/$folder"
  else
    s3Path="s3://$s3BucketName/$folder"
  fi

  if [ ! -d "$localPath" ]; then
    echo "The user-provided path $localPath does not exist."
    uploadStatuses+=("Failed")
    uploadPaths+=("$localPath (Path not found)")
    return
  fi

  echo "Processing folder: $localPath"

  # -----------------------
  # EXOME: upload full folder
  # -----------------------
  if [[ "$data_type" =~ ^([Ee]xome)$ ]]; then
    echo "Attempting to upload full folder from $localPath to $s3Path..."
    if aws s3 cp "$localPath" "$s3Path" --profile praged --recursive; then
      echo "Upload successful: $localPath to $s3Path"
      uploadStatuses+=("Success")
    else
      echo "Upload failed: $localPath to $s3Path"
      uploadStatuses+=("Failed")
    fi
    uploadPaths+=("$localPath -> $s3Path")
    return
  fi

  # -----------------------
  # GENOME: selective upload
  # -----------------------
  echo "Uploading selected genome files from $localPath to $s3Path..."

  local filesToUpload=(
    "*recal.bai"
    "*_hardfiltered.vcf.gz.tbi"
    "*_Filtered.tsv"
    "*_hardfiltered.vcf.gz"
    "*_R1.fastq.gz"
    "*_R2.fastq.gz"
    "*recal.bam"
  )

  local successCount=0
  local failCount=0

  # Upload selected files
  for pattern in "${filesToUpload[@]}"; do
    matches=($(find "$localPath" -maxdepth 1 -type f -name "$pattern"))
    if [ ${#matches[@]} -eq 0 ]; then
      echo "No files matching $pattern found in $localPath"
      continue
    fi
    for file in "${matches[@]}"; do
      if aws s3 cp "$file" "$s3Path/" --profile praged; then
        echo "Uploaded: $file"
        ((successCount++))
      else
        echo "Failed to upload: $file"
        ((failCount++))
      fi
    done
  done

  # Upload QC folder (if exists)
  if [ -d "$localPath/QC" ]; then
    echo "Uploading QC folder..."
    if aws s3 cp "$localPath/QC" "$s3Path/QC" --recursive --profile praged; then
      echo "Uploaded QC folder successfully."
      ((successCount++))
    else
      echo "Failed to upload QC folder."
      ((failCount++))
    fi
  else
    echo "No QC folder found."
  fi

  # Track status
  if [ $failCount -eq 0 ]; then
    uploadStatuses+=("Success")
  else
    uploadStatuses+=("Partial/Failed")
  fi
  uploadPaths+=("$localPath -> $s3Path")
}

# Read the folder names from user
echo "Enter the folder names separated by spaces (e.g., 010065 160041/160041P-Blood):"
read -a folderNames

# Call the function for each folder name
for folderName in "${folderNames[@]}"; do
    upload_to_s3 "$folderName"
done

# Summarize upload results
completedUploads=0
for i in "${!uploadStatuses[@]}"; do
    if [ "${uploadStatuses[$i]}" == "Success" ]; then
        ((completedUploads++))
    fi
    echo "Upload ${i}: ${uploadPaths[$i]} - ${uploadStatuses[$i]}"
done

if [ $completedUploads -eq ${#folderNames[@]} ]; then
    echo -e "\n"
    echo " ######   ##      ##    ######    ######   ######   ######    ######"
    echo "##    ##  ##      ##  ##        ##        ##       ##    ##  ##    ##"
    echo "##        ##      ##  ##        ##        ##       ##        ##      "
    echo "  ####    ##      ##  ##        ##        ####       ####      ####  "
    echo "      ##  ##      ##  ##        ##        ##             ##        ##"
    echo "##    ##  ##      ##  ##        ##        ##       ##    ##   ##    ##"
    echo " ######    #######      ######    ######   ######    ######    ###### "
    echo -e "\n"
else
    echo "Some uploads did not complete successfully. Please check the logs above."
fi
