#!/usr/bin/env bash
set -euo pipefail

############################################
# 1. ASK DATA TYPE
############################################
ask_data_type() {
  echo "Do you want to upload exome or genome data? (exome/genome)"
  read -r data_type

  case "$data_type" in
    exome|Exome|EXOME)
      DATA_TYPE="exome"
      SOURCE_BASE="/home/diag_ashwin/exome/Mission_project"
      FRANKLIN_BASE="/home/diag_ashwin/exome/Mission_project/franklin-upload"
      ;;
    genome|Genome|GENOME)
      DATA_TYPE="genome"
      SOURCE_BASE="/home/diag_ashwin/genome/Mission_project"
      FRANKLIN_BASE="/home/diag_ashwin/genome/Mission_project/franklin-upload"
      ;;
    *)
      echo "❌ Invalid option. Enter exome or genome."
      ask_data_type
      ;;
  esac
}
ask_data_type

############################################
# 2. DATE-BASED TARGET DIRECTORY
############################################
TODAY=$(date +"%d-%m-%Y")
TARGET_DIR="${FRANKLIN_BASE}/${TODAY}"

mkdir -p "${TARGET_DIR}"
echo "📁 Target directory created:"
echo "   ${TARGET_DIR}"
echo

############################################
# 3. PROJECT MAP
############################################
declare -A projectMap=(
  ["01"]="01-AIIMS-Delhi"
  ["02"]="02-BGSBU-Jammu"
  ["03"]="03-BHU-Varanasi"
  ["04"]="04-CUP-Punjab"
  ["05"]="05-CDFD-Hyderabad"
  ["06"]="06-CHG-Bengaluru"
  ["07"]="07-NIRTH-Jabalpur"
  ["08"]="08-ILS-Bhubaneshwar"
  ["09"]="09-NIMHANS-Bengaluru"
  ["10"]="10-NIRRH-Mumbai"
  ["11"]="11-NIMS-Hyderabad"
  ["12"]="12-RGCB-Trivandrum"
  ["13"]="13-SGPIMS-Lucknow"
  ["14"]="14-UOJ-Jammu"
  ["15"]="15-UOM-Chennai"
  ["16"]="16-Others"
)

############################################
# 4. READ SAMPLE IDS
############################################
echo "Enter space-separated sample IDs (e.g. 160029 109998):"
read -r -a SAMPLE_IDS
echo

############################################
# 5. COPY FILES (STEP 1)
############################################
echo "================ COPYING FILES ===================="

COPIED_ANY=false

for SAMPLE in "${SAMPLE_IDS[@]}"; do
  PROJECT_ID="${SAMPLE:0:2}"
  PROJECT_FOLDER="${projectMap[$PROJECT_ID]}"

  if [[ -z "${PROJECT_FOLDER}" ]]; then
    echo "⚠️ Unknown project for sample ${SAMPLE}, skipping"
    continue
  fi

  SOURCE_DIR="${SOURCE_BASE}/${PROJECT_FOLDER}/${SAMPLE}"

  if [[ ! -d "${SOURCE_DIR}" ]]; then
    echo "❌ Directory not found: ${SOURCE_DIR}"
    continue
  fi

  echo "🔍 Searching in: ${SOURCE_DIR}"

  shopt -s nullglob
  if [[ "${DATA_TYPE}" == "genome" ]]; then
    FILES=("${SOURCE_DIR}"/*.g.vcf.gz)
  else
    FILES=("${SOURCE_DIR}"/*.vcf)
  fi
  shopt -u nullglob

  if [[ "${#FILES[@]}" -eq 0 ]]; then
    echo "⚠️ No matching files found for ${SAMPLE}"
    continue
  fi

  cp "${FILES[@]}" "${TARGET_DIR}/"
  echo "✅ Copied ${#FILES[@]} file(s) from ${SAMPLE}"
  COPIED_ANY=true
done

if [[ "${COPIED_ANY}" != true ]]; then
  echo "❌ No files were copied. Exiting."
  exit 1
fi

echo "================ COPY STEP COMPLETE ================"
echo

############################################
# 6. GUNZIP (GENOME ONLY) – STEP 2
############################################
if [[ "${DATA_TYPE}" == "genome" ]]; then
  echo "================ DECOMPRESSING gVCFs ==============="
  cd "${TARGET_DIR}"

  shopt -s nullglob
  GZ_FILES=(*.g.vcf.gz)
  shopt -u nullglob

  if [[ "${#GZ_FILES[@]}" -eq 0 ]]; then
    echo "❌ No .g.vcf.gz files found to decompress. Exiting."
    exit 1
  fi

  gunzip -f "${GZ_FILES[@]}"
  echo "✅ Decompressed ${#GZ_FILES[@]} gVCF files"
  echo "================ GUNZIP STEP COMPLETE =============="
  echo
fi

############################################
# 7. AWS CREDENTIALS (INLINE – TEMP)
############################################
AWS_ACCESS_KEY_ID="AKIA2ZHDMMUW4CUWJW7N"
AWS_SECRET_ACCESS_KEY="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"
AWS_DEFAULT_REGION="eu-west-1"

export AWS_ACCESS_KEY_ID
export AWS_SECRET_ACCESS_KEY
export AWS_DEFAULT_REGION

############################################
# 8. FRANKLIN S3 DESTINATIONS
############################################
S3_PATH_1="s3://genoox-upload-cdfd/samples/WES-VCF/"
S3_PATH_2="s3://genoox-upload-cdfd/samples/genome-fastq/"
S3_PATH_3="s3://genoox-upload-cdfd/samples/WES-Fastq/"
S3_PATH_4="s3://genoox-upload-cdfd/samples/WGS-VCF/"

echo "================ FRANKLIN UPLOAD =================="
echo "Uploading directory:"
echo "  ${TARGET_DIR}"
echo
echo "Choose the S3 destination:"
echo "  1) ${S3_PATH_1}"
echo "  2) ${S3_PATH_2}"
echo "  3) ${S3_PATH_3}"
echo "  4) ${S3_PATH_4}"
echo "=================================================="
read -rp "Enter the number (1-4): " PATH_OPTION

case "${PATH_OPTION}" in
  1) S3_PATH="${S3_PATH_1}" ;;
  2) S3_PATH="${S3_PATH_2}" ;;
  3) S3_PATH="${S3_PATH_3}" ;;
  4) S3_PATH="${S3_PATH_4}" ;;
  *) echo "❌ Invalid option. Exiting."; exit 1 ;;
esac

############################################
# 9. UPLOAD – STEP 3
############################################
echo
echo "🚀 Uploading to Franklin..."
echo "Local : ${TARGET_DIR}"
echo "Remote: ${S3_PATH}${TODAY}/"
echo

aws s3 cp "${TARGET_DIR}" "${S3_PATH}${TODAY}/" --recursive

echo
echo "✅ Upload complete. Pipeline finished successfully."

