# PRAGED
__A complete guide for analysing PRAGED data on a routine basis__

Please use the latest versions of these scripts on HPC. All the scripts are stored in the  below mentioned directory:
```
/home/diag_ashwin/NGS_Databases_Tools/scripts
```

## WGS and WES data Preprocessing

__Please compile in order:__
1. [MissionProjectFileRename.sh](MissionProjectFileRename.sh) - The renaming script enforces filename conventions according to PRAGED standards.
2. [MissionProjectFileTransfer.sh](MissionProjectFileTransfer.sh) - This script moves files into Mission project directories and separates PRAGED data from unrelated files.
3. [exomescript.sh](exomescript.sh)/[genomescript.sh](genomescript.sh) - Performs WES and WGS data analysis

__In place of exome_script you may use the following script as per your requirement:__

- [genomescript_withadapterseq.sh](genomescript_withadapterseq.sh) -
- [genomescript_withouttrimgalore.sh](genomescript_withouttrimgalore.sh) -
- [MissionExomeBatchRun.sh](MissionExomeBatchRun.sh) -
- [MissionExomeBatchRunwithadapterseq.sh](MissionExomeBatchRunwithadapterseq.sh) -
- [MissionExomeBatchRun_withouttrimgalore.sh](MissionExomeBatchRun_withouttrimgalore.sh) -

## Data upload to AWS

4. [awstransfer.sh](awstransfer.sh) - Transferring processed WES and WGS data

__The files to be uploaded for each sample are:__

- *.fastq.gz
- *recal.bam
- *recal.bai
- *.vcf
- *.tsv (Annotation files from Annovar)

## Data upload to Franklin


5. [aws-franklin.sh](aws-franklin.sh) - Uploads Fastq and vcf files to Franklin.



