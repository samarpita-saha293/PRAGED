#!/bin/bash
#srun --nodes=1 --cpus-per-task 40 --pty bash -i
#SBATCH --job-name="GENOME from openmpi 4.1.4"
#SBATCH --output="GENOME_414.%j.%N.out"
#SBATCH --error="GENOME_414.%j.%N.err"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=samarpitasaha@cdfd.org.in

module load FastQC
module load anaconda3-multiqc
module load bwa-0.7.17
module load samtools-1.16.1
module load gatk-4.3
module load bcftools-1.16

Fasta="/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/Homo_sapiens_assembly38.fasta"
dbsnp='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/Homo_sapiens_assembly38.dbsnp138.vcf'
GoldIndels='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
hapmap='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/hapmap_3.3.hg38.vcf.gz'
Omni25='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/1000G_omni2.5.hg38.vcf.gz'
OneKgvcf='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/1000G_omni2.5.hg38.vcf.gz'
Annovar='/home/diag_ashwin/NGS_Databases_Tools/Hg38/annovar/'
AnnovarDb="/home/diag_ashwin/NGS_Databases_Tools/Hg38/annovar/humandb_hg38/"
Scripts='/home/diag_ashwin/NGS_Databases_Tools/scripts/'

##Bed file; Depending on the target capture kit bed file need to change
Exome_bed='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/merged_probe_file_shareable_CDFD_CH-LSD_TE-91704033_hg38_Pad30bp_250110174100.bed'
Exome_Vcfbed='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/merged_probe_file_shareable_CDFD_CH-LSD_TE-91704033_hg38_Pad30bp_250110174100.bed'

#Softwares
picard='java -jar /home/diag_ashwin/NGS_Databases_Tools/picard/picard.jar'
Mosedepth_plots='python /home/diag_ashwin/NGS_Databases_Tools/Mosedepth_plots/plot-dist.py'

N=40 ; #Number of jobs can be run parallel
threads=40

#Creating Sample names based on fastq File
ls *.fastq.gz > fasqfile.txt && rm Sample_Names.txt temp.txt
cat fasqfile.txt

while read Fastq; do basename "$Fastq" | sed -E 's/_?R[12].*\.fastq\.gz//' >> temp.txt ;done < fasqfile.txt
cat temp.txt |sort| uniq |sed '/^$/d' > Sample_Names.txt
cat -n Sample_Names.txt

#####################----------------------------------Analysis started---------------------------------------------------------------------
while read -r Samplename; do ((i = i % N)); ((i++ == 0)) && wait

  unset array && array=($(ls "${Samplename}"*fastq.gz))  # Use quotes to handle spaces in filenames

  echo "Step1: Prealignment QC"

  echo "1.A) QC checking of Raw FastQ data using FastQC"
  mkdir -p QC/FastQC/Raw_Fastq_QC/"${Samplename}"
  fastqc -t ${threads} -o QC/FastQC/Raw_Fastq_QC/ "${array[0]}" & ##If fasq is not ending like '_R1.fq.gz'  it will not work
  fastqc -t ${threads} -o QC/FastQC/Raw_Fastq_QC/ "${array[1]}" & ##If fasq is not ending like '_R2.fq.gz'  it will not work

done < Sample_Names.txt
wait

mkdir -p MultiQC/FastQC/RawFastq_QC
multiqc -o MultiQC/FastQC/RawFastq_QC/ QC/FastQC/Raw_Fastq_QC/

echo "1.B) Remove low-quality reads, adaptor, and bases using TrimGalore or FASTP"

mkdir -p QC/TrimGalore.stdout_stats
while read -r Samplename; do
  ((i = i % N)); ((i++ == 0)) && wait

  unset array
  array=($(ls "${Samplename}"_*fastq.gz))

  trim_galore --paired "${array[0]}" "${array[1]}" \
  --length 14 --length_1 15 --length_2 15 \
  --retain_unpaired --basename "${Samplename}" --cores ${threads} >& QC/TrimGalore.stdout_stats/"${Samplename}"_TrimGalore.stdout_stats.txt &

done < Sample_Names.txt
wait

mkdir -p QC/FastQC/TrimGalore_Fastq_QC/

echo "1.C) FastQC of cleaned fastq file  "

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait

mkdir -p QC/FastQC/TrimGalore_Fastq_QC/
fastqc -t ${threads} -o QC/FastQC/TrimGalore_Fastq_QC/ ${Samplename}_val_1.fq.gz &
fastqc -t ${threads} -o QC/FastQC/TrimGalore_Fastq_QC/ ${Samplename}_val_2.fq.gz &
done < Sample_Names.txt ;wait

echo "1.D) Moving the TrimGalore Filtered/Trimed fastq files "

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait

mkdir -p CleanedFastQ_Files/TrimGalore_QC/
mv  ${Samplename}_val_1.fq.gz  CleanedFastQ_Files/TrimGalore_QC/
mv  ${Samplename}_val_2.fq.gz   CleanedFastQ_Files/TrimGalore_QC/

mkdir -p QC_FailedFastQ

mv ${Samplename}*{1,2}.fastq.gz_trimming_report.txt QC/TrimGalore.stdout_stats/
mv ${Samplename}*{1,2}_unpaired_{1,2}.fq.gz QC_FailedFastQ

done < Sample_Names.txt ;wait # for all the something with stuff


echo "Step2: Alignment and Post Alignment Processing"

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
   echo "2.A) Alignment ;"
ID="HTYYCBBXX.1.${Samplename}"
SM="${Samplename}"
LB=Lib-1
PU=HTYYCBBXX.1.TCAATCCG+TTCGCAGT
PL=ILLUMINA

#SM : give the name of the sample, CH_plate1_A01, the read is from,
#LB : gives the name of the library prep, Lib-1, the read is from,
#PU : gives the name of the flowcell, HTYYCBBXX, and the lane, 1, the read is from.

mkdir -p BamFiles
bwa mem  -R "@RG\tID:${Samplename}\tSM:${Samplename}\tLB:Lib\tPU:PU\tPL:Illumina" \
-M -t ${threads} \
${Fasta} \
CleanedFastQ_Files/TrimGalore_QC/${Samplename}_val_1.fq.gz \
CleanedFastQ_Files/TrimGalore_QC/${Samplename}_val_2.fq.gz | samtools sort -@ ${threads} -o BamFiles/${Samplename}.bam  \
 & done < Sample_Names.txt ; wait

#2.B)

echo "mark duplicates ; https://gatk.broadinstitute.org/hc/en-us/articles/360057438771-MarkDuplicatesSpark; https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-"

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
   #gatk MarkDuplicatesSpark \
   #-I BamFiles/${Samplename}.bam \
   #-O BamFiles/${Samplename}_Mduplicates.bam \
   #-M BamFiles/${Samplename}_marked_dup_metrics.txt \
   #--create-output-bam-index true & done < Sample_Names.txt ; wait

${picard} MarkDuplicates I=BamFiles/${Samplename}.bam \
O=BamFiles/${Samplename}_Mduplicates.bam \
M=BamFiles/${Samplename}_marked_dup_metrics.txt  & done < Sample_Names.txt ; wait

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
${picard} BuildBamIndex I=BamFiles/${Samplename}_Mduplicates.bam O=BamFiles/${Samplename}_Mduplicates.bai & done < Sample_Names.txt ; wait


#echo "2.C) Base (Quality Score) Recalibration ; https://gatk.broadinstitute.org/hc/en-us/articles/360056969412-BaseRecalibrator"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
mkdir -p Recal_Files
gatk BaseRecalibrator \
-I BamFiles/${Samplename}_Mduplicates.bam \
-R ${Fasta} \
--known-sites ${dbsnp} \
--known-sites ${GoldIndels} \
-O Recal_Files/${Samplename}_recal_data.table & done < Sample_Names.txt ; wait


echo "2.D) ApplyBQSR ; https://gatk.broadinstitute.org/hc/en-us/articles/360056968652-ApplyBQSR"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk ApplyBQSR \
-R ${Fasta} \
-I BamFiles/${Samplename}_Mduplicates.bam \
--bqsr-recal-file Recal_Files/${Samplename}_recal_data.table \
-O BamFiles/${Samplename}_recal.bam \
--create-output-bam-index true & done < Sample_Names.txt ; wait

echo "2.E) Create Post recal table #https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk BaseRecalibrator \
-I BamFiles/${Samplename}_recal.bam \
-R ${Fasta} \
--known-sites ${dbsnp} \
--known-sites ${GoldIndels} \
--output  Recal_Files/${Samplename}_post_recal_data.table & done < Sample_Names.txt ; wait

echo "2.F) AnalyzeCovariates ; https://gatk.broadinstitute.org/hc/en-us/articles/360056967752-AnalyzeCovariates"
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
gatk AnalyzeCovariates \
-before Recal_Files/${Samplename}_recal_data.table \
-after Recal_Files/${Samplename}_post_recal_data.table \
-plots Recal_Files/${Samplename}_AnalyzeCovariates.pdf & done < Sample_Names.txt ; wait

#-------------------------------------------------------------------------------------------------------------------------------------------------------------

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait

echo "3.C) Mos depth for LSD panel (ultra-deep coverage)"

mkdir -p QC/AlignmentQC/Mosdepth/panel/

mosdepth -x -t 20 --by ${genome_bed} \
QC/AlignmentQC/Mosdepth/panel/${Samplename}_panel \
BamFiles/${Samplename}_Mduplicates.bam \
--thresholds 50,100,200,300,500,800,1000 &

mkdir -p QC/AlignmentQC/flagstat
samtools flagstat -@ 35 BamFiles/${Samplename}_Mduplicates.bam > QC/AlignmentQC/flagstat/${Samplename}_flagstat.txt &

done < Sample_Names.txt ; wait

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
${Mosedepth_plots}  QC/AlignmentQC/Mosdepth/genome/*.dist.txt --output QC/AlignmentQC/Mosdepth/genome/${Samplename}_MosdepthCoverage.html
done < Sample_Names.txt ; wait

#Include AlphaMissense column next to REVEL in tsv 

mv *bamout.ba* ~/genome/haplobam-genome_all/
mv *g.vcf* ~/genome/gvcf-genome_all/
mv BamFiles/*recal* .
mv QC/FastQC/Raw_Fastq_QC QC/
mv MultiQC/FastQC/RawFastq_QC/ QC

rm -rf *.avinput *hg38_multianno.txt CleanedFastQ_Files/ GENOME_414* ${SAMPLENAME}_raw.vcf.gz* ${SAMPLENAME}_INDEL.vcf.gz* ${SAMPLENAME}_INDEL_filtered.vcf.gz* ${SAMPLENAME}_INDEL.vcf.gz ${SAMPLENAME}_Laligned.vcf ${SAMPLENAME}_snps.vcf.gz* ${SAMPLENAME}_snps_filtered.vcf.gz* Recal_Files QC_FailedFastQ CleanedFastQ_Files/ custom-out.log  fasqfile.txt  Sample_Names.txt  temp.txt  QC_FailedFastQ/ *hardfiltered_Target.vcf.idx QC/TrimGalore.stdout_stats QC/AlignmentQC 
