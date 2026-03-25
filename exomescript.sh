#!/bin/bash
#srun --nodes=1 --cpus-per-task 40 --pty bash -i
#SBATCH --job-name="EXOME from openmpi 4.1.4"
#SBATCH --output="EXOME_414.%j.%N.out"
#SBATCH --error="EXOME_414.%j.%N.err"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=pragnapcu@gmail.com

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
Exome_bed='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated.bed'
Exome_Vcfbed='/home/diag_ashwin/NGS_Databases_Tools/Hg38/GATK_hg38_Resources/hg38_exome_comp_spikein_v2.0.2_targets_sorted.re_annotated_10.bed'

#Softwares
picard='java -jar /home/diag_ashwin/NGS_Databases_Tools/picard/picard.jar'
Mosedepth_plots='python /home/diag_ashwin/NGS_Databases_Tools/Mosedepth_plots/plot-dist.py'

N=15  # Number of jobs that can run in parallel
threads=25

# Step 1: Creating sample names based on FASTQ files
ls *.fastq.gz > fasqfile.txt && rm -f Sample_Names.txt temp.txt  # Use -f to avoid errors if the files don't exist
cat fasqfile.txt

while read -r Fastq; do echo "$Fastq" | cut -d "_" -f 1 >> temp.txt ; done < fasqfile.txt

cat temp.txt | sort | uniq | sed '/^$/d' > Sample_Names.txt
cat -n Sample_Names.txt

#####################----------------------------------Analysis started---------------------------------------------------------------------

while read -r Samplename; do ((i = i % N)); ((i++ == 0)) && wait

  unset array && array=($(ls "${Samplename}"_*fastq.gz))  # Use quotes to handle spaces in filenames

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

echo "1.C) FastQC of cleaned FASTQ files"

while read -r Samplename; do
  ((i = i % N)); ((i++ == 0)) && wait

  fastqc -t "${threads}" -o QC/FastQC/TrimGalore_Fastq_QC/ "${Samplename}"_val_1.fq.gz &
  fastqc -t "${threads}" -o QC/FastQC/TrimGalore_Fastq_QC/ "${Samplename}"_val_2.fq.gz &
done < Sample_Names.txt
wait

mkdir -p MultiQC/FastQC/TrimGalore_Fastq_QC
multiqc -o MultiQC/FastQC/TrimGalore_Fastq_QC/ QC/FastQC/TrimGalore_Fastq_QC/

echo "1.D) Moving the TrimGalore Filtered/Trimmed FASTQ files"

while read -r Samplename; do
  ((i = i % N)); ((i++ == 0)) && wait

  mkdir -p CleanedFastQ_Files/TrimGalore_QC/
  mv "${Samplename}"_val_1.fq.gz CleanedFastQ_Files/TrimGalore_QC/
  mv "${Samplename}"_val_2.fq.gz CleanedFastQ_Files/TrimGalore_QC/

  mkdir -p QC_FailedFastQ

  mv "${Samplename}"*trimming_report.txt QC/TrimGalore.stdout_stats/
  mv "${Samplename}"*_unpaired_*.fq.gz QC_FailedFastQ

done < Sample_Names.txt
wait

echo "Step2: Alignment and Post Alignment Processing"

while read -r Samplename; do
  ((i = i % N)); ((i++ == 0)) && wait

  echo "2.A) Alignment"

  ID="HTYYCBBXX.1.${Samplename}"
  SM="${Samplename}"
  LB="Lib-1"
  PU="HTYYCBBXX.1.TCAATCCG+TTCGCAGT"
  PL="ILLUMINA"

  mkdir -p BamFiles
  bwa mem -R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPU:${PU}\tPL:${PL}" \
  -M -t "${threads}" \
  "${Fasta}" \
  CleanedFastQ_Files/TrimGalore_QC/"${Samplename}"_val_1.fq.gz \
  CleanedFastQ_Files/TrimGalore_QC/"${Samplename}"_val_2.fq.gz | samtools sort -@ "${threads}" -o BamFiles/"${Samplename}".bam &
done < Sample_Names.txt
wait

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

echo "3.C) Mos depth global and target region and gene exon"

mkdir -p QC/AlignmentQC/Mosdepth/Exome/


mosdepth -x --by ${Exome_bed} QC/AlignmentQC/Mosdepth/Exome/${Samplename}_Exome BamFiles/${Samplename}_Mduplicates.bam \
--thresholds 1,10,20,30,50,100 &


mkdir -p QC/AlignmentQC/flagstat
samtools flagstat BamFiles/${Samplename}_Mduplicates.bam > QC/AlignmentQC/flagstat/${Samplename}_flagstat.txt &

done < Sample_Names.txt ; wait


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
${Mosedepth_plots}  QC/AlignmentQC/Mosdepth/Exome/*.dist.txt --output QC/AlignmentQC/Mosdepth/Exome/${Samplename}_MosdepthCoverage.html
done < Sample_Names.txt ;wait



mkdir -p MultiQC/AlignmentQC/Mosdepth/Exome
multiqc -o MultiQC/QC/AlignmentQC/Mosdepth/Exome/ QC/AlignmentQC/Mosdepth/Exome


mkdir -p MultiQC/CompleteQC
multiqc -o MultiQC/CompleteQC/ QC


#
#
#
###---------------------------------------Variant calling-----------------------------------------------------------------------"
#
echo "#Step4 ;  Variant calling ; "

while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
    echo "4.A) HaplotypeCaller; gvcf creation"
 gatk --java-options "-Xmx16g -XX:ParallelGCThreads=12"  HaplotypeCaller  \
   -R ${Fasta} \
   -I BamFiles/${Samplename}_recal.bam \
   -O ${Samplename}.g.vcf.gz \
   -bamout ${Samplename}_bamout.bam \
   -ERC GVCF \
   -A StrandBiasBySample \
   --native-pair-hmm-threads 8 \
   -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation &

done < Sample_Names.txt ; wait



while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
echo "4.B) GenotypeGVCFs ; create single sample vcf file ; https://gatk.broadinstitute.org/hc/en-us/articles/360036899732-GenotypeGVCFs"
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=12"  GenotypeGVCFs \
 -R ${Fasta} \
 -V ${Samplename}.g.vcf.gz \
 -O ${Samplename}_raw.vcf.gz \
 -A StrandBiasBySample \
 -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation &
 done < Sample_Names.txt ; wait
#
###---------------------------------------Variant recalibration----------------------------------------------------------------------
#Step5 ;VariantRecalibrator  ; Not performing ; insted hard filter will do ; https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering




#normalize variants ;
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        bcftools norm \
        -f ${Fasta} \
        ${Samplename}_raw.vcf.gz \
        --multiallelics -any \
        -o ${Samplename}_Laligned.vcf -O v  &
done < Sample_Names.txt ; wait

Subset to SNPs-only callset with SelectVariants
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        gatk SelectVariants \
        -V ${Samplename}_Laligned.vcf \
        -select-type SNP \
        -O ${Samplename}_snps.vcf.gz  & done < Sample_Names.txt ; wait


Subset to indels-only callset with SelectVariants
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        gatk SelectVariants \
        -V ${Samplename}_Laligned.vcf \
        -select-type INDEL \
        -O ${Samplename}_INDEL.vcf.gz  &  done < Sample_Names.txt ; wait




i#Hard-filter SNPs on multiple expressions using VariantFiltration
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        gatk VariantFiltration \
            -V ${Samplename}_snps.vcf.gz \
           -filter "QD < 2.0" --filter-name "QD2" \
           -filter "QUAL < 30.0" --filter-name "QUAL30" \
           -filter "SOR > 3.0" --filter-name "SOR3" \
           -filter "FS > 60.0" --filter-name "FS60" \
           -filter "MQ < 40.0" --filter-name "MQ40" \
           -filter "DP < 10.0" --filter-name "DP10" \
           -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
           -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O ${Samplename}_snps_filtered.vcf.gz & done < Sample_Names.txt ; wait


while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        gatk VariantFiltration \
            -V ${Samplename}_INDEL.vcf.gz \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "DP < 10.0" --filter-name "DP10" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
           -O ${Samplename}_INDEL_filtered.vcf.gz & done < Sample_Names.txt ; wait


#STEP16: Combine SNP and INDEL
        while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
                gatk MergeVcfs  \
                -R ${Fasta} \
                -I "${Samplename}_"snps_filtered.vcf.gz \
                -I "${Samplename}_"INDEL_filtered.vcf.gz \
                 --CREATE_INDEX true \                -O "${Samplename}_"hardfiltered.vcf.gz & done < Sample_Names.txt ; wait




#STEP16: Combine SNP and INDEL
        while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
                gatk MergeVcfs  \
                -R ${Fasta} \
                -I "${Samplename}_"snps_filtered.vcf.gz \
                -I "${Samplename}_"INDEL_filtered.vcf.gz \
                --CREATE_INDEX true \
                -O "${Samplename}_"hardfiltered.vcf.gz & done < Sample_Names.txt ; wait



#Remove the intermediate files
        while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        rm ${Samplename}_Laligned.vcf ${Samplename}_snps.vcf.gz*  ${Samplename}_INDEL.vcf.gz*
        rm ${Samplename}_snps_filtered.vcf.gz* ${Samplename}_INDEL_filtered.vcf.gz* ${Samplename}_raw.vcf.gz* & done < Sample_Names.txt ; wait




#Select variants from target region
while read Samplename ; do ((i=i%N)); ((i++==0)) && wait
        gatk SelectVariants \
       -V "${Samplename}_"hardfiltered.vcf.gz \
        -L ${Exome_Vcfbed} \
        -O "${Samplename}_"hardfiltered_Target.vcf
        done < Sample_Names.txt ; wait

###ANNOTATION OF VCF FILE USING aNNOVAR
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait     

        perl ${Annovar}convert2annovar.pl \
                -format vcf4 -allsample -withfreq ${SAMPLENAME}_hardfiltered_Target.vcf \
                -outfile ${SAMPLENAME}.avinput \
                -includeinfo  &
done < Sample_Names.txt 
wait

while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait

perl ${Annovar}table_annovar.pl \
    ${SAMPLENAME}.avinput ${AnnovarDb} \
    -buildver hg38 \
    -out ${SAMPLENAME} \
    -remove \
    -otherinfo \
    -protocol refGene,gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad41_exome,gnomad40_exome,gnomad30_genome,gnomad40_genome,gnomad41_genome,kaviar_20150923,avsnp150,avsnp151,dbnsfp42a,dbnsfp47a,dbscsnv11,intervar_20250721,revel,clinvar_20250715,clinvar_20250721,genomicSuperDups,dgvMerged,gwasCatalog \
    -operation gx,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
    -nastring . \
    -polish \
    -xref ${AnnovarDb}gene_fullxref.txt

done < Sample_Names.txt
wait

##Filter annovar output file based on the defined parameters
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait
        python ${Annovar}AnnovarFiltering_v4.py  --AnnovarFile ${SAMPLENAME}.hg38_multianno.txt &
done < Sample_Names.txt ; wait

# Adding AlphaMissense to the tsv output files
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait
        python ${Annovar}alphamissense-v2.py --annovar_file ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv --filtered_file ${SAMPLENAME}.hg38_multiannoDesiredColumns_Filtered.tsv &
done < Sample_Names.txt ; wait

#Adding LLR score to the tsv output files
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait
	python ${Annovar}/ESMscore-v2.py --annovar_file ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv --filtered_file ${SAMPLENAME}.hg38_multiannoDesiredColumns_Filtered.tsv &
done < Sample_Names.txt ; wait


#Edit column headers and Incorporate additional columns
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait       
        find="Otherinfo4"
        replace=$(grep "#CHROM"  ${SAMPLENAME}_hardfiltered_Target.vcf)
        sed -i "s+${find}+${replace}+g" ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv
        sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv

        find="Otherinfo4"
        replace=$(grep "#CHROM" ${SAMPLENAME}_hardfiltered_Target.vcf)
        sed -i "s+${find}+${replace}+g" ${SAMPLENAME}.hg38_multiannoDesiredColumns_Filtered.tsv
        sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}.hg38_multiannoDesiredColumns_Filtered.tsv

done < Sample_Names.txt ; wait

##Additional Analysis
while read SAMPLENAME ; do ((i=i%N)); ((i++==0)) && wait

#Exome region based coverage %
python /home/diag_ashwin/NGS_Databases_Tools/Hg38/MosedepthTargetComparison.py \
        -Mfolder QC/AlignmentQC/Mosdepth/Exome/ \
        -Outputfolder  QC/AlignmentQC/Mosdepth/Exome/

##dEFINE gENOTYPE(rEFINED GENOTYPE ) BASED ON THE DEFINED PARAMETERS
python /home/diag_ashwin/NGS_Databases_Tools/Hg38/MosedepthTargetRegionBed_toPrecentage.py \
        -Mfolder QC/AlignmentQC/Mosdepth/Exome/

python /home/diag_ashwin/NGS_Databases_Tools/Hg38/RefinedGenotype_v2.py \
         --AnnovarOutputFile ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv

python /home/diag_ashwin/NGS_Databases_Tools/Hg38/RefinedGenotype_v2.py \
                 --AnnovarOutputFile ${SAMPLENAME}.hg38_multiannoDesiredColumns_Filtered.tsv

done < Sample_Names.txt ; wait

mv *bamout.ba* /home/diag_ashwin/exome/haplobam_all/
mv *g.vcf* /home/diag_ashwin/exome/gvcf_all/
mv BamFiles/*recal* .
mv QC/FastQC/Raw_Fastq_QC QC/
mv QC/AlignmentQC/Mosdepth QC/
mv QC/Mosdepth/Exome/* QC/Mosdepth/
rm -r *.avinput *hardfiltered.vcf.gz *hardfiltered.vcf.gz.tbi *hg38_multianno.txt CleanedFastQ_Files/ custom-out.log  fasqfile.txt  Sample_Names.txt  temp.txt  MultiQC/ QC_FailedFastQ/ Recal_Files/ *hardfiltered_Target.vcf.idx QC/TrimGalore.stdout_stats QC/Mosdepth/Exome QC/AlignmentQC QC/FastQC BamFiles EXOME*
