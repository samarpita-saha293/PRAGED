#!/bin/bash
#srun --nodes=1 --cpus-per-task 40 --pty bash -i
#SBATCH --job-name="annotation from openmpi 4.1.4"
#SBATCH --output="EXOME_414.%j.%N.out"
#SBATCH --error="EXOME_414.%j.%N.err"
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


