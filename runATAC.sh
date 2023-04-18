#!/bin/bash
#SBATCH --partition=%partition%
#SBATCH --nodes=%nodes%
#SBATCH --cpus-per-task=20
#SBATCH --job-name=%jobname%
#SBATCH --output=%outputlog%
#SBATCH --error=%errorlog%
#SBATCH --mail-user=%mail%
#SBATCH --mail-type=FAIL
#SBATCH --time=99:00:00
#SBATCH --mem=%mem%

########################
##BULK ATAC SEQ PIPELINE
########################
templatefile="/home/alghurmf/bulkATAC/scripts/templateATAC.sh"
inputdir="bulkATAQ/original_data/Nichefastq/191108_R355"   #path where is the FASTQ, new samples
outputdir="bulkATAQ/original_data/Nichefastq/191108_R355"   # path of output folder


#sample information
fastqid=$1 # if it is PE would be fastqid_1.fastq.gz, fastqid_2.fastq.gz, for SE fastqid.fastq.gz
organism="human" #human or mouse
if [[ ${organism} = "human" ]]; then
    referencedir="/home/alghurmf/bulkATAC/reference_genome/bowtie2" # indexed hg38 genome by bowtie2
	blacklist="/home/alghurmf/bulkATAC/hg38-blacklist.v2.bed"
fi


#library information
layout="PE" #PE=pair-end, SE=single-end
adapters="/home/alghurmf/bulkATAC/adapters/adapters/NexteraPE-PE.fa" #remove adapters in trimmomatic if needed

#software
fastqc="module load fastqc/0.11.8"
trimmomatic="module load trimmomatic/0.38"
bwa="module load bwa/0.7.17"
bowtie2="module load bowtie2/2.3.5"
samtools="module load samtools/1.8"
multiqc="module load multiqc/1.9"
picard="module load picard/2.20.4"
macs2="module load macs2/2.1.1.2/anaconda2-2.5.0"
bedtools="module load bedtools/2.29.0"
bamCovarage="module load deeptools/python2.7/3.3.1"

#slurm configuration
partition="batch"
mem="100G"
nodes="1"
mail="YOUR EMAIL"
jobname=${fastqid}_bulkATAC
outputlog=${outputdir}/${fastqid}/00script/${fastqid}_bulkATAC.out
errorlog=${outputdir}/${fastqid}/00script/${fastqid}_bulkATAC.err

########
#STEPS
########
###################
#1-CHECK: template file, input dir, output dir and original fastq
###################
if [ ! -d ${inputdir} ]; then
    echo "INPUT DIR: ${inputdir} NOT FOUND."
    exit 0
fi
if [ ! -d ${outputdir} ]; then
    echo "OUTPUT DIR: ${outputdir} NOT FOUND."
    exit 0
fi
if [ -d ${outputdir}/${fastqid} ]; then
    echo "OUTPUT RESULTS DIR: ${outputdir} ALREADY EXISTS."
    exit 0
fi
if [ ! -f ${templatefile} ]; then
    echo "TEMPLATE FILE: ${templatefile} NOT FOUND."
    exit 0
fi
echo "1-PATHS CHECKED SUCCESFULLY"

#################
#2-CREATE folders
#################
mkdir ${outputdir}/${fastqid}
mkdir ${outputdir}/${fastqid}/00script
mkdir ${outputdir}/${fastqid}/01fastqc
mkdir ${outputdir}/${fastqid}/02trimmed
mkdir ${outputdir}/${fastqid}/03trimmed_fastqc
mkdir ${outputdir}/${fastqid}/04bam
mkdir ${outputdir}/${fastqid}/05PeakCalling
echo "2-FOLDERS CREATED SUCCESFULLY"

################
#3-CREATE script
################
outputresults=${outputdir}/${fastqid}
outputfastqc=${outputresults}/01fastqc/
if [[ ${layout} = "PE" ]]; then
    inputfiles="${inputdir}/${fastqid}_R1_001.fastq.gz ${inputdir}/${fastqid}_R2_001.fastq.gz"
    outputtrim="${outputresults}/02trimmed/${fastqid}_trim1.fastq.gz ${outputresults}/02trimmed/${fastqid}_U1.fastq.gz ${outputresults}/02trimmed/${fastqid}_trim2.fastq.gz ${outputresults}/02trimmed/${fastqid}_U2.fastq.gz"
    trimmed="-1 ${outputresults}/02trimmed/${fastqid}_trim1.fastq.gz -2 ${outputresults}/02trimmed/${fastqid}_trim2.fastq.gz"
    trimmedU="${outputresults}/02trimmed/${fastqid}_U1.fastq.gz ${outputresults}/02trimmed/${fastqid}_U2.fastq.gz"
    layoutpeak="BAMPE"
else
    inputfiles="${inputdir}/${fastqid}.fastq.gz"
    outputtrim="${outputresults}/02trimmed/${fastqid}_trim.fastq.gz"
    trimmed=${outputtrim}
    layoutpeak="BAM"
fi



sed "s,%fastqid%,${fastqid},g
    s,%partition%,${partition},g
    s,%nodes%,${nodes},g
    s,%jobname%,${jobname},g
    s,%outputlog%,${outputlog},g
    s,%errorlog%,${errorlog},g
    s,%mail%,${mail},g
    s,%mem%,${mem},g
    s,%modulefastqc%,${fastqc},g
    s,%inputfiles%,${inputfiles},g
    s,%moduletrimmomatic%,${trimmomatic},g
    s,%outputtrim%,${outputtrim},g
    s,%layout%,${layout},g
    s,%adapters%,${adapters},g
    s,%trimmed%,${trimmed},g
    s,%trimmedU%,${trimmedU},g
    s,%modulebwa%,${bwa},g
    s,%modulesamtools%,${samtools},g
    s,%modulemultiqc%,${multiqc},g
    s,%modulepicard%,${picard},g
    s,%modulemacs2%,${macs2},g
    s,%modulebowtie2%,${bowtie2},g
    s,%modulebedtools%,${bedtools},g
    s,%blacklist%,${blacklist},g
    s,%layoutpeak%,${layoutpeak},g
    s,%referencedir%,${referencedir},g
    s,%modulebamCovarage%,${bamCovarage},g
    s,%outputresults%,${outputresults},g" ${templatefile} > ${outputresults}/00script/${fastqid}_bulkATAC.sh

echo "3-SCRIPT CREATED SUCCESFULLY"

################
#4-EXECUTE script
################
echo "4-RUN SCRIPT"
sbatch ${outputresults}/00script/${fastqid}_bulkATAC.sh
