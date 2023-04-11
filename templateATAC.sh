#!/bin/bash
#SBATCH --partition=%partition%
#SBATCH --nodes=%nodes%
#SBATCH --cpus-per-task=20
#SBATCH --job-name=%jobname%
#SBATCH --output=%outputlog%
#SBATCH --error=%errorlog%
#SBATCH --mail-user=%mail%
#SBATCH --mail-type=FAIL
#SBATCH --time=24:00:00
#SBATCH --mem=%mem%



#################
#1-FASTQC OF ORIGINAL DATA
#################
echo "1-FASTQC start: $(date)"

%modulefastqc%
fastqc %inputfiles% -o %outputresults%/01fastqc/

echo "1-FASTQC finishes: $(date)"

#################
#2-TRIMMOMATIC OF ORIGINAL DATA
#################
echo "2-TRIMMING starts: $(date)"

%moduletrimmomatic%

# With primer
#java -jar $TRIMMOMATIC_JAR %layout% -threads 20 -phred33 %inputfiles% %outputtrim% ILLUMINACLIP:%adapters%:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

# Without primer
java -jar $TRIMMOMATIC_JAR %layout% -threads 20 -phred33 %inputfiles% %outputtrim% LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "2-TRIMMING finishes: $(date)"

################
#3-FASTQC OF TRIMMED FASTQ
################
echo "3-TRIMMING FASTQC starts: $(date)"

fastqc %trimmed%  -o %outputresults%/03trimmed_fastqc/

echo "3-TRIMMING FASTQC finishes: $(date)"

################
#4-MAQPPING WITH BWA
###############

#echo "4-MAPPING starts: $(date)"
#%modulebwa%
#%modulesamtools%
################
###4.1.-MAPPING
###############
#bwa mem -t 20 %referencedir% %trimmed% | samtools sort -@20 -o %outputresults%/04bam/%fastqid%.sorted.bam - 


################
###4.3.-MAPPING QC
###############
#samtools flagstat  %outputresults%/04bam/%fastqid%.sorted.bam > %outputresults%/04bam/%fastqid%.flagstat.txt

#echo "4-MAPPING finishes: $(date)"

################
#4-MAQPPING WITH bowtie2
###############
echo "4-MAPPING starts: $(date)"
%modulebowtie2%
%modulesamtools%

## Indexing
#bowtie2-build home/alghurmf/bulkATAC/reference_genome/hg38.fa.gz bowtie

################
###4.1.-MAPPING
###############
#bowtie2 --very-fast-local -x %referencedir% %trimmed% -S %outputresults%/04bam/%fastqid%.sam
bowtie2 -X 1000 --no-discordant --no-mixed --very-sensitive -x %referencedir% %trimmed%  -S %outputresults%/04bam/%fastqid%.sam
#bowtie2 -p 20 -X 1000 --no-discordant --no-mixed --very-sensitive -x %referencedir% %trimmed% -U ${outputresults}/02trimmed/${fastqid}_U1.fastq.gz,${outputresults}/02trimmed/${fastqid}_U2.fastq.gz -S %outputresults%/04bam/%fastqid%.sam

################
###4.2.-SAM to BAM 
###############
samtools view -S -b %outputresults%/04bam/%fastqid%.sam > %outputresults%/04bam/%fastqid%.bam

################
###4.3.-Sorting
###############
samtools sort %outputresults%/04bam/%fastqid%.bam -o %outputresults%/04bam/%fastqid%.sorted.bam

################
###4.4.-Index Mapped reads 
###############
samtools index %outputresults%/04bam/%fastqid%.sorted.bam

echo "4-MAPPING finishes: $(date)"

################
#5-Mark duplicate by picard
###############
echo "5-Mark duplicate starts: $(date)"
%modulepicard%
################
###5.1.-Mark duplicate
###############
picard MarkDuplicates I=%outputresults%/04bam/%fastqid%.sorted.bam  O=%outputresults%/04bam/%fastqid%.markdup.bam METRICS_FILE=%outputresults%/04bam/%fastqid%.markdup.txt REMOVE_DUPLICATES=FALSE ASSUME_SORTED=true CREATE_INDEX=true

################
###5.2.-Mark duplicate QC
###############
samtools flagstat %outputresults%/04bam/%fastqid%.markdup.bam > %outputresults%/04bam/%fastqid%.flagstat.markdup.txt
echo "5-Mark duplicate finishes: $(date)"

################
#6-Filtering to remove:
###############
echo "6-Filtering starts: $(date)"

################
###6.1.-Remove pcr duplicates, reads unmapped, mapq quality below 10 and blacklist 
###############

samtools view -b -F 4 -F 1024 -q 10  %outputresults%/04bam/%fastqid%.markdup.bam | samtools view -b -o %outputresults%/04bam/%fastqid%.blacklist.bam -U %outputresults%/04bam/%fastqid%.filtered.bam -L %blacklist%

################
###6.2.-Index filtered bam file
###############
samtools index %outputresults%/04bam/%fastqid%.filtered.bam

################
###6.3.-Keep just known chromosomes, remove random, unknown and chrM
###############
samtools view -b %outputresults%/04bam/%fastqid%.filtered.bam  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > %outputresults%/04bam/%fastqid%.clean.bam

################
###6.4.-Index filtered bam file
###############
samtools index %outputresults%/04bam/%fastqid%.clean.bam

################
###6.5.-Picard statistics
###############
#picard MarkDuplicates I=%outputresults%/04bam/%fastqid%.clean.bam O=%outputresults%/04bam/%fastqid%.clean.markdup.bam REMOVE_DUPLICATES=FALSE M=%outputresults%/04bam/%fastqid%.clean.markdup.txt 

################
###6.6.-samtools flagstats QC
###############
#samtools flagstat %outputresults%/04bam/%fastqid%.clean.markdup.bam >  %outputresults%/04bam/%fastqid%.clean.markdup.flagstat.txt
samtools flagstat %outputresults%/04bam/%fastqid%.clean.bam >  %outputresults%/04bam/%fastqid%.clean.flagstat.txt


echo "6-Filtering finishes: $(date)"


################
#7-Peack calling by MACS2
###############
echo "7-Peack calling starts: $(date)"
%modulemacs2%
%modulebedtools%

macs2 callpeak -t %outputresults%/04bam/%fastqid%.clean.bam -n %fastqid% --outdir %outputresults%/05PeakCalling/ -g hs -f %layoutpeak% -q 0.01 --nomodel --shift -100 --extsize 200 -B --SPMR 
total_reads=$(samtools view -c %outputresults%/04bam/%fastqid%.clean.bam)

################
###7.1.-FRiP score QC
###############
# reads in peaks

reads_in_peaks=$(bedtools sort -i %outputresults%/05PeakCalling/%fastqid%_peaks.narrowPeak \
  | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
  -a %outputresults%/04bam/%fastqid%.clean.bam -b stdin -ubam | samtools view -c)
  
# FRiP score QC
FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
echo -e "${total_reads} \t ${reads_in_peaks} \t ${FRiP}" >> %outputresults%/05PeakCalling/%fastqid%.frip.txt

echo "7-Peack calling finishes: $(date)"

################
#8-bigwig for visualization
###############
echo "9-bigwig starts: $(date)"

%modulebamCovarage%
bamCoverage --binSize 20 --normalizeUsing RPKM --effectiveGenomeSize 2913022398 -b  %outputresults%/04bam/%fastqid%.clean.bam -of bigwig -o %outputresults%/04bam/%fastqid%.coverage.bw

echo "9-bigwig finishes: $(date)"

################
#9-MULTIQC
###############
echo "8-MULTIQC starts: $(date)"

%modulemultiqc%
multiqc %outputresults% -o %outputresults%

echo "8-MULTIQC finishes: $(date)"

