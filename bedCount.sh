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

#slurm configuration
partition="batch"
mem="100G"
nodes="1"
mail="mashael.alghuraybi@kaust.edu.sa"
jobname=bedCount
outputlog=${Out_Dir}/bedCount.out
errorlog=${Out_Dir}/bedCount.err

#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 10_bedCount.sh            #
#***********************************#



module load bedtools/2.29.0

# 10. Count the number of reads mapping to each feature (OCR)

#filenames= "/ibex/scratch/alghurmf/bulAtac/*/04bam/*.clean.bam | grep -w HSC\|CLP\|proB\|preB\|PreProB\|ImmatureB\|neutrophils\|PC"

# for all cell types: 
filenames=$(ls /ibex/scratch/alghurmf/bulkAtac/*/04bam/*.clean.bam | grep "HSC\|CLP\|proB\|preB\|PreProB\|ImmatureB\|neutrophils\|PC")
bed_file="/ibex/scratch/alghurmf/bulkAtac/BedCount/Allconsensus_cell_types.bed"

Out_Dir="/ibex/scratch/alghurmf/bulkAtac/BedCount"   # path of output folder

##/bedtools2/bin/bedtools multicov -bams $filenames -bed $bed_file -D > $Out_Dir/count_table.txt
bedtools multicov -bams $filenames -bed $bed_file  -D > $Out_Dir/count_table.txt


# for neutrophils only
#filenames=$(ls /ibex/scratch/alghurmf/bulkAtac/*/04bam/*.clean.bam | grep "neutrophils")
#bed_file="/ibex/scratch/alghurmf/bulkAtac/BedCount/consensusReduce_Nutrophils.bed"

#Out_Dir="/ibex/scratch/alghurmf/bulkAtac/BedCount"   # path of output folder

#bedtools multicov -bams $filenames -bed $bed_file  -D > $Out_Dir/count_table_nuetrophils.txt
