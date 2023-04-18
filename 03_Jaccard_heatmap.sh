
#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 11_ccard_heatmap.sh   #
# script: 11_jaccard_heatmap.R      #
#***********************************#

# 9. Quality Control of Consensus Peaks - Measuring dataset similarity
# Similarity between cell types


#Input directory BAM files
In_Dir="project/consensus/bed"
#Output directory
Out_Dir="project/consensus/bed"

# Measuring dataset similarity (Jaccard)
# https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp6--measuring-dataset-similarity


file_labels=`ls $In_Dir | grep .bed`


    echo name" "$file_labels > $Out_Dir/pairwise_jaccard.txt
    for file1 in `ls $In_Dir | grep .bed`
    do
        echo -n $file1 >> $Out_Dir/pairwise_jaccard.txt

        for file2 in `ls $In_Dir | grep .bed`;
        do
            echo $file2
            # compute the jaccard stat for these two files.
            jaccard=`/opt/homebrew/bin/bedtools jaccard -a $In_Dir/$file1 -b $In_Dir/$file2`
            echo -n $jaccard
            # report the jaccard stat for these two files
            value_jaccard=$(echo $jaccard | cut -d " " -f 7)
            echo -n " "$value_jaccard >> $Out_Dir/pairwise_jaccard.txt
        done
        echo "\n" >> $Out_Dir/pairwise_jaccard.txt
    done

