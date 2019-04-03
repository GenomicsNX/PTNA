#!/bin/bash 


# load global variables
# source ../globs.sh


REF="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38/genome"


for i in {1..63}; do
    SAMPLE="S$i"
	echo $SAMPLE
    sbatch --cpus-per-task=8 \
	   --time=96:00:00 \
	   --mem=16g \
	   -J $SAMPLE \
	   -o $SAMPLE.hisat.out \
	   sbatch-hisat2.sh $SAMPLE  /data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/01-fastqs   $REF

done



