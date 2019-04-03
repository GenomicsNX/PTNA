#!/bin/bash 


# load global variables
# source ../globs.sh

module load fastqc

function do_fastqc
{
for i in {1..63}; do
    SAMPLE="S$i"
	echo $SAMPLE
	sbatch --mem=10g -J $SAMPLE \
	       -o $SAMPLE.log \
	       sbatch-fastqc.sh $SAMPLE \
	       /data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/01-fastqs
done
}

do_fastqc


function do_unzip
{
	# unzip them
	while read SAMPLE; do
		cd $SAMPLE		
		unzip ${SAMPLE}.R1_val_1_fastqc.zip
		unzip ${SAMPLE}.R2_val_2_fastqc.zip
		wait
		cd -
	done < $SAMPLE_LIST

}



#do_unzip


