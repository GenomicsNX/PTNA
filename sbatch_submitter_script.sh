#!/bin/sh -l
#$ -cwd

. /data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/rnaseq_box/config.sh

BAMDIR="/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/03-alignment"

while read SAMPLEID
do
	mkdir "$OUTPUTDIR_STRINGTIE_SACK_V2/$SAMPLEID"
	sbatch -J $SAMPLEID -o $SAMPLEID.log "$OUTPUTDIR_STRINGTIE_SACK_V2/run_stringtie_human_only.sh $SAMPLEID $BAMDIR/$SAMPLEID/$SAMPLEID.unique.bam $OUTPUTDIR_STRINGTIE_SACK_V2/$SAMPLEID"
done < "/data/NHLBI_BCB/Sack_Kim-Han_RNAseq_v2/04-featureCounts/sampleids.txt"

