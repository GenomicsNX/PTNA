#!/bin/sh -l
#$ -cwd

#SBATCH --job-name="StringTie"
#SBATCH --partition=largemem
#SBATCH --time=12:00:00
#SBATCH --mem=12g
#SBATCH --cpus-per-task=4

. /data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/rnaseq_box/config.sh

SAMPLEID=$1
BAMFILE=$2
OUTPUTDIR=$3

#assemble transcripts non-novel mode only
stringtie -p $NUMCPUS \
  -o $OUTPUTDIR/$SAMPLEID.gtf \
  -G /data/NHLBI_BCB/Fayaz/02-RNA-seq-transcript-level/data_source/GENCODE_GRCh38_p7_human/gencode.v25.annotation.mod.chrname.gtf \
  -l $SAMPLEID \
  -C $OUTPUTDIR/$SAMPLEID.coverage \
  -b $OUTPUTDIR \
  -e \
  -j 5 \
  -c 5 \
  --rf \
  -v \
  $BAMFILE

