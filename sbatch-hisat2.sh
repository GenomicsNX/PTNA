#!/bin/bash



set -e        # stop the script if a command fails


function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}



SAMPLE=$1
DATAPATH=$2
REF=$3


module load hisat
module load samtools

mkdir -p $SAMPLE

#hisat2 -p 8 -x $REF \
#       --dta \
#       -U $DATAPATH/${SAMPLE}.R1_trimmed.fq.gz \
#       -S $SAMPLE/$SAMPLE.out.sam \
#       --rg-id $SAMPLE --rg SM:$SAMPLE 

hisat2 -p 8 -x $REF \
		--dta-cufflinks \
		-1 $DATAPATH/${SAMPLE}_R1.fastq.gz \
		-2 $DATAPATH/${SAMPLE}_R2.fastq.gz \
		-S $SAMPLE/$SAMPLE.out.sam \
		--rg-id $SAMPLE --rg SM:$SAMPLE

samtools sort -@ 8 -o $SAMPLE/$SAMPLE.out.sorted.bam $SAMPLE/$SAMPLE.out.sam

# exclude unmapped (-F 4) : mapped
# proper mapped -f 3
# exclude multiple-mapped (grep -v "XS:") 
# NH:i:n # to get the uniq : n=1

samtools view -H  $SAMPLE/$SAMPLE.out.sorted.bam > $SAMPLE/$SAMPLE.header.sam
samtools view -f 3 $SAMPLE/$SAMPLE.out.sorted.bam \
    | grep  -w "NH:i:1" | cat $SAMPLE/$SAMPLE.header.sam - \
    | samtools view -b - > $SAMPLE/$SAMPLE.unique.bam

samtools index $SAMPLE/$SAMPLE.unique.bam


rm $SAMPLE/$SAMPLE.header.sam
rm $SAMPLE/$SAMPLE.out.sam 


date
echo "hisat2 done"




