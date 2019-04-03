#!/bin/bash

set -e   

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}


SAMPLE=$1
DATAPATH=$2

module load fastqc

# fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

mkdir -p $SAMPLE 

fastqc -o $SAMPLE  \
  $DATAPATH/${SAMPLE}_R1.fastq.gz  \
  $DATAPATH/${SAMPLE}_R2.fastq.gz  \
  || fail "fastqc failed"


date
echo "fastqc done"
