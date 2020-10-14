#!/bin/bash

# Remember to activate your environment.

conda activate rnaseq

# Unload system perl for compatability with entrez-direct.

module unload perl

# Variables.

BIOPROJECT='PRJNA388952'

##########################
## Retrieve fastq Files ##
##########################

# Get csv with accession numbers.

esearch -db bioproject -query $BIOPROJECT | \
elink -target sra | \
efetch -format runinfo > \
run_info.csv

# Select two example Dmel files.

DMEL=($(grep 'melanogaster' run_info.csv | cut -f1 -d, | head -n 2))

# Make a directory to download the fastq files to.

mkdir -p sequences

# Download fastq files.

for SRA in ${DMEL[@]}; do
  fasterq-dump -O sequences $SRA
done

# Sample 100k reads.

FASTQS=($(find ./sequences -name "*\.fastq"))

for FASTQ in ${FASTQS[@]}; do
  NEWFILE=sequences/$(basename $FASTQ .fastq)_sampled.fastq
  seqtk sample $FASTQ 100000 > $NEWFILE
done

###########################
## fastq Quality Control ##
###########################

# Make output directory for quality control files.

mkdir -p results/fastqc_reports

## Run fastqc.

FASTQS=$(find ./sequences -name "*sampled\.fastq")

fastqc -o results/fastqc_reports $FASTQS
