#!/bin/bash

# Remember to activate your environment.

conda activate rnaseq

# Unload system perl for compatability with entrez-direct.

module unload perl

# Variables.

BIOPROJECT='PRJNA388952'

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

# sample 500k reads.

for FASTQ in $(find ./sequences -name "*\.fastq"); do
  NEWFILE=sequences/$(basename $FASTQ .fastq)_sampled.fastq
  seqtk sample -s 100 $FASTQ 500000 > $NEWFILE
done
