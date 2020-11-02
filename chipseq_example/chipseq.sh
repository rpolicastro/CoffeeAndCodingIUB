#!/bin/bash

# Remember to activate your environment.

conda activate rnaseq

# Unload system perl for compatability with entrez-direct.

module unload perl

# Variables.

BIOPROJECT='PRJNA279980'
ASSEMBLY='ftp://ftp.ensembl.org/pub/release-101/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz'
ANNOTATION='ftp://ftp.ensembl.org/pub/release-101/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.101.gtf.gz'

##########################
## Retrieve fastq Files ##
##########################

# Get csv with accession numbers.

esearch -db bioproject -query $BIOPROJECT | \
  elink -target sra | \
  efetch -format docsum | \
  xtract -pattern DocumentSummary -element Experiment@name,Run@acc |
  egrep "Reb1.*30s_1|FLAG.*30s" > \
  run_info.tsv

# Make a directory to download the fastq files to.

mkdir -p sequences

# Download fastq files.

ACCESSIONS=($(cut -f2 run_info.tsv))

for ACCESSION in ${ACCESSIONS[@]}; do
  fasterq-dump -S -e 8 -O ./sequences $ACCESSION
done

# Sample 1m reads.

for ACCESSION in ${ACCESSIONS[@]}; do
  R1=sequences/${ACCESSION}_1.fastq
  R2=sequences/${ACCESSION}_2.fastq

  seqtk sample -s100 $R1 1000000 > sequences/${ACCESSION}_1_sampled.fastq
  seqtk sample -s100 $R2 1000000 > sequences/${ACCESSION}_2_sampled.fastq
done

###########################
## fastq Quality Control ##
###########################

# Make output directory for quality control files.

mkdir -p results/fastqc_reports

# Run fastqc.

FASTQS=$(find ./sequences -name "*sampled\.fastq")

fastqc -o results/fastqc_reports $FASTQS

###################################
## Generate Bowtie2 Genome Index ##
###################################

# Make a directory to store the genome files.

mkdir -p genome

# Download and unpack the genome assembly.

curl $ASSEMBLY | gunzip > ./genome/assembly.fasta

# Generate the Bowtie2 index.

mkdir -p genome/index

bowtie2-build --threads 4 -f genome/assembly.fasta genome/index/scer

##############################
## Align Reads with Bowtie2 ##
##############################


