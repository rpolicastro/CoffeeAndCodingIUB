#!/bin/bash

# Remember to activate your environment.

conda activate rnaseq

# Unload system perl for compatability with entrez-direct.

module unload perl

# Variables.

BIOPROJECT='PRJNA388952'
ASSEMBLY='ftp://ftp.ensembl.org/pub/release-101/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna_sm.toplevel.fa.gz'
ANNOTATION='ftp://ftp.ensembl.org/pub/release-101/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.101.chr.gtf.gz'

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

# Run fastqc.

FASTQS=$(find ./sequences -name "*sampled\.fastq")

fastqc -o results/fastqc_reports $FASTQS

################################
## Generate STAR Genome Index ##
################################

# Make a directory to store the genome files.

mkdir -p genome

# Download and unpack the genome assembly.

curl $ASSEMBLY | gunzip > ./genome/assembly.fasta

# Download and unpack the genome annotation.

curl $ANNOTATION | gunzip > ./genome/annotation.gtf

# Create a directory to store the index.

mkdir -p genome/index

# Create the STAR genome index.
# --genomeSAindexNbases 12 was recommended by software.

STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir genome/index \
  --genomeFastaFiles genome/assembly.fasta \
  --sjdbGTFfile genome/annotation.gtf \
  --genomeSAindexNbases 12

###########################
## Align Reads to Genome ##
###########################

# Create an output directory for aligned reads.

mkdir -p results/aligned

# Align the reads.

FASTQS=($(find ./sequences -name "*\sampled.fastq"))

for FASTQ in ${FASTQS[@]}; do
  PREFIX=results/aligned/$(basename $FASTQ .fastq)_
  STAR \
    --runThreadN 4 \
    --genomeDir genome/index \
    --readFilesIn $FASTQ \
    --outFileNamePrefix $PREFIX \
    --outSAMtype BAM SortedByCoordinate
done

# Indexing the BAM files.

BAMS=($(find ./results/aligned -name "*\.bam"))

for BAM in ${BAMS[@]}; do
  samtools index $BAM
done

####################
## Count Features ##
####################

# Create an output directory for read counts.

mkdir -p results/counts

# Count reads.

BAMS=$(find ./results/aligned -name "*\.bam")

featureCounts \
  -a genome/annotation.gtf \
  -o results/counts/counts.tsv \
  -t exon \
  -g gene_name \
  --largestOverlap \
  --readExtension3 150 \
  --primary \
  -s 0 \
  -T 4 \
  ${BAMS}
