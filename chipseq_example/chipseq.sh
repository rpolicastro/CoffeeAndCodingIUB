#!/bin/bash

# Remember to activate your environment.

conda activate chipseq

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

FASTQS=$(find ./sequences -name "*sampled*fastq")

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

# Make a directory to store the aligned reads.

mkdir -p results/aligned

# Align the reads.

for ACCESSION in ${ACCESSIONS[@]}; do
  R1=sequences/${ACCESSION}_1_sampled.fastq
  R2=sequences/${ACCESSION}_2_sampled.fastq

  bowtie2 \
    -x genome/index/scer \
    -1 $R1 \
    -2 $R2 \
    --no-mixed \
    --no-unal \
    -p 4 \
    -S results/aligned/${ACCESSION}.sam
done

# Convert sams to coordinate sorted and indexed bam.

for ACCESSION in ${ACCESSIONS[@]}; do
  samtools sort \
    -o results/aligned/${ACCESSION}.bam \
    -O bam \
    -@ 4 \
    results/aligned/${ACCESSION}.sam

  samtools index results/aligned/${ACCESSION}.bam
done

#############################
## Peak Calling with MACS2 ##
#############################

# Make a directory to store the results.

mkdir -p results/peaks

# Call peaks using MACS2.

macs2 callpeak \
  -t results/aligned/SRR1947819.bam \
  -c results/aligned/SRR1947777.bam \
  -f BAMPE \
  -g 1.2e+7 \
  --outdir results/peaks \
  -n Reb1_30s

##########################
## Bigwig File Creation ##
##########################

# Create an output directory for bigwig files.

mkdir -p results/bigwigs

# Convert to bigwigs using deeptools.

for ACCESSION in ${ACCESSIONS[@]}; do
  bamCoverage \
    -b results/aligned/${ACCESSION}.bam \
    -o results/bigwigs/${ACCESSION}.bigwig \
    -of bigwig \
    -p 4 \
    --normalizeUsing CPM
done

########################
## Generate a Heatmap ##
########################

# Download the genome annotation.

curl $ANNOTATION | gunzip > ./genome/annotation.gtf

# Create a bed file from the annotation file.

conda activate agat

agat_convert_sp_gff2bed.pl \
  --gff genome/annotation.gtf \
  -o genome/annotation.bed

conda deactivate

# Create an output directory for matrix and heatmap.

mkdir -p results/heatmaps

# Compute a count matrix.

computeMatrix reference-point \
  -R genome/annotation.bed \
  -S results/bigwigs/SRR1947819.bigwig results/bigwigs/SRR1947777.bigwig \
  -o results/heatmaps/Reb1_30s.matrix \
  --referencePoint TSS \
  -b 500 \
  -a 500 \
  -bs 50 \
  --sortRegions descend \
  --sortUsingSamples 1 \
  --missingDataAsZero \
  --samplesLabel Reb1_30s Free_MNase_30s \
  -p 4

# Create the heatmap.

plotHeatmap \
  -m results/heatmaps/Reb1_30s.matrix \
  -o results/heatmaps/Reb1_30s.pdf \
  --colorMap viridis \
  --heatmapHeight 12 \
  --heatmapWidth 3 \
  --plotFileFormat pdf
