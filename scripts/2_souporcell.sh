#!/bin/bash
#
#SBATCH --job-name=demux
#SBATCH --export=All 
#SBATCH -t 24:00:00
#SBATCH --mem=100G
#SBATCH -N 1
#SBATCH -c 16   
#SBATCH --gres=tmpspace:100G 

## snp-based demultiplexing 10X samples using souporcell 
# folder structure and file location matter for the souporcell to work 
# 10X files (bam and barcodes.tsv), souporcell.sif and genome file need to be within the (sub)folder where script is run. 

## example:
# sbatch souporcell.sh MER-PC-g012 souporcell/files 2

VAR=${1:?"Specify sample name"}
VAR=${2:?"Specify folder location of the files"} 
VAR=${3:?"Specify nr of expected clusters (i.e nr pooled samples per library)"}


singularity exec -B $PWD souporcell/souporcell_latest.sif souporcell_pipeline.py \
-i $2/$1-bam.bam \
-b $2/$1-filtered-feature-bc-matrix/barcodes.tsv.gz \
-f $2/genome.fa \
-t 8 \
-o $1 \
-k $3