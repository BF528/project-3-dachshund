# author: Sheila Yee

#!/bin/bash -l

#$ -P bf528
#$ -cwd
#$ -pe omp 16

source /etc/bashrc

# load STAR into the environment
module load gcc star/2.6.0c

# GENOMEDIR is where the genome index is stored
GENOMEDIR=/project/bf528/project_3/reference/rn4_STAR

# perform STAR alignments against the rat genome
# readFilesCommand zcat is necessary since the fastq files are compressed/gzipped as fastq.gz
# readFilesIn $1 $2 signifies that for all of the fastq files, go through read1 and then read2
 
STAR --genomeDir $GENOMEDIR --runThreadN 16 --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --readFilesIn $1 $2 \
  --outFileNamePrefix ./STAR_output/$3 \
  --outFilterType BySJout \
  --alignSJDBoverhangMin 1 \
  --alignSJoverhangMin 8 \
 

