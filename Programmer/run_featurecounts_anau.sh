#!/bin/bash
#$ -P bf528
#$ -cwd
#$ -pe omp 16
source /etc/bashrc
#$ -o run_featurecounts_anau.stdout
#$ -e run_featurecounts_anau.stderr

# This script will run featureCounts for samples stored in specified folder.

# To run: qsub run_featurecounts_anau.sh
# Check status: qstat -u anau

# Star output files:
# /projectnb/bf528/users/dachshund/project_3/samples/STAR_output


# Load module:
module load subread/1.6.2

# GTF annotation file:
GTF=/project/bf528/project_3/reference/rn4_refGene_20180308.gtf

# For every applicable file in our folder:
for file in /projectnb2/bf528/users/dachshund/project_3/samples/STAR_output/*Aligned.sortedByCoord.out.bam

do

# Create a shortened filename:
f=$(basename “$file”)
f=${f:0:10}

# Run the following command:
# Some of the default arguments used:
## minOverlap: 1 overlapping base by default
echo $f
OUT=/projectnb2/bf528/users/dachshund/project_3/programmer/featureCounts_ourSamples/featureCounts_output$f.txt
featureCounts -T 16 -a $GTF -o $OUT $file

# End for loop
done;