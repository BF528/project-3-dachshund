#!/bin/bash -l

#$ -P bf528
#$ -cwd
#$ -pe mpi_16_tasks_per_node 16

echo "Running job $JOB_ID"
echo "Job started: $(date +%F)"
echo "Currently running in directory: $PWD"

module load fastqc

cd /projectnb/bf528/users/dachshund/project_3/samples

# run fastqc on all of the fastq files and specify output name
fastqc *.fastq.gz

echo "Job finished: $(date +%F)"
