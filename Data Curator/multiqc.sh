#!bin/bash -l

#$ -P bf528
#4 -cwd

source /etc/bashrc

module load python2
module load multiqc

multiqc STAR_output/


