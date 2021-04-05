# author: Sheila Yee

#!/bin/bash -l

#$ -P bf528
#$ -cwd
#$ -o multiqc.log
#$ -m bae

source /etc/bashrc

module load python2
module load multiqc

#cd /projectnb/bf528/users/dachshund/project_3/samples

#run multiqc using folder path containing results from fastqc and STAR alignment output files
#supply multiqc with as many directories and folders to scan 
#multiqc /projectnb/bf528/users/dachshund/project_3/samples /projectnb/bf528/users/dachshund/project_3/samples/STAR_ouput

multiqc /projectnb/bf528/users/dachshund/project_3/samples/
