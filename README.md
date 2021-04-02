# Project Description

A brief description of what this repository is for and what it contains

# Contributors

Data Curator: Sheila Yee
Programmer: Allison Nau
Analyst: Abhishek Thakar
Biologist: Mae Rose Gott

# Repository Contents

Provide a brief description of each script/code file in this repo, what it does, and how to execute it

### Programs from Data Curator, Sheila Yee:

### Programs from Programmer, Allison Nau:

#### run_featurecounts_anau.sh #### 
Counts genes by running featureCounts for “*Aligned.sortedByCoord.out.bam” files stored in 
/projectnb2/bf528/users/dachshund/project_3/samples/STAR_output
with the reference GTF annotation file stored in:
/project/bf528/project_3/reference/rn4_refGene_20180308.gtf

To run, submit as a job on the cluster:
```
qsub run_featurecounts_anau.sh
```

#### run_multiqc_anau.sh ####
Runs the quality control multiqc on feature counts stored in:
/projectnb2/bf528/users/dachshund/project_3/programmer/featureCounts_ourSamples

To run, submit as a job on the cluster:
```
qsub run_multiqc_anau.sh
```

#### make_csv_anau.R  ####
Combines featureCounts_output files into one csv file. Feature count files must be in the current working directory. If packages are note already installed, within R script change install packages to TRUE. Recommend running through Rstudio. To run through command line instead:
```
module load R/4.0.2
Rscript make_csv_anau.R
```

#### count_dist_anau.R ####
Creates boxplot of count distribution for samples. “featureCounts_combined.csv” must already be made using make_csv_anau.R and in the current working directory. Meta data “toxgroup_1_rna_info.csv” must also be in current working directory. Recommend running through Rstudio. To run through command line instead:
```
module load R/4.0.2
Rscript count_dist_anau.R
```

#### deseq2_anau.R ####
Runs differential expression analysis using DESeq2. Experimental count file “featureCounts_combined.csv”, control count file “control_counts.csv”, and meta data file “toxgroup_1_rna_info.csv” must be in current working directory. If packages require installation, within R script change install packages to TRUE. Recommend running through Rstuido. To run through command line instead:
```
module load R/4.0.2
Rscript deseq2_anau.R
```



### Programs from Analyst, Abhishek Thakar:

### Programs from Biologist, Mae Rose Gott:
