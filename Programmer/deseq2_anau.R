# Runs deseq2_anau.R

# To run from command line:
# module load R/4.0.2
# Rscript deseq2_anau.R

# Install packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")

#Install hardhat:
install.packages("hardhat", repos = "http://cran.us.r-project.org")

# Import libraries:
library(tidyverse)
library(DESeq2)
library("apeglm")
library(hardhat)
library(stringr)

# Read in csv as dataframes:
temp_ds_exp <- read.csv("featureCounts_combined.csv", header=TRUE)
temp_ds_ctrl <- read.csv("control_counts.csv", header=TRUE)

# Convert to tibbles:
temp_tib_exp <- as_tibble(temp_ds_exp)
temp_tib_ctrl <- as_tibble(temp_ds_ctrl)

##To add prefix to control samples:
#colnames(temp_tib_ctrl) <- paste("ctrl", colnames(temp_tib_ctrl), sep=".")
##Remove prefix from Geneid col name:
#names(temp_tib_ctrl)[names(temp_tib_ctrl)=="ctrl.Geneid"] <- "Geneid"

# Join the two tibbles, keeping all the rows that are present in both controls and exp
tib_comb1 <- inner_join(temp_tib_exp, temp_tib_ctrl, by="Geneid")

# filter out rows that have any zeros for funzies
tib_comb2 <- subset(tib_comb1,rowSums(tib_comb1==0)==0)

#Read in Metadata dataframe:
meta <- read.csv("toxgroup_1_rna_info.csv", header=TRUE)
# Convert mode_of_action to factor variable:
meta$mode_of_action <- factor(meta$mode_of_action)

# Keep only the samples in meta:
control_list <- meta$Run[meta$mode_of_action=="Control"] 
print("Control sample list:")
control_list
exp_list <- meta$Run[meta$mode_of_action!="Control"]
print("Experimental sample list:")
exp_list
full_list <- meta$Run
full_list

# Grab mode of action factor levels:
modes_of_action <- get_levels(meta)$mode_of_action
modes_of_action
exp_modes_of_action <- modes_of_action[modes_of_action!="Control"]
exp_modes_of_action

# Subset big dataframe and only keep ones I want
tib_comb_subset <- subset(tib_comb2, select=c("Geneid", full_list))

# Convert Geneids to rownames:
tib_comb_subset <- column_to_rownames(tib_comb_subset, var="Geneid")

# Run DESeq for each group:
for (exp_group in exp_modes_of_action){
  cat("===============================================================\n\n\n")
  print("Experimental group:")
  print(exp_group)
  
  # Get the experimental group sample names:
  temp_exp_sample_list <- meta$Run[meta$mode_of_action==exp_group]
  print("Sample list:")
  print(temp_exp_sample_list)
  
  # Subset dataframe:
  temp_tib <- tib_comb_subset[,c(temp_exp_sample_list, control_list)]
  # Names of columns remaining:
  print("Columns kept:")
  print(names(temp_tib))
  #Subset meta data:
  temp_meta <- meta[meta$mode_of_action==exp_group | meta$mode_of_action=="Control",]
  print("Meta data kept:")
  print(temp_meta)
  
  # create the DESeq object
  dds1 <- DESeqDataSetFromMatrix(
    countData = temp_tib,
    colData = temp_meta,
    design= ~ mode_of_action
  )
  
  # relevel mode_of_action as factor
  dds1$mode_of_action <- relevel(dds1$mode_of_action, ref='Control')
  
  # run DESeq
  dds1 <- DESeq(dds1)
  # Adjust for our actual samples:
  res1 <- results(dds1, contrast=c('mode_of_action', exp_group,'Control'))
  res1 <- lfcShrink(dds1, coef=2)

  # Strip out weird characters from exp groups:
  exp_group_str <- str_replace_all(exp_group, "[[:punct:]]", "")
  print(exp_group)
  print("Striped group name:")
  print(exp_group_str)
  
  # write out DE results
  write.csv(res1, paste(exp_group_str, 'deseq_results1.csv', sep="_"))
  
  # write out matrix of normalized counts
  write.csv(counts(dds1,normalized=TRUE),
            paste(exp_group_str, 'deseq_norm_counts1.csv', sep="_"))
}







