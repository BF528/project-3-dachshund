# Install packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")


# Import libraries:
library(tidyverse)
library(DESeq2)
library("apeglm")

# Read in csv as dataframes:
temp_ds_exp <- read.csv("featureCounts_combined.csv", header=TRUE)
temp_ds_ctrl <- read.csv("control_counts.csv", header=TRUE)

# Convert to tibbles:
temp_tib_exp <- as_tibble(temp_ds_exp)
temp_tib_ctrl <- as_tibble(temp_ds_ctrl)

#Add prefix to control samples:
colnames(temp_tib_ctrl) <- paste("ctrl", colnames(temp_tib_ctrl), sep=".")
#Remove prefix from Geneid col name:
names(temp_tib_ctrl)[names(temp_tib_ctrl)=="ctrl.Geneid"] <- "Geneid"

# Join the two tibbles, keeping all the rows in the sample dataset (tib1)
# but not all rows in the control dataset (tib2)
#TODO: do full_join instead?
tib_comb1 <- left_join(temp_tib_exp, temp_tib_ctrl, by="Geneid")

# filter out rows that have any zeros for funzies
tib_comb2 <- subset(tib_comb1,rowSums(tib_comb1==0)==0)

#Read in Metadata dataframe:
# TODO: replace with our actual file:
meta <- read.csv("group_EX_rna_info.csv", header=TRUE)
# Convert mode_of_action to factor variable:
meta$mode_of_action <- factor(meta$mode_of_action)

# Keep only the samples in meta:
# TODO: modify for the multiple groups!!!!
control_list <- meta$Run[meta$mode_of_action=="Control"] 
print("Control sample list:")
control_list
exp_list <- meta$Run[meta$mode_of_action!="Control"]
print("Experimental sample list:")
exp_list
full_list <- meta$Run

# Subset big dataframe and only keep ones I want
tib_comb_subset <- subset(tib_comb2, select=c("Geneid", full_list))

# Convert Geneids to rownames:
tib_comb_subset <- column_to_rownames(tib_comb_subset, var="Geneid")


# create the DESeq object
dds1 <- DESeqDataSetFromMatrix(
  countData = tib_comb_subset,
  colData = meta,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds1$mode_of_action <- relevel(dds1$mode_of_action, ref='Control')

# run DESeq
dds1 <- DESeq(dds1)
# Adjust for our actual samples:
res1 <- results(dds1, contrast=c('mode_of_action','HMGCOA','Control'))
res1 <- lfcShrink(dds1, coef=2)

# write out DE results
write.csv(res1,'example_deseq_results1.csv')

# write out matrix of normalized counts
write.csv(counts(dds1,normalized=TRUE),'example_deseq_norm_counts1.csv')

