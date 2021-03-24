# Install packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")


# Import libraries:
library(tidyverse)
library(DESeq2)

# Read in csv as dataframes:
temp_ds_exp <- read.csv("featureCounts_combined.csv", header=TRUE)
temp_ds_ctrl <- read.csv("control_counts.csv", header=TRUE)

# Convert to tibbles:
temp_tib_exp <- as.tibble(temp_ds_exp)
temp_tib_ctrl <- as.tibble(temp_ds_ctrl)

#Add prefix to control samples:
colnames(temp_tib_ctrl) <- paste("ctrl", colnames(temp_tib_ctrl), sep=".")
#Remove prefix from Geneid col name:
names(temp_tib_ctrl)[names(temp_tib_ctrl)=="ctrl.Geneid"] <- "Geneid"

# Join the two tibbles, keeping all the rows in the sample dataset (tib1)
# but not all rows in the control dataset (tib2)
tib_comb1 <- left_join(temp_tib_exp, temp_tib_ctrl, by="Geneid")

# filter out rows that have any zeros for funzies
tib_comb2 <- subset(tib_comb1,rowSums(tib_comb1==0)==0)



# # # Deseq example shell script:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#R



# load counts
cnts <- read.csv('groups/group_EX_rna_counts.csv',row.names=1)

# filter out rows that have any zeros for funzies
cnts <- subset(cnts,rowSums(cnts==0)==0)

# sample information
info <- read.csv('groups/group_EX_rna_info.csv')

# create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = cnts,
  colData = info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','HMGCOA','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'example_deseq_results.csv')

# write out matrix of normalized counts
write.csv(counts(dds,normalized=TRUE),'example_deseq_norm_counts.csv')