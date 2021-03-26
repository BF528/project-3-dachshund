# Runs deseq2_anau.R

# To run from command line:
# module load R/4.0.2
# Rscript deseq2_anau.R

# If you want to install packages switch to TRUE: 
if (FALSE){
  # Install packages:
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("DESeq2")
  BiocManager::install("apeglm")
  
  #Install hardhat:
  install.packages("hardhat", repos = "http://cran.us.r-project.org")
}

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

# Initialize list to store deseq results:
deseq_list <- list()  # To store deseq objects
deseq_tib_list  <- list()  # To store results converted to tibble:

# Run DESeq for each group: ---------------------------------------------------
for (exp_group in exp_modes_of_action){
  cat("===============================================================\n\n\n")
  print("Experimental group:")
  print(exp_group)
  
  # Get the experimental group sample names:
  temp_exp_sample_list <- meta$Run[meta$mode_of_action==exp_group]
  print("Sample experiment list:")
  print(temp_exp_sample_list)
  
  #Grab list of applicable control samples (based on vehicle value):
  exp_vehicle <- meta$vehicle[meta$Run==temp_exp_sample_list][1]
  print("exp_vehicle:")
  print(exp_vehicle)
  control_list2 <- meta$Run[meta$vehicle==exp_vehicle & meta$mode_of_action=="Control"]
  print("Control list:")
  print(control_list2)
  
  # Subset dataframe:
  temp_tib <- tib_comb_subset[,c(temp_exp_sample_list, control_list2)]
  # Names of columns remaining:
  print("Columns kept:")
  print(names(temp_tib))
  #Subset meta data:
  temp_meta <- meta[meta$mode_of_action==exp_group | (meta$mode_of_action=="Control" & meta$vehicle==exp_vehicle),]
  print("Meta data kept:")
  print(temp_meta)
  
  # Check order matches:
  print("DOES ORDER MATCH?")
  does_it_match <- length(names(temp_tib)) == sum(names(temp_tib)==temp_meta$Run)
  print(does_it_match)
  if(does_it_match){
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
  
    # Order by adjusted p-val
    res2 <- as_tibble(res1, rownames=NA)
    res2 <- res2 %>% arrange(padj)
    
    # Append results to lists:
    deseq_list[[exp_group]] <- res1
    deseq_tib_list[[exp_group]] <- res2
    
    # Strip out weird characters from exp groups:
    exp_group_str <- str_replace_all(exp_group, "[[:punct:]]", "")
    print(exp_group)
    print("Striped group name:")
    print(exp_group_str)
    
    # write out DE results
    write.csv(res2, paste(exp_group_str, 'deseq_results1.csv', sep="_"))
    
    # write out matrix of normalized counts
    write.csv(counts(dds1,normalized=TRUE),
              paste(exp_group_str, 'deseq_norm_counts1.csv', sep="_"))
    
  } else {
      print("There is an issue with sample order prior to deseq analysis")
    }
}

# Write table with top 10 results: --------------------------------------------

# Initialize list to store reformatted tibbles:
deseq_tib_list2 <- list()

for (exp_group in exp_modes_of_action){
  # Convert rownames to own column:
  deseq_tib_list2[[exp_group]] <- rownames_to_column(deseq_tib_list[[exp_group]], var="Geneid")

  # Change column names:
  colnames(deseq_tib_list2[[exp_group]]) <- paste(exp_group, colnames(deseq_tib_list2[[exp_group]]), sep="_")
} 

# Print head of each tibble in list:
print(deseq_tib_list2)

# Create CSV with top 10 results:
top_ten_tib <- bind_cols(deseq_tib_list2$AhR[c(1:10), c(1, 6)],
                          deseq_tib_list2$`CAR/PXR`[c(1:10), c(1, 6)], 
                          deseq_tib_list2$Cytotoxic[c(1:10), c(1, 6)])
write.csv(top_ten_tib, "top_ten_diff_exp.csv", row.names=FALSE)


# Determine # of significant genes for each:----------------------------------

# Initialize list:
num_sig <- list()
deseq_tib_list_sig_only <- list()

for (exp_group in exp_modes_of_action){
  deseq_tib_list_sig_only[[exp_group]] <- deseq_tib_list[[exp_group]][deseq_tib_list[[exp_group]]$padj < 0.05,]
  num_sig[[exp_group]] <- length(deseq_tib_list_sig_only[[exp_group]]$padj)
} 

# Number of significant genes for each:
num_sig

# Write to CSV:
write.csv(data.frame(num_sig),"num_sig_genes_padj.csv", row.names=FALSE)

# Create histograms for significant genes: -------------------------------------
for (exp_group in exp_modes_of_action){
  hist(deseq_tib_list_sig_only[[exp_group]]$log2FoldChange, 
       breaks=30, xlim=c(-6, 6),
       main=exp_group, xlab="log2 Fold Change")
} 

# Create scatter plots fold change vs nominal p-value for significant genes:----
for (exp_group in exp_modes_of_action){
 plot(deseq_tib_list_sig_only[[exp_group]]$log2FoldChange, 
      deseq_tib_list_sig_only[[exp_group]]$pvalue)

} 
