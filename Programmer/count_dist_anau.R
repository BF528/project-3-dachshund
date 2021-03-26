# Creates boxplot of count distribution for samples

# Import libraries:
library(tidyverse)
library(ggplot2)

# Read in csv as tibbles:
ds1 <- as_tibble(read.csv("featureCounts_combined.csv", header=TRUE))

# Read in meta data:
meta1 <- read.csv("toxgroup_1_rna_info.csv", header=TRUE)
# Convert mode_of_action to factor variable:
meta1$mode_of_action <- factor(meta1$mode_of_action)
print(meta1[,c("Run", "mode_of_action")])

# Save current default pars:
previouspar<-par(no.readonly=TRUE)

#------------------------------------------------------------------------------

# Create boxplot:
# Create boxplot, with labels horizontal (1) and increased left margin
boxplot(ds1[,c(2:10)], horizontal=TRUE, ylim=c(0, 1000), las=1, 
        pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green")) 

# -----------------------------------------------------------------------------
# Produce boxplot on log-scale

# Add 0.0001 to each cell so log scale will work:
ds_shifted <- ds1
ds_shifted[,c(2:10)] <- ds1[,c(2:10)] + 0.0001

# Create boxplot logscale
boxplot(ds_shifted[,c(2:10)], horizontal=TRUE, las=1, log="x", 
        ylim=c(0.9, 1000000), pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green")) 




# TODO delete all below???-----------------------------------------------------
boxplot(ds1$SRR1177987)

boxplot(ds1$SRR1177987, ds1$SRR1177988, ds1$SRR1177989)
boxplot(ds1[,c("SRR1177987", "SRR1177988", "SRR1177989")])


boxplot(ds1[,c(2:10)], horizontal=TRUE, ylim=c(0, 1000), las=1) 

boxplot(ds1[,c(2:10)], horizontal=FALSE, ylim=c(0, 1000))

ggplot(data=ds1[,c(2:10)], mapping=aes(x))




# Add count of 1 to each count so can plot on log-scale:

