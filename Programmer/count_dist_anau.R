# Creates boxplot of count distribution for samples

# Clean R session:
rm(list=ls())

# Import libraries:
library(tidyverse)

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
jpeg("boxplot_panel1.jpg")
boxplot(ds1[,c(2:10)], horizontal=TRUE, ylim=c(0, 1000), las=1, 
        pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green"),
        main="Distribution of Gene Counts", xlab="Gene Counts") 
legend("bottomright", legend=c("CAR/PXR", "AhR", "Cytotoxic"), 
       fill=c("green", "red", "blue"), cex=0.75)
dev.off()

# -----------------------------------------------------------------------------
# Produce boxplot on log-scale

# Add 0.0001 to each cell so log scale will work:
ds_shifted <- ds1
ds_shifted[,c(2:10)] <- ds1[,c(2:10)] + 0.0001

# Create boxplot logscale
jpeg("boxplot_panel2.jpg")
boxplot(ds_shifted[,c(2:10)], horizontal=TRUE, las=1, log="x", 
        ylim=c(0.9, 1000000), pars=list(par(mar=c(4,7,4,4))),
        col=c("blue", "blue", "blue", "red", "red", "red", "green", "green", "green"),
        main="Distribution of Gene Counts", xlab="Gene Counts (log)") 
dev.off()


# Restore par defaults:
par(previouspar)


