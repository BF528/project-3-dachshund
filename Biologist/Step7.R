#Imports
library(readr)
library(dplyr)

#set working directory for reading csv
setwd("/projectnb/bf528/users/dachshund/project_3/programmer/results")

#read for filtering and heatmapping
ahr <- read.csv("AhR_deseq_norm_counts1.csv")
cxr<- read.csv("CARPXR_deseq_norm_counts1.csv")
ctc <- read.csv("Cytotoxic_deseq_norm_counts1.csv")
top <- read.csv("top_ten_diff_exp.csv")

setwd("/projectnb/bf528/users/dachshund/project_3/biologist")


#This it to find the DE genes
ahr_de <- read.csv("AhR_diff_exp.csv")
cxr_de <- read.csv("CAR_diff_exp.csv")
ctc_de <- read.csv("CYT_diff_exp.csv")


colnames(ahr)[2:7] <- paste(colnames(ahr)[2:7], "AhR")
colnames(ctc)[2:7] <- paste(colnames(ctc)[2:7], "Cytotoxic")
colnames(cxr)[2:7] <- paste(colnames(cxr)[2:7], "CARPXR")


#Heatmapping

#funciton for cov filter
cv <- function(row){sd(row)/mean(row)}
cv.cutoff <- 0.182

#Merge dataframes for heatmap

df <- merge(ahr, ctc, by="X", all = TRUE)
df <- merge(df, cxr, by="X", all=TRUE)
rownames(df) <- df$X
df <- df[,-1]


#filtering
df.m <- median(apply(df, 1,FUN=mean, na.rm=T))
df.mean <- df[apply(df, 1,FUN=mean, na.rm=T)<df.m,]
#Higher than average variance only
#df.med <- median(apply(df.mean, 1, FUN=var, na.rm=T))
#df.var <- df.mean[apply(df.mean, 1, FUN=var, na.rm=T)>df.med,]pri

df.cov <- df[apply(df, 1, FUN=cv)<cv.cutoff,]
df.m.cov <- (apply(df.cov, 1,FUN=mean, na.rm=T))
df.mean.cov <- df.cov[apply(df.cov, 1,FUN=mean, na.rm=T)<df.m.cov[2],]

df.m.cov <- apply(df.cov, 1,FUN=mean, na.rm=T)
#Median is third
df.mean.m <- df.cov[apply(df.cov, 1,FUN=mean, na.rm=T)<df.m.cov[3],]

#Create Heatmap
heatmap(data.matrix(df.mean.m), margins=c(5,5))
