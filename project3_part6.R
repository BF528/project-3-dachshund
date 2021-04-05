#input files from correct directories
mapping <- read.csv('/project/bf528/project_3/refseq_affy_map.csv',as.is=TRUE)
#microarray files
mic_methylcholanthrene<- read.csv("/projectnb2/bf528/users/dachshund/project_3/analyst/3Methylcholanthrene_limma_results.csv", as.is = TRUE)
mic_clotramizole <-read.csv("/projectnb2/bf528/users/dachshund/project_3/analyst/Clotrimazole_limma_results.csv", as.is = TRUE)
mic_chloroform <- read.csv("/projectnb2/bf528/users/dachshund/project_3/analyst/Chloroform_limma_results.csv", as.is = TRUE)
#rnaseq files
rna_methylcholanthrene <-read.csv("/projectnb2/bf528/users/dachshund/project_3/programmer/results/AhR_deseq_results1.csv", as.is = TRUE)
rna_clotramizole <-read.csv("/projectnb2/bf528/users/dachshund/project_3/programmer/results/CARPXR_deseq_results1.csv", as.is = TRUE)
rna_chloroform <-read.csv("/projectnb2/bf528/users/dachshund/project_3/programmer/results/Cytotoxic_deseq_results1.csv", as.is = TRUE)


#filter steps for rna seq and microarray (THESE ARE THE MAIN DATAFRAMES TO BE WORKED ON "SIGNIFICANT RNA SEQ AND MICROARRAY DEGs"
mic_methylcholanthrene1 <- mic_methylcholanthrene[mic_methylcholanthrene$P.Value<0.05,]
rna_methylcholanthrene1 <-rna_methylcholanthrene[rna_methylcholanthrene$pvalue<0.05,]


rna_clotramizole1 <- rna_clotramizole[rna_clotramizole$pvalue<0.05,]
mic_clotramizole1 <- mic_clotramizole[mic_clotramizole$P.Value<0.05,]

rna_chloroform1 <- rna_chloroform[rna_chloroform$pvalue<0.05,]
mic_chloroform1 <- mic_chloroform[mic_chloroform$P.Value<0.05,]

#probe column for microarray 
mic_methyl_probes <- mic_methylcholanthrene1$X
mic_clotr_probes <- mic_clotramizole1$X
mic_chloro_probes <- mic_chloroform1$X

#matching limma data to mapping reference
matched_limma_methylcholanthrene <- mapping$PROBEID %in% mic_methyl_probes
matched_limma_clotramizole <- mapping$PROBEID %in% mic_clotr_probes
matched_limma_chloroform <- mapping$PROBEID %in% mic_chloro_probes

#matching rna-seq data to mapping reference
matched_deseq_methylcholanthrene <- mapping$REFSEQ %in% rna_methylcholanthrene1$X
matched_deseq_clotramizole <- mapping$REFSEQ %in% rna_clotramizole1$X
matched_deseq_chloroform <- mapping$REFSEQ %in% rna_chloroform1$X

#3-methylcholanthrene index 
mic_methylcholanthrene_index <- mic_methylcholanthrene1$logFC[match(mapping$PROBEID,mic_methylcholanthrene1$X)]
rna_methylcholanthrene_index <- rna_methylcholanthrene1$log2FoldChange[match(mapping$REFSEQ,rna_methylcholanthrene1$X)]
same_FC_direction <- sign(mic_methylcholanthrene_index) == sign(rna_methylcholanthrene_index)

#Clotramizole
mic_clotramizole_index <- mic_clotramizole1$logFC[match(mapping$PROBEID,mic_clotramizole1$X)]
rna_clotramizole_index <- rna_clotramizole1$log2FoldChange[match(mapping$REFSEQ,rna_clotramizole1$X)]
same_FC_direction2 <- sign(mic_clotramizole_index) == sign(rna_clotramizole_index)

#Chloroform
mic_chloroform_index <- mic_chloroform1$logFC[match(mapping$PROBEID,mic_chloroform1$X)]
rna_chloroform_index <- rna_chloroform1$log2FoldChange[match(mapping$REFSEQ,rna_chloroform1$X)]
same_FC_direction3 <- sign(mic_chloroform_index) == sign(rna_chloroform_index)




#matches between microarray, RNA-SEQ, and FC direction by positive or negative sign
#3-methylcholanthrene
intersection1 <- matched_limma_methylcholanthrene & matched_deseq_methylcholanthrene & same_FC_direction
#clotramizole
intersection2 <- matched_limma_clotramizole & matched_deseq_clotramizole & same_FC_direction2
#chloroform
intersection3 <- matched_limma_chloroform & matched_deseq_chloroform & same_FC_direction3

#how many observed interactions between same fold change and RNA-seq and microarray
intersection1_sum <- sum(intersection1)
intersection2_sum <- sum(intersection2)
intersection3_sum <- sum(intersection3)
#total number of matched differentially expressed genes using the mapping reference for microarray and RNA-seq)
total1 <- sum(matched_limma_methylcholanthrene | matched_deseq_methylcholanthrene)
total2 <- sum(matched_limma_clotramizole | matched_deseq_clotramizole)
total3 <- sum(matched_limma_chloroform | matched_deseq_chloroform)
#Using cross platform concordance analysis equation
concordance_3methylcholanthrene <- 2*intersection1_sum/total1
concordance_clotramizole <- 2*intersection2_sum/total2
concordance_chloroform <- 2*intersection3_sum/total3
#set n0 as intersection1_sum
n0 <- intersection1_sum
n0_2 <- intersection2_sum
n0_3 <- intersection3_sum
#number of items in independent set 1 (microarray for 3-methylcholanthrene)
n1 <- dim(rna_methylcholanthrene1)[1]
n1_2 <- dim(rna_clotramizole1)[1]
n1_3 <- dim(rna_chloroform1)[1]
#number of items in independent set 2 (RNA-seq for 3-methylcholanthrene)
n2 <- dim(mic_methylcholanthrene1)[1]
n2_2 <- dim(mic_clotramizole1)[1]
n2_3 <- dim(mic_chloroform1)[1]
#dimension of mapping reference
N <- 25225
#background corrected intersection 
new_intersection1 <- (n0*N-n1*n2)/(n0+N-n1-n2)
new_intersection2 <- (n0_2*N-n1_2*n2_2)/(n0_2+N-n1_2-n2_2)
new_intersection3 <- (n0_3*N-n1_3*n2_3)/(n0_3+N-n1_3-n2_3)
#corrected concordance
corrected_concordance_3methylcholanthrene <- 2*new_intersection1/total1
corrected_concordance_clotramizole <- 2*new_intersection2/total2
corrected_concordance_chloroform <- 2*new_intersection3/total3

#datapoints for MOAs
rna_seq_DEG <- c(dim(rna_methylcholanthrene1)[1],dim(rna_clotramizole1)[1], dim(rna_chloroform1)[1])
mic_DEG <- c(dim(mic_methylcholanthrene1)[1],dim(mic_clotramizole1)[1],dim(mic_chloroform1)[1])
concordance <- c(corrected_concordance_3methylcholanthrene,corrected_concordance_clotramizole,corrected_concordance_chloroform)
name <- c("3ME","CLO", "CHR")
plot(x=rna_seq_DEG, y=concordance*100,type="p",pch=19, col="red", main="Concordance plot between microarray and RNA-Seq",xlab="Treatment effect (number of DEGs from RNA-Seq)",ylab="Concordance of DEG (%)",ylim=c(0,70), xlim=c(1000,3500))
text(x=rna_seq_DEG,y=concordance*100,labels=name, pos=3)

plot(x=mic_DEG, y=concordance*100, type="p", pch=20, col="blue", main="Concordance plot between microarray and RNA-Seq",xlab="Treatment effect (number of DEGs from Microarray)",ylab="Concordance of DEG (%)", ylim=c(0,70),xlim=c(2000,15000))
text(x=mic_DEG,y=concordance*100,labels=name, pos=3)


####
#compute medians of deseq and microarray significant DEGs data
#remove X.1 column (Filled with NAs and na.omit is screwing it up
rna_methylcholanthrene1 <- rna_methylcholanthrene1[1:6]
rna_clotramizole1 <- rna_clotramizole1[1:6]
rna_chloroform1 <- rna_chloroform1[1:6]
mic_methylcholanthrene1 <- mic_methylcholanthrene1[1:6]
mic_clotramizole1 <- mic_methylcholanthrene1[1:6]
mic_chloroform1 <- mic_chloroform1[1:6]

#remove any NA rows
rna_methylcholanthrene1 <- na.omit(rna_methylcholanthrene1)
mic_methylcholanthrene1 <- na.omit(mic_methylcholanthrene1)
rna_clotramizole1 <- na.omit(rna_clotramizole1)
mic_clotramizole1 <- na.omit(mic_clotramizole1)
rna_chloroform1 <- na.omit(rna_chloroform1)
mic_chloroform1 <- na.omit(mic_chloroform1)

#filter significant DEGs by a median in both RNA seq and microarray data
med_0 <- median(rna_methylcholanthrene1$baseMean)
rnaseq_above_median_methyl <- rna_methylcholanthrene1[rna_methylcholanthrene1$baseMean>=med_0,]
rnaseq_below_median_methyl <- rna_methylcholanthrene1[rna_methylcholanthrene1$baseMean<med_0,]

med_1 <- median(mic_methylcholanthrene1$AveExpr)
mic_above_median_methyl <- mic_methylcholanthrene1[mic_methylcholanthrene1$AveExp>=med_1,]
mic_below_median_methyl <- mic_methylcholanthrene1[mic_methylcholanthrene1$AveExpr<med_1,]

med_2 <- median(rna_clotramizole1$baseMean)
rnaseq_above_median_clot <- rna_clotramizole1[rna_clotramizole1$baseMean>=med_2,]
rnaseq_below_median_clot <- rna_clotramizole1[rna_clotramizole1$baseMean<med_2,]

med_3 <- median(mic_clotramizole1$AveExpr)
mic_above_median_clot <- mic_clotramizole1[mic_clotramizole1$AveExpr>=med_3,]
mic_below_median_clot <- mic_clotramizole1[mic_clotramizole1$AveExpr<med_3,]

med_4 <- median(rna_chloroform1$baseMean)
rnaseq_above_median_chlo <- rna_chloroform1[rna_chloroform1$baseMean>=med_4,]
rnaseq_below_median_chlo <- rna_chloroform1[rna_chloroform1$baseMean<med_4,]

med_5 <- median(mic_chloroform1$AveExpr)
mic_above_median_chlo <- mic_chloroform1[mic_chloroform1$AveExpr>=med_5,]
mic_below_median_chlo <- mic_chloroform1[mic_chloroform1$AveExpr<med_5,]


#merge rnaseq and microarray using the mapping reference, you can find matches between above and below data frames
# and remove duplicates as it may affect the results with multiple probes for each gene? 
match_methyl_above <- mapping[which(mapping$PROBEID %in% mic_above_median_methyl$X & mapping$REFSEQ %in% rnaseq_above_median_methyl$X),]
match_methyl_above <-match_methyl_above[!duplicated(match_methyl_above$SYMBOL),]

match_methyl_below <- mapping[which(mapping$PROBEID %in% mic_below_median_methyl$X & mapping$REFSEQ %in% rnaseq_below_median_methyl$X),]
match_methyl_below <- match_methyl_below[!duplicated(match_methyl_below$SYMBOL),]

match_clot_above <- mapping[which(mapping$PROBEID %in% mic_above_median_clot$X & mapping$REFSEQ %in% rnaseq_above_median_clot$X),]
match_clot_above <- match_clot_above[!duplicated(match_clot_above$SYMBOL),]

match_clot_below <- mapping[which(mapping$PROBEID %in% mic_below_median_clot$X & mapping$REFSEQ %in% rnaseq_below_median_clot$X),]
match_clot_below <-match_clot_below[!duplicated(match_clot_below$SYMBOL),]

match_chlo_above <- mapping[which(mapping$PROBEID %in% mic_above_median_chlo$X & mapping$REFSEQ %in% rnaseq_above_median_chlo$X),]
match_chlo_above <- match_chlo_above[!duplicated(match_chlo_above$SYMBOL),]

match_chlo_below<- mapping[which(mapping$PROBEID %in% mic_below_median_chlo$X & mapping$REFSEQ %in% rnaseq_below_median_chlo$X),]
match_chlo_below <- match_chlo_below[!duplicated(match_chlo_below$SYMBOL),]

# you have the matched cases in above and below median per chemical 
# now we have to find the individual sets for microarray and rnaseq for both above and below sets per chemical (n1 & n2)
# remove duplicates as was done before (I wonder what the effect would be if duplicates weren't removed, but I'm too lazy to figure that out
#the "a" stands for "above" and the "b" stands for "below"
#only use the corrected concordances for the above and below cases

#3methylcholanthrene
n_1a <- mapping[which(mapping$PROBEID %in% mic_above_median_methyl$X),]  
n_1a <- n_1a[!duplicated(n_1a$SYMBOL),]
n_1a <- nrow(n_1a)
n_2a <- mapping[which(mapping$REFSEQ %in% rnaseq_above_median_methyl$X),]  
n_2a <- n_2a[!duplicated(n_2a$SYMBOL),]
n_2a <- nrow(n_2a)
concordance_3me_above <- 2*nrow(match_methyl_above)/(n_1a+n_2a)
above_x_methyl <- (nrow(match_methyl_above)*N-n_1a*n_2a)/(nrow(match_methyl_above)+N-n_1a-n_2a)
corrected_concordance_3me_above <- 2*above_x_methyl/(n_1a +n_2a)
#below
n_1b <-mapping[which(mapping$PROBEID %in% mic_below_median_methyl$X),]
n_1b <-n_1b[!duplicated(n_1b$SYMBOL),]
n_1b <-nrow(n_1b)
n_2b <-mapping[which(mapping$REFSEQ %in% rnaseq_below_median_methyl$X),]
n_2b <-n_2b[!duplicated(n_2b$SYMBOL),]
n_2b <-nrow(n_2b)
concordance_3me_below <-2*nrow(match_methyl_below)/(n_1b+n_2b)
below_x_methyl <- (nrow(match_methyl_below)*N-n_1b*n_2b)/(nrow(match_methyl_below)+N-n_1b-n_2b)
corrected_concordance_3me_below <- 2*below_x_methyl/(n_1b +n_2b)
#clotramizaole
n_3a <- mapping[which(mapping$PROBEID %in% mic_above_median_clot$X),]  
n_3a <- n_3a[!duplicated(n_3a$SYMBOL),]
n_3a <- nrow(n_3a)
n_4a <- mapping[which(mapping$REFSEQ %in% rnaseq_above_median_clot$X),]  
n_4a <- n_4a[!duplicated(n_4a$SYMBOL),]
n_4a <- nrow(n_4a)
concordance_clot_above <- 2*nrow(match_clot_above)/(n_3a+n_4a) 
above_x_clot <- (nrow(match_clot_above)*N-n_3a*n_4a)/(nrow(match_clot_above)+N-n_3a-n_4a)
corrected_concordance_clot_above <- 2*above_x_clot/(n_3a +n_4a)
#below
n_3b <-mapping[which(mapping$PROBEID %in% mic_below_median_clot$X),]
n_3b <-n_3b[!duplicated(n_3b$SYMBOL),]
n_3b <-nrow(n_3b)
n_4b <-mapping[which(mapping$REFSEQ %in% rnaseq_below_median_clot$X),]
n_4b <-n_4b[!duplicated(n_4b$SYMBOL),]
n_4b <-nrow(n_4b)
concordance_clot_below <-2*nrow(match_clot_below)/(n_3b+n_4b)
below_x_clot <- (nrow(match_clot_below)*N-n_3b*n_4b)/(nrow(match_clot_below)+N-n_3b-n_4b)
corrected_concordance_clot_below <- 2*below_x_clot/(n_3b +n_4b)
#chloroform
n_5a <- mapping[which(mapping$PROBEID %in% mic_above_median_chlo$X),]  
n_5a <- n_5a[!duplicated(n_5a$SYMBOL),]
n_5a <- nrow(n_5a)
n_6a <- mapping[which(mapping$REFSEQ %in% rnaseq_above_median_chlo$X),]  
n_6a <- n_6a[!duplicated(n_6a$SYMBOL),]
n_6a <- nrow(n_6a)
concordance_chlo_above <- 2*nrow(match_chlo_above)/(n_5a+n_6a)
above_x_chlo <- (nrow(match_chlo_above)*N-n_5a*n_6a)/(nrow(match_chlo_above)+N-n_5a-n_6a)
corrected_concordance_chlo_above <- 2*above_x_chlo/(n_5a +n_6a)
#below
n_5b <-mapping[which(mapping$PROBEID %in% mic_below_median_chlo$X),]
n_5b <-n_5b[!duplicated(n_5b$SYMBOL),]
n_5b <-nrow(n_5b)
n_6b <-mapping[which(mapping$REFSEQ %in% rnaseq_below_median_chlo$X),]
n_6b <-n_6b[!duplicated(n_6b$SYMBOL),]
n_6b <-nrow(n_6b)
concordance_chlo_below <-2*nrow(match_chlo_below)/(n_5b+n_6b)
below_x_chlo <- (nrow(match_chlo_below)*N-n_5b*n_6b)/(nrow(match_chlo_below)+N-n_5b-n_6b)
corrected_concordance_chlo_below <- 2*below_x_chlo/(n_5b +n_6b)




#combine the chemical concordances for all, above, and below median and format into barplot
all_conc_methyl <- 100*c(corrected_concordance_3me_above, corrected_concordance_3methylcholanthrene, corrected_concordance_3me_below)
all_conc_clot <- 100*c(corrected_concordance_clot_above, corrected_concordance_clotramizole, corrected_concordance_clot_below)
all_conc_chloro <- 100*c(corrected_concordance_chlo_above,corrected_concordance_chloroform,corrected_concordance_chlo_below)

bar <- as.data.frame(rbind(all_conc_methyl,all_conc_clot,all_conc_chloro))
row.names(bar) <- c("AhR: 3ME", "CAR/PXR: CLO", "Cytotoxic: CHR")
colnames(bar) <- c("Above median", "All", "Below median")

bar <- t(bar)

barplot(bar, main="Concordance for 3 MOAs", xlab="MOA", col=c("darkblue", "red", "white"), legend=TRUE, args.legend = list(x="top",cex=0.7),ylim=c(0,60),beside=TRUE, ylab="Concordance (%)")
