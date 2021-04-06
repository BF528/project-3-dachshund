library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_1_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb2/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical =='3-METHYLCHOLANTHRENE'|samples$chemical=='Control'])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design) <- c('Intercept','3-METHYLCHOLANTHRENE') #CLOTRIMAZOLE CHLOROFORM swap out below

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t1 <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t1,'3Methylcholanthrene_limma_results.csv')

rma.subset <- rma[paste0('X',samples$array_id[samples$chemical =='CLOTRIMAZOLE'|samples$chemical=='Control'])]

#CLOTRIMAZOLE limma results
# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','CLOTRIMAZOLE')
  )
)
colnames(design) <- c('Intercept','CLOTRIMAZOLE') 

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t2 <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t2,'Clotrimazole_limma_results.csv')

#Chloroform limma results 
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical =='CHLOROFORM'|samples$chemical=='Control'])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','CHLOROFORM')
  )
)
colnames(design) <- c('Intercept','CHLOROFORM') 

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t3 <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t3,'Chloroform_limma_results.csv')


#how many significant results
#3-METHYLCHOLANTHRENE
t1_significant <- t1[which(t1$adj.P.Val<0.05),]
#CLOTRAMIZOLE
t2_significant <- t2[which(t2$adj.P.Val<0.05),]
#CHLOROFORM
t3_significant <- t3[which(t3$adj.P.Val<0.05),]


#Create histogram of fold change values 
library(ggplot2)

qplot(t1_significant$logFC, geom="histogram",main="Histogram for FoldChange: 3-Methylcholanthrene microarray", xlab="log FC", ylab= "Frequency", fill=I("blue"),col=I("black"),xlim=c(-4,4))
qplot(t2_significant$logFC, geom="histogram",main="Histogram for FoldChange: Clotramizole microarray", xlab="log FC", ylab = "Frequency", fill=I("blue"),col=I("black"),xlim=c(-4,4))
qplot(t3_significant$logFC, geom="histogram",main="Histogram for FoldChange: Chloroform microarray", xlab="log FC", ylab = "Frequency", fill=I("blue"),col=I("black"),xlim=c(-4,4))


# Creating volcano plots with a logFC threshold on the significant DE genes to help identify the top ones in the list
#add an up regulated and down regulated section
t1_significant$regulated[t1_significant$logFC > 0.6 & -log(t1_significant$P.Value) > 18] <- "UP"
t1_significant$regulated[t1_significant$logFC < -0.6 & -log(t1_significant$P.Value) > 18] <- "DOWN"

ggplot(t1_significant, aes(x=logFC, y=-log(P.Value),col=regulated)) + geom_point() +theme_minimal() +xlim(-4,4) +ggtitle("3-Methylcholanthrene Scatter") +theme(plot.title=element_text(hjust = 0.5))


#Clotramizole
t2_significant$regulated[t2_significant$logFC > 1.2 & -log(t2_significant$P.Value) > 25] <- "UP"
t2_significant$regulated[t2_significant$logFC < -1.2 & -log(t2_significant$P.Value) > 25] <- "DOWN"
ggplot(t2_significant, aes(x=logFC, y=-log(P.Value),col=regulated)) + geom_point() +theme_minimal() +xlim(-5,5) +ggtitle("Clotramizole Scatter") +theme(plot.title=element_text(hjust = 0.5))

#Chloroform
t3_significant$regulated[t3_significant$logFC > 4 & -log(t3_significant$P.Value) > 30] <- "UP"
t3_significant$regulated[t3_significant$logFC < -4 & -log(t3_significant$P.Value) > 30] <- "DOWN"
ggplot(t3_significant, aes(x=logFC, y=-log(P.Value),col=regulated)) + geom_point() +theme_minimal() +xlim(-8,8) +ggtitle("Chloroform Scatter") +theme(plot.title=element_text(hjust = 0.5))


