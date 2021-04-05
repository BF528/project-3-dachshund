library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('\\Users\\abhis\\OneDrive\\Desktop\\group_1_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('\\Users\\abhis\\OneDrive\\Desktop\\liver-normalization-rma.txt',
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
colnames(design) <- c('Intercept','3-METHYLCHOLANTHRENE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t1 <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t1,'3-METHYLCHOLANTHRENE_limma_results.csv')

#CLOTRIMAZOLE

rma.subset <- rma[paste0('X',samples$array_id[samples$chemical =='CLOTRIMAZOLE'|samples$chemical=='Control'])]

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
write.csv(t2,'CLOTRIMAZOLE_limma_results.csv')

#Chloroform
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
write.csv(t3,'CHLOROFORM_limma_results.csv')

#how many significant results
#3-METHYLCHOLANTHRENE
t1_significant <- t1[which(t1$adj.P.Val<0.05),]
#CLOTRAMIZOLE
t2_significant <- t2[which(t2$adj.P.Val<0.05),]
#CHLOROFORM
t3_significant <- t3[which(t3$adj.P.Val<0.05),]

#Create histogram of fold change values 
library(ggplot2)
library(grid)
require(gridExtra)
plot1 <- qplot(t1_significant$logFC, geom="histogram",main="3ME foldchange", xlab="log FC", ylab= "Frequency", fill=I("blue"),col=I("black"),xlim=c(-4,4))
plot2 <- qplot(t2_significant$logFC, geom="histogram",main="CLO foldchange", xlab="log FC", ylab = "Frequency", fill=I("blue"),col=I("black"),xlim=c(-4,4))
plot3 <- qplot(t3_significant$logFC, geom="histogram",main="CHR foldchange", xlab="log FC", ylab = "Frequency", fill=I("blue"),col=I("black"),xlim=c(-4,4))
grid.arrange(plot1,plot2,plot3,nrow=1,top=textGrob("Histograms for microarray"),widths=c(3,3,3))

# Creating volcano plots with a logFC threshold on the significant DE genes to help identify the top ones in the list
#add an up regulated and down regulated section
t1_significant$regulated[t1_significant$logFC > 0.6 & -log(t1_significant$P.Value) > 18] <- "UP"
t1_significant$regulated[t1_significant$logFC < -0.6 & -log(t1_significant$P.Value) > 18] <- "DOWN"


plot4 <- ggplot(t1_significant, aes(x=logFC, y=-log(P.Value),col=regulated)) + geom_point() +theme_minimal() +xlim(-4,4) +ggtitle("3ME Scatter") +theme(plot.title=element_text(hjust = 0.5))+theme(legend.position = "none")


#Clotramizole
t2_significant$regulated[t2_significant$logFC > 1.2 & -log(t2_significant$P.Value) > 25] <- "UP"
t2_significant$regulated[t2_significant$logFC < -1.2 & -log(t2_significant$P.Value) > 25] <- "DOWN"
plot5 <- ggplot(t2_significant, aes(x=logFC, y=-log(P.Value),col=regulated)) + geom_point() +theme_minimal() +xlim(-5,5) +ggtitle("CLO Scatter") +theme(plot.title=element_text(hjust = 0.5)) +theme(legend.position = "none")


t3_significant$regulated[t3_significant$logFC > 4 & -log(t3_significant$P.Value) > 30] <- "UP"
t3_significant$regulated[t3_significant$logFC < -4 & -log(t3_significant$P.Value) > 30] <- "DOWN"
plot6 <- ggplot(t3_significant, aes(x=logFC, y=-log(P.Value),col=regulated)) + geom_point()+theme_minimal() +xlim(-8,8) +ggtitle("CHR Scatter") +theme(plot.title=element_text(hjust = 0.5))+theme(legend.position = "right")

grid.arrange(plot1, plot2, plot3,plot4,plot5,plot6,ncol=3,widths=c(1,1,1),top="Fold Change and Scatter Plots for Microarray")
#write.csv(t1_significant, file="3M_significant.csv")
#write.csv(t2_significant, file="Clot_significant.csv")
#write.csv(t3_significant, file="chlo_significant.csv")


