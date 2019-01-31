
#Transcriptomics analysis Deseq2 
#After have the count table (Genes vs raw counts)

#Install DEseq2
source("https://bioconductor.org/biocLite.R")
biocLite()  
biocLite("DESeq2")

#Load the library
library('DESeq2') 

#Other libraries
library('RColorBrewer')
library (ggplot2) 
library(dplyr)


#
#ANALYSIS FOR BROCADIA
#
#Import the count table
countData = read.table(file = "~/Input files/BRO_counttable.txt", header = TRUE, sep = '\t', row.names = 1) 

#Check that Cout Table was imported properly
dim(countData)
head(countData)

#Create experiments conditions labels
colData = read.table(file = "~/Input files/BRO.colData_2conditions.txt", header = TRUE, sep = '\t', row.names = 1) 

#Create DEseq input matrix
dds <- DESeqDataSetFromMatrix(countData, colData = colData, design = ~ condition)

#Run DEseq algorith
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
head(res)

#Optional, order by adjusted p-value
res <- res[order(res$padj), ]
head(res)

#Export results of differential analysis expression in a table 
write.csv(res, file = "BRO_Expression_results_total.csv")  

#Export just genes with significant change in expression Padj <0.05
sig <- subset(res, padj < 0.05)
write.csv(sig, file = "BRO_Expression_results_significant.csv")  

#
#
#
#
#
#
#
#Other scripts
#Histogram p-values
hist(res$pvalue, breaks=50, col="grey")

#Dispersion plot 
plotDispEsts(dds)

#PCA analysis of the samples (group according expression profile between samples)
#The differential expression analysis above operates on the raw (normalized) count data. 
#But for visualizing or clustering data  you ned to work with transformed versions of the data. 
#First, use a regularlized log trasnformation while re-estimating the dispersion ignoring any information
#you have about the samples (blind=TRUE). Perform a principal components analysis and hierarchical clustering.
#rld <- rlogTransformation(dds, blind=TRUE)
#plotPCA(rld, intgroup= 'condition' ) 


#PCA with ggplot2
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1,PC2, color=condition,shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


#visualize differentially expressed genes in MA plots
plotMA(dds)


#Optional MA plot
#Could do with built-in DESeq2 function:
#I like mine better:

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}

maplot(resdata, main="MA Plot")



# Hierarchical clustering analysis let's get the actual values for the first
# few genes
head(assay(rld))
## now transpose those
t(head(assay(rld)))
## now get the sample distances from the transpose of the whole thing
dist(t(assay(rld)))
sampledist <- dist(t(assay(rld)))
plot(hclust(sampledist))


#Heatmap with hierachical clustering

# ?heatmap for help
sampledist
as.matrix(sampledist)
sampledistmat <- as.matrix(sampledist)
heatmap(sampledistmat)


#Correlation plots
plotFun <- function(x,y){ dns <- densCols(x,y); points(x,y, col=dns, pch=".") }
pairs(log2(countData  + 1), panel=plotFun, lower.panel = NULL)

#Correlation test with lower panel
#Code for lower panel
panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

#Correlation test with lower panel
plotFun <- function(x,y){ dns <- densCols(x,y); points(x,y, col=dns, pch=".") }
pairs(log2(countData  + 1), panel=plotFun, lower.panel = panel.cor)



#Visualize the reads on the genome 

# download the Integrative Genome Viewer from the Broad Institute. Download all your .bam files from your 
#AWS instance, and load them into IGV. Try navigating to regions around differentially expressed genes to view
#how reads map to genes differently in the controls versus the irradiated samples.





 
  
