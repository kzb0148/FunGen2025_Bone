
####  DESeq2 Script for Differential Gene Expression Analysis in 
# Functional Genomics BIOL: 6850
### Resources and Citations:
# Love et al 2016 DESeq2 GenomeBiology
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


### You will need to set your working directory to the location you have your data.
# You can do this by using  the Session menu to set working directory To Source File Directory

#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("DESeq2")

## Load the DESeq2 library 
library(DESeq2)

## Use the Session menu to set working directory To Source File Directory

##########   1.3 Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
### How you inport this will depend on what your final output was from the mapper/counter that you used.
## this works with output from PrepDE.py from Ballgown folder.
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)

#### IF necessary,depending on what program made the count matrix. Remove the unwanted row (the * and zero row)
#countdata<- countdata[-1,c(-1)]
# OR   Remove the unwanted column (length)
#countdata<- countdata[,c(-1)]
#dim(countdata)
#head(countdata)

### Input the meta data or phenotype data
# Note: The PHENO_DATA file contains information on each sample, e.g., sex or population. The exact way to import this depends on the format of the file.
##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table("PHENO_DATA.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~Group)
#look at it
dds



#####   Prefiltering    Manual - starting at  1.3.6 
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
## You can play around with this number to see how it affects your results!
dds <- dds[ rowSums(counts(dds)) > 20, ]
# look.  How many genes were filtered out?
dds

## set factors for statistical analyses
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: Ad_lib is the control, Caloric_Restriction is the treatment group
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treated"))
dds$condition <- factor(dds$Group, levels=c("Female","Male"))


######     1.4 Differential expression analysis
### Question 2. Look at the manual - what is happening at this point?
dds <- DESeq(dds)
res <- results(dds)
res


###  Question 3. What does each column mean?
# We can order our results table by the smallest adjusted p value:
resOrdered <- res[order(res$padj),]
resOrdered
# We can summarize some basic tallies using the summary function the default is p<0.1.
summary(res)
#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)





###    1.5.1 MA-plot
##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
## Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down
plotMA(res, main="DESeq2", ylim=c(-8,8))

#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
# One can then recover the gene identiers by saving the resulting indices:
idx <- identify(res$baseMean, res$log2FoldChange)
# after selecting a gene. You need to press escape to move on
rownames(res)[idx]


##  1.5.2 Plot counts - sanity check!

# You can select the gene to plot by rowname or by numeric index.
plotCounts(dds, gene="FUN_025750", intgroup="treatment")
# You can plot the gene with th lowest adjuated P-value
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")
dds

##  Write your results to a file 
write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  

## 2.1.2 Extracting transformed values
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)

### Heatmap of the count matrix
#library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

library("pheatmap")
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group", "type")])
df <- as.data.frame(colData(dds)[,c("Group","type")])
pheatmap(mat, annotation_col = anno)

library("pheatmap")
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group", "type")])

# Adjusting font size for gene names
pheatmap(mat, 
         annotation_col = anno, 
         fontsize_row = 6,    # Smaller font size for row labels (gene names)
         fontsize_col = 10    # Font size for column labels
)



## 2.2.1 Heatmap of the count matrix
#  library("pheatmap")
#  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
#  nt <- normTransform(dds) # defaults to log2(x+1)
#  log2.norm.counts <- assay(nt)[select,]
#   df <- as.data.frame(colData(dds)[,c("treatment","type")])
#  pheatmap(assay(vsd)[mat,], cluster_rows=TRUE, show_rownames=TRUE,
#          cluster_cols=TRUE, annotation_col=df)

# pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,  cluster_cols=FALSE, annotation_col=df) 


#2.2.2 Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Group)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# 2.2.3 Principal component plot of the samples
plotPCA(rld, intgroup=c("Group"))




#### Preparing Data For GSEA and Cytoscape ####

### Import the DGE results file and make sure the gene model name is "gene_id."
DGEresults <- read.csv("DGESeq_results.csv", stringsAsFactors = FALSE)
summary(DGEresults)
dim(DGEresults)

### Rename first column so it matches "gene_id."
names(DGEresults)[1]<- "gene_id" 

#### Make Ranked List For GSEA ####

DGE_Anno_Rank <-  within(DGEresults, rank <- sign(log2FoldChange) * -log10(pvalue))
DGE_Anno_Rank 

### Subset the results with only gene name and rank.
DGErank = subset(DGE_Anno_Rank, select = c(gene_id,rank) )
DGErank

DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)

### Remove the "gene-" from row names.

gene <-gsub("^[^-]+-", "", DGErank_withName$gene_id)
DGErankIDs  <-cbind(gene,DGErank_withName)
head(DGErankIDs)
summary(DGErankIDs)  

DGErank_onlyName = subset(DGErankIDs, select = c(gene,rank) )

### Remove duplicate gene names and any genes with infinite LFC values.
DGErank_onlyName$gene <- sub("\\|.*", "", DGErank_onlyName$gene)
DGErank_onlyName <- subset(DGErank_onlyName, rank != "-Inf")
DGErank_onlyName <- subset(DGErank_onlyName, rank != "Inf")

write.table(as.data.frame(DGErank_onlyName), file="DGErankName.rnk", quote=FALSE, row.names=FALSE, sep = "\t")  

#### Normalized Expression Data ####

### Obtain the transformed normalized count matrix.
nt <- normTransform(dds)
head(assay(nt))

### Make the transformed normalized count matrix a new dataframe.
NormTransExp<-assay(nt)
summary(NormTransExp)
head(NormTransExp)

gene_id <-gsub("^[^-]+-", "", rownames(NormTransExp))
NormTransExpIDs  <-cbind(gene_id,NormTransExp)
head(NormTransExpIDs)

### Rename first column so it matches "gene_id" in annotation file.
names(NormTransExp)[1]<- "gene_id"
head(NormTransExp)

# Combine the gene IDs with the normalized expression data
NormTransExp_withName <- cbind(gene_id, NormTransExp)

# Now, write the transformed normalized count matrix to a tab-delimited text file
write.table(as.data.frame(NormTransExp_withName), file="NormTransExp.txt", quote=FALSE, row.names=FALSE, sep = "\t")




