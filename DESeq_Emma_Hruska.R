#### Individual Component DESeq2 Script ####
## Group: Bonne
## Member: Emma Hruska
## Adapted from: Schwartz, T. (2025, March 6). 4_DESeq.teachingTOenrichment_2025.R [R script]. 
## GitHub repository, https://github.com/Schwartz-Lab-at-Auburn/FunGen2025/tree/main. Accessed April 10, 2025.
## Purpose: Differential gene expression analysis.

### You will need to set your working directory to the location of your data.

### Install the DESeq2 package if you have not already.
  install.packages("BiocManager")
  BiocManager::install("DESeq2")

### Load the DESeq2 library.
library(DESeq2)

#### Input Data ####

### Input the count data, gene count matrix, and labels.
countdata <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)

### Input the metadata or phenotype data.
coldata <-(read.table("pheno_data.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)

### Check that all sample IDs in coldata are also in countdata and match their orders.
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))

### Create the DESeq2 dataset and define the statistical model.
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~Tissue)

### View it.
dds

### Set factors for statistical analyses.
dds$Tissue <- factor(dds$Tissue, levels=c("Bone_marrow","Liver"))

#### Differential Expression Analysis ####
dds <- DESeq(dds)
res <- results(dds)
res

### Order results table by the smallest adjusted p-value.
resOrdered <- res[order(res$padj),]
resOrdered

### Summary statistics at p-value < 0.05.
  res05 <- results(dds, alpha=0.05)
  summary(res05)
  sum(res05$padj < 0.05, na.rm=TRUE)

#### MA Plot ####
  
  plotMA(res, main="DESeq2", ylim=c(-8,8))
  idx <- identify(res$baseMean, res$log2FoldChange)
  rownames(res)[idx]

### Write your results to a file.
  write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  
  
### Extract transformed values.
  rld <- rlog(dds)
  vsd <- varianceStabilizingTransformation(dds)
  head(assay(rld), 3)
  
### Heatmap of the count matrix.
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
  
  library("pheatmap")
  mat  <- assay(vsd)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  anno <- as.data.frame(colData(vsd)[, c("Tissue", "type")])
  df <- as.data.frame(colData(dds)[,c("Tissue","type")])
    pheatmap(mat, annotation_col = anno)
  
### Heatmap of the sample-to-sample distances.
  sampleDists <- dist(t(assay(rld)))
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$treatment)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
### Principal component plot of the samples.
  plotPCA(rld, intgroup=c("Tissue"))
  
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

### Write the transformed normalized count matrix with gene names to a tab-delimited text file that can be imported into Cytoscape.
write.table(as.data.frame(NormTransExp_withName), file="NormTransExp.txt", quote=FALSE, row.names=FALSE, sep = "\t")  
