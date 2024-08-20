# load the packages 
library(DESeq2)

# load the dataset - 'contig' 'length' 'Mapped bases' 'Mean coverage' 'SD' 
# the outcomes of qualimap in text file
# we changed the file names 
AG_F <- read.table("../AG_F.txt", header = F, row.name=1)
head(AG_F)

AG_M <- read.table("../AG_M.txt", header = F, row.name=1)
head(AG_M)

PG_F <- read.table("../PG_F.txt", header = F, row.name=1)
head(PG_F)

PG_M <- read.table("../PG_M.txt", header = F, row.name=1)
head(PG_M)

# making the dataset for mean coverage 
countdat <- data.frame(AG_F = c(AG_F[,3]),AG_M = c(AG_M[,3]),
                       PG_F = c(PG_F[,3]),PG_M = c(PG_M[,3]), 
                       row.names = rownames(AG_F))
head(countdat)                         
write.csv(countdat, 'mean.coverage.csv')

# make the sample file 
sample <- data.frame(ID=c('AG_F','AG_M','PG_F','PG_M'),
                     Tissue = c('AG','AG','PG','PG'))
sample$Tissue <- factor(sample$Tissue)
sample

# Create Deseq2 obj
## Round the values in countdata to the nearest integer
countdat <- round(countdat)

dds <- DESeqDataSetFromMatrix(countData=countdat, 
                              colData=sample, design=~Tissue)

# calculate fold changes 
dds <- DESeq(dds)

## see the levels 
resultsNames(dds)

## extracting the output of logfoldchanges and test 
results <- results(dds, contrast=c("Tissue", "PG", "AG"))
head(results)

# Sort results by padj in ascending order
sorted_results <- results[order(results$padj),]

# View top differentially expressed genes
head(sorted_results)

# Subset results to include only genes with log2FoldChange > 1 and padj < 0.1
sig_results <- subset(sorted_results, log2FoldChange > 1 & padj < 0.1)

# View top significant genes
head(sig_results)
dim(sig_results) # 2244 from 38822

# Write out the outcome 
write.csv(sig_results, 'DE_PGgene.csv')


