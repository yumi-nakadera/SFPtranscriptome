#load the packages (DEseq2, apeglm, edgeR are from bioconductor) 

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library('tools')
library("apeglm")
library("dplyr")
library('edgeR') # for cpm
library(stringr) #needed on M1 
library('ashr')
library(VennDiagram)
library(ggplotify)
library(ggpubr)

# Set some directories for input/output files
count_dir <- "../Data/quant"
list.files(count_dir)

# Read in count data
# create a empty dataframe called co to merge the data into
co <- data.frame()

# using for loop read all the count files in the count_dir path
# iterate over each "quant.sf" file
# merge each file, i, into a single data frame

count_files <- list.files(path = count_dir,recursive=TRUE,pattern="quant.sf",full.names=TRUE)
count_files

for (i in count_files) {
  
  # print the file that is being loaded
  print(paste0("reading file: ", i))
  
  # read file i as a data frame. note the reads info is in 5th column (Salmon output)
  f <- read.table(i, sep = "\t", header = TRUE)[,c(1,5)]
  
  # extract the sample ID from the file name:
  sam <- str_extract(i,regex("Ls_.."))
  
  # rename the columns
  colnames(f) <- c("gene_id", sam)
  
  #if the counts object is empty just copy the f to m
  if(length(co) == 0){
    co <- f
    
  } else 
  {
    #if the dataframe is not empty then merge the data
    co <- merge(co, f, by.x = "gene_id", by.y = "gene_id")
  }
  rm(f)
}

#grab the rows from the 1st column and use it as the row-names in the dataframe
rownames(co) <- co[,1]
co <- co[,-1]

# create a vector containing transcript lengths
# read in length
df_length <- read.table(count_files[1], sep = "\t", header = TRUE)[,1:2]

# create named vector
mylength <- df_length$Length
names(mylength) <- df_length[,1]

#get metadata
sample = read.table("../Data/quant/sample.txt", header = T)
sample

#create data frame with each count file treatment
sampletable = data.frame(
  sampleName = sample$Sample,
  fileName = count_files,
  Treatment = sample$Treatment,
  sample = sample$sample)

sampletable

# Preliminary data exploration and filtering
# how many transcripts did we quantify expression for?
dim(co) #38,820 genes (2 missing somehow?)

# what is the distribution of total expression across samples?
rowSums(co) %>% log(.,10) %>% hist(., breaks = 100, main ='', ylab='',xlab='',
                                   xlim=c(-2,8), ylim=c(0,2500),axes = F)

# Before moving forward, we're going to exclude genes that are most likely to be useless in the analysis
# let's throw out the transcripts with genes with an expression lower than 1 counts per million (CPM)
keep <- rowSums(cpm(co) >= 1) >= 1
table(keep) # TRUE= 18653

co2 <- co[keep,]
mylength <- mylength[keep]

par(new=T)
rowSums(co2) %>% log(.,10) %>% hist(., breaks = 50, main ='', 
                                    xlab='ln(Total expression across samples per gene)', 
                                    xlim=c(-2,8), ylim=c(0,2500),col='orange')

##############################################################

#create DESeq data set for analysis
ddsHTSeq <- DESeqDataSetFromMatrix(
  countData = round(co2),
  design = ~ Treatment,
  colData = sampletable
)
#check the Treatment
ddsHTSeq$Treatment

#to see the counts
counts(ddsHTSeq) %>% head()

# To replace the order with one of your choosing, create a vector with the order you want:
treatments <- c('Isolated', 'Paired','Grouped')

# Then reset the factor levels:
ddsHTSeq$Treatment <- factor(ddsHTSeq$Treatment, levels = treatments)

# verify the order
ddsHTSeq$Treatment

######################################################
# Run the statistical analysis
dds <- DESeq(ddsHTSeq)

#dispersion plot 
plotDispEsts(dds, ylab='Dispersion', xlab= 'Mean of normalized counts')

# get results table
resDE <- results(dds)

# get a quick summary of the table
summary(resDE)

#out of 18653 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 60, 0.32%
#LFC < 0 (down)     : 25, 0.13%
#outliers [1]       : 78, 0.42%
#low counts [2]     : 1052, 5.6%
#(mean count < 4)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results


# check out the first few lines
head(resDE, 10)
######################################################
# Get a table of shrunken log2 fold changes
######################################################

# see coefficient names:
resultsNames(dds)

# get shrunken log fold changes, specifying the coefficient to I vs P 
res_shrink <- lfcShrink(dds, coef="Treatment_Paired_vs_Isolated", type="ashr") 

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=resDE$log2FoldChange,
  y=res_shrink$log2FoldChange,pch=20,
  cex=.2,
  col=1+(resDE$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change",
  xlim=c(-5,5),
  ylim=c(-5,5)
)
abline(0,1)

points(
  x=resDE$log2FoldChange,
  y=res_shrink$log2FoldChange,pch=20,
  cex=.2,
  col=4 * is.na(resDE$padj)
)

# get the top 20 genes by shrunken log2 fold change
top50 <- order(-abs(res_shrink$log2FoldChange))[1:50]
res_shrink[top50,]

######################################################
# Data visualization
######################################################

# MA plot, specify it is DEseq2 one .. 
par(mfrow=c(1,2))
DESeq2::plotMA(resDE, ylim=c(-4,4), cex=0.7, 
               xlab='Mean of normalized count',
               ylab='Log Fold Change')

DESeq2::plotMA(res_shrink, ylim=c(-4,4), cex=0.7, 
               xlab='Mean of normalized count',
               ylab='Shrunk log Fold Change')


# distribution of log2 fold changes:
# there should be a peak at 0
par(mfrow=c(1,1))
hist(res_shrink$log2FoldChange,breaks=200, xlim=c(-1,1))
abline(v=0,col="red",lwd=2) 

##############
#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=1,
     xlim=c(-10,10),
     col=(log_padj > 1)+1, # color padj < 0.1 red
     ylab="Negative log-scaled adjusted p-value",
     xlab="Shrunken log2 fold changes")

#############
# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="Treatment")

# alternatively, using ggplot
dat <- plotPCA(vsd, intgroup="Treatment",returnData=TRUE)
percentVar <- round(100 * attr(dat, "percentVar"))

## pca plot with sample names 
p <- ggplot(dat,aes(x=PC1,y=PC2,col=Treatment)) +
  geom_point() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel(aes(label=name)) + theme_test() + theme(legend.title=element_blank()) 
p 

## pca plot without sample names 
colors <- c('Isolated' = '#00BA38', 'Paired' = '#F8766D', 'Grouped' = '#619CFF')

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group)) +
  geom_point(size =5) + 
  scale_color_manual(values = colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_test() + theme(legend.title=element_blank(),text = element_text(size = 12)) 
p

##############
# heatmap of DE genes

# get shrunken log fold changes, specifying the coefficient to I vs P 
res_shrink.ip <- lfcShrink(dds, coef="Treatment_Paired_vs_Isolated", type="ashr") #some warnings 

# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# get top 50 log fold change genes (excluding cook's cutoff outliers)
## I vs P
top50 <- data.frame(res_shrink.ip) 
top50 <- filter(top50, !is.na(top50$padj))
top50 <- arrange(top50, -abs(log2FoldChange))
top50 <- rownames(top50)
top50 <- head(top50,n=50)

df <- data.frame(colData(dds)[,"Treatment"])
rownames(df) <- colnames(dds)
colnames(df) <- "Treatment"

sampleorder <- c('Ls_01','Ls_02','Ls_06','Ls_12','Ls_13',
                 'Ls_03','Ls_04','Ls_07','Ls_08','Ls_14',
                 'Ls_05','Ls_09','Ls_10','Ls_11','Ls_15')
ordered.rld <- rld[,sampleorder]

p1 <- pheatmap(
  assay(ordered.rld)[top50,], 
  cluster_rows=T, 
  show_rownames=T,
  cluster_cols=F,
  annotation_col=df, 
  show_colnames = F,
  fontsize_row	= 7,
  main = 'Top 50 from Isolated vs Paired contrast', 
)
p1

## I vs G
# see coefficient names:
resultsNames(dds)

# get shrunken log fold changes, specifying the coefficient to i vd G
res_shrink.ig <- lfcShrink(dds, coef="Treatment_Grouped_vs_Isolated", type="ashr") #some warnings 

top50 <- data.frame(res_shrink.ig) 
top50 <- filter(top50, !is.na(top50$padj))
top50 <- arrange(top50, -abs(log2FoldChange))
top50 <- rownames(top50)
top50 <- head(top50,n=50)

df <- data.frame(colData(dds)[,"Treatment"])
rownames(df) <- colnames(dds)
colnames(df) <- "Treatment"

sampleorder <- c('Ls_01','Ls_02','Ls_06','Ls_12','Ls_13',
                 'Ls_03','Ls_04','Ls_07','Ls_08','Ls_14',
                 'Ls_05','Ls_09','Ls_10','Ls_11','Ls_15')
ordered.rld <- rld[,sampleorder]

p2 <-pheatmap(
  assay(ordered.rld)[top50,], 
  cluster_rows=T, 
  show_rownames=T,
  cluster_cols=F,
  annotation_col=df,
  show_colnames = F,
  fontsize_row	= 7,
  main = 'Top 50 from Isolated vs Grouped contrast'
)
p2

p1 <-  as.ggplot(p1)
p2 <-  as.ggplot(p2)

ggarrange(p1, p2, common.legend = F )

####################################
# Filtering for putative SFP genes #
####################################
summary(res_shrink.ip)
summary(res_shrink.ig)

# get significant log fold change genes (excluding cook's cutoff outliers)
# alpha = 0.1, I vs P 
sfp.ip <- data.frame(res_shrink.ip) 
sfp.ip <- filter(sfp.ip, !is.na(sfp.ip$padj)) #
sfp.ip <- filter(sfp.ip, log2FoldChange > 0) # excluding the case that I > P 
sfp.ip <- arrange(sfp.ip, -abs(log2FoldChange))
head(sfp.ip)

thr <- sum(sfp.ip$padj < 0.1, na.rm=TRUE) #138
sfp.ip <- rownames(sfp.ip)
sfp.ip <- head(sfp.ip,n=thr)
length(sfp.ip)

## I vs G 
sfp.ig <- data.frame(res_shrink.ig) 
sfp.ig <- filter(sfp.ig, !is.na(sfp.ig$padj))
sfp.ig <- filter(sfp.ig, log2FoldChange > 0) # excluding the case that I > G
sfp.ig <- arrange(sfp.ig, -abs(log2FoldChange))
head(sfp.ig)

thr <- sum(sfp.ig$padj < 0.1, na.rm=TRUE) # 57
sfp.ig <- rownames(sfp.ig)
sfp.ig <- head(sfp.ig,n=thr) 

length(sfp.ig)

## making the list of gene names with three category (A,B,C)
test <- sfp.ip %in% sfp.ig # which genes are overlapped? 
table(test) # f = 95, t = 46 

str(test)
test <- as.numeric(test) # convert logi to numbers
test <- replace(test, test==1, 'A') # overlapped, I < P = G
test <- replace(test, test==0, 'C') # Gray, I < P > G 

test2 <- sfp.ig %in% sfp.ip
table(test2) # f = 14, t = 46 
test2 <- as.numeric(test2)
test2 <- replace(test2, test2 == 0, 'B') 
dat <- data.frame(sfp.ig, test2)
head(dat)
dat <- subset(dat, test2 == 'B') # 8 genes for senario B, I = P > G
dat

sfp.gene <- c(sfp.ip, dat$sfp.ig) # combine the vectors to make the data frame 
Pattern <- c(test, dat$test2)
sfp.gene <- data.frame(sfp.gene, Pattern)
sfp.gene$Pattern <- as.factor(sfp.gene$Pattern)

summary(sfp.gene) # n=150

## test if these candidates are PG specific
pg <- read.csv('../Data/DE_PGgene.csv', header = T)

summary(pg)
head(pg)

test <- sfp.gene$sfp.gene %in% pg$X # counting how many candidate SFP genes are PG-specific expressed gene
table(test) 

# FALSE  TRUE 
# 73    77 

sfp.gene2 <- sfp.gene[test,]
summary(sfp.gene2)


## visualization for candidate SFP genes 
## venn diagram to show how many genes are overlapped
cbPalette <- c("#999999", "#E69F00")

summary(sfp.gene2)

test <- sfp.ip %in% pg$X
table(test)
sfp.ip2 <- sfp.ip[test]
summary(sfp.ip2)

test <- sfp.ig %in% pg$X
table(test)
sfp.ig2 <- sfp.ig[test]
summary(sfp.ig2)

plt <- venn.diagram(x= list(sfp.ip2, sfp.ig2), category.names =c('IvsP','IvsG'),
                    fill = cbPalette,  filename = NULL, cex = 2, cat.cex = 2, lwd = 2,ext.text = F)
grid.newpage()
grid::grid.draw(plt)


## making the plot for each gene in each scenario 
PatternA <- subset(sfp.gene2, Pattern == 'A')
PatternB <- subset(sfp.gene2, Pattern == 'B')
PatternC <- subset(sfp.gene2, Pattern == 'C')

# sc A
gene_set <- PatternA$sfp.gene
names(gene_set) <- gene_set

df <- lapply(gene_set, \(x) {
  y <- plotCounts(dds, x, c("Treatment"), returnData=TRUE)
  y$feature <- x
  return(y)
})

df <- do.call(rbind, df)
head(df)
dim(df)

ggplot(df, aes(x=Treatment, y=count)) +
  geom_jitter(width = 0.05)+ facet_wrap(~feature, scales = "free", ncol = 6) + 
  labs(subtitle ='Scenario A') + theme_test()

# sc B
gene_set <- PatternB$sfp.gene
names(gene_set) <- gene_set

df <- lapply(gene_set, \(x) {
  y <- plotCounts(dds, x, c("Treatment"), returnData=TRUE)
  y$feature <- x
  return(y)
})

df <- do.call(rbind, df)
head(df)
dim(df)

ggplot(df, aes(x=Treatment, y=count)) +
  geom_jitter(width = 0.05)+ facet_wrap(~feature, scales = "free", ncol = 4) + 
  labs(subtitle ='Scenario B') + theme_test()

# sc C
gene_set <- PatternC$sfp.gene
names(gene_set) <- gene_set

df <- lapply(gene_set, \(x) {
  y <- plotCounts(dds, x, c("Treatment"), returnData=TRUE)
  y$feature <- x
  return(y)
})

df <- do.call(rbind, df)
head(df)
dim(df)

ggplot(df, aes(x=Treatment, y=count)) +
  geom_jitter(width = 0.05)+ facet_wrap(~feature, scales = "free", ncol = 8) + 
  labs(subtitle ='Scenario C') + theme_test()


## write out the file for candidate SFP gene names 
summary(sfp.gene2)
write.csv(sfp.gene2, file = '../Data/candidate.sfp.gene.name.csv')

## test if these list includes known SFP genes 
ly5 <- '03_TRINITY_DN8105_c1_g1_i1.p1'
plotCounts(dds, ly5, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly5, c("Treatment"))
ly5 %in% sfp.gene2$sfp.gene #F
ly5 %in% sfp.gene$sfp.gene #F
ly5 %in% rownames(res_shrink.ip) #T
ly5 %in% rownames(res_shrink.ig) #T
ly5 %in% sfp.ip
ly5 %in% sfp.ig


ly7a1 <- '08_TRINITY_DN216_c0_g1_i7.p1'
plotCounts(dds, ly7a1, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly7a1, c("Treatment"))
ly7a1 %in% sfp.gene2$sfp.gene #F
ly7a1 %in% sfp.gene$sfp.gene #F
ly7a1 %in% rownames(res_shrink.ip) #T
ly7a1 %in% rownames(res_shrink.ig) #T
ly7a1 %in% sfp.ip
ly7a1 %in% sfp.ig


ly7a2 <- '11_TRINITY_DN343_c0_g1_i2.p2'
plotCounts(dds, ly7a2, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly7a2, c("Treatment"))
ly7a2 %in% sfp.gene2$sfp.gene #F
ly7a2 %in% sfp.gene$sfp.gene #F
ly7a2 %in% rownames(res_shrink.ip) #T
ly7a2 %in% rownames(res_shrink.ig) #T
ly7a2 %in% sfp.ip
ly7a2 %in% sfp.ig

ly7b <- '02_TRINITY_DN508_c0_g4_i2.p1'
plotCounts(dds, ly7b, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly7b, c("Treatment"))
ly7b %in% sfp.gene2$sfp.gene #F
ly7b %in% sfp.gene$sfp.gene #F
ly7b %in% rownames(res_shrink.ip) #T
ly7b %in% rownames(res_shrink.ig) #T
ly7b %in% sfp.ip
ly7b %in% sfp.ig

ly8a <- '01_TRINITY_DN4850_c1_g1_i1.p1'
plotCounts(dds, ly8a, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly8a, c("Treatment"))
ly8a %in% sfp.gene2$sfp.gene #F
ly8a %in% sfp.gene$sfp.gene #F
ly8a %in% rownames(res_shrink.ip) #T
ly8a %in% rownames(res_shrink.ig) #T
ly8a %in% sfp.ip
ly8a %in% sfp.ig

ly8b <- '01_TRINITY_DN271_c0_g1_i5.p1'
plotCounts(dds, ly8b, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly8b, c("Treatment"))
ly8b %in% sfp.gene2$sfp.gene #T
ly8b %in% sfp.gene$sfp.gene #T
ly8b %in% rownames(res_shrink.ip) #T
ly8b %in% rownames(res_shrink.ig) #T
ly8b %in% sfp.ip
ly8b %in% sfp.ig

ly10 <- '02_TRINITY_DN114_c0_g1_i3.p1'
plotCounts(dds, ly10, c("Treatment"), returnData=TRUE)
plotCounts(dds, ly10, c("Treatment"))
ly10 %in% sfp.gene2$sfp.gene #F
ly10 %in% sfp.gene$sfp.gene #F
ly10 %in% rownames(res_shrink.ip) #T
ly10 %in% rownames(res_shrink.ig) #T
ly10 %in% sfp.ip
ly10 %in% sfp.ig


### only ly8b is present in the candidate SFP gene list 
### outcome did not change if I use alpha = 0.05 or 0.1
### I think it is due to the thresholds, because all the pattern of LFC look not that bad (except LyAcp5) 

par(mfrow=c(3,2))
plotCounts(dds, ly5, c("Treatment"), sub = ly5, main = 'LyAcp5', xlab='', pch =16) 
plotCounts(dds, ly7a1, c("Treatment"), sub = ly5, main = 'LyAcp7a1', xlab='', pch =16)
plotCounts(dds, ly7b, c("Treatment"),sub = ly5, main = 'LyAcp7b', xlab='', pch =16)
plotCounts(dds, ly8a, c("Treatment"),sub = ly5, main = 'LyAcp8a', xlab='', pch =16)
plotCounts(dds, ly8b, c("Treatment"),sub = ly5, main = 'LyAcp8b', xlab='', pch =16)
plotCounts(dds, ly10, c("Treatment"),sub = ly5, main = 'LyAcp10', xlab='', pch =16)


