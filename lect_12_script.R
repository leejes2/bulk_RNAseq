#read in metadata for DE
#lecture 12

library(stringr)
library(DESeq2)

myMetaData <- read.csv(file = "~/bulk_RNAseq/bulk_RNAseq_data_qc_metrics.csv", header = T)

#read in counts metrics
load(file = "~/bulk_RNAseq/myCounts.rdat")
geneCounts <- myCounts$counts

#check rownames are gene names in geneCounts
rownames(geneCounts)

#check colnames are sample names
colnames(geneCounts)

#changing the row names using column 1
rownames(myMetaData) <- myMetaData$short_sample_ID
rownames(myMetaData)

#check to see that row names and column names are the same
all(rownames(myMetaData) == colnames(geneCounts))

#remove _sorted.bam from geneCounts name, rename the columns
colnames(geneCounts) <- stringr::str_remove(string = colnames(geneCounts), pattern = "_sorted.bam")

#check names are updated
colnames(geneCounts)
geneCounts <- geneCounts[,rownames(myMetaData)]

#check to see that row names and column names are the same after changing
all(rownames(myMetaData) == colnames(geneCounts))

#check that both dataframes have the same samples, check both ways
all(rownames(myMetaData) %in% colnames(geneCounts))
all(colnames(geneCounts) %in% rownames(myMetaData))




#Step5a: geneCount stats
summary(geneCounts)
dim(geneCounts)


#STEP5b: filtering raw gene counts on low/no expressed gene
keep <- rowSums(geneCounts) >= 15
filteredGeneCounts <- geneCounts[keep, ]


#do it with tidyverse, this is not finished
#geneCounts %>% mutate(sumOfRead = rowSums()) %>% filter(sumofRead>=15)


#Put the data into DESeq2 object
deObj <- DESeqDataSetFromMatrix(countData = filteredGeneCounts, 
                                colData = myMetaData, 
                                design = ~biological_groups)

#STEP5c: normalize counts
deObj <- DESeq(deObj)

#view normalized counts in a new dataframe
normCountsDf <- counts(deObj, normalized = TRUE)
deObj$sizeFactor

#you can always go back to your raw counts by setting the normalized parameter to FALSE
rawCountsDf <- counts(deObj, normalized = FALSE)



#STEP5d: PCA 
pcs <- prcomp(t(normCountsDf), scale. = TRUE)








