#read in metadata for DE
#lecture 12

library(stringr)

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

