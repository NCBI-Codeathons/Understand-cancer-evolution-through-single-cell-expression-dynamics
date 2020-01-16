############################################################
# 
# author: Ludwig Geistlinger
# date: 2020-01-13 11:02:06
# 
# descr: obtain GSE103322 (head and neck cancer scRNA-seq) 
#        from GEO and store in a SingleCellExperiment 
# 
############################################################

library(GEOquery)
library(SummarizedExperiment)
library(SingleCellExperiment)

# Download GSE103322_HNSCC_all_data.txt.gz from GEO
# gunzip GSE103322_HNSCC_all_data.txt.gz
gse <- read.delim("GSE103322_HNSCC_all_data.txt", as.is = TRUE, row.names = 1L)
rownames(gse) <- gsub("\'", "", rownames(gse))

# store in a SingleCellExperiment 
coldat <- DataFrame(t(gse[1:5,]))
for(i in 1:4) coldat[,i] <- as.logical(as.integer(coldat[,i]))
sce <- SingleCellExperiment(assays = list(TPM = gse[6:nrow(gse),]),
                            colData = coldat)

# information about tumors
gse.info <- getGEO("GSE103322")[[1]]
gse.info <- as(gse.info, "SummarizedExperiment")
colData(gse.info)[,11]

sce <- sce[,as.character(colData(gse.info)[,1])]
colData(sce)$tumor_site <- sub("^tumor site: ", "", as.character(colData(gse.info)[,11]))
