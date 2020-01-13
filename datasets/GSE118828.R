############################################################
# 
# author: Ludwig Geistlinger
# date: 2020-01-13 11:02:06
# 
# descr: obtain GSE118828 (ovarian cancer scRNA-seq dataset) 
#        from GEO and store in a SingleCellExperiment 
# 
############################################################

library(GEOquery)
library(SummarizedExperiment)
library(SingleCellExperiment)

# Download GSE118828_RAW.tar from GEO
# untar -xvf GSE118828_RAW.tar

gsm.files <- list.files(".", pattern = "^GSM")

# helper function that reads a GSM file 
# ... and returns a SingleCellExperiment 
readGSM <- function(gsm.file)
{
    gsm <- read.csv(gsm.file, as.is = TRUE, row.names = 1L)
    gsm <- t(gsm)
    SingleCellExperiment(assays = list(counts = gsm))    
}

# obtain scRNA-seq data for a specific tumor as a SingleCellExperiment
sce <- readGSM(gsm.files[1])

# information about tumors
gse <- getGEO("GSE118828")[[1]]
gse <- as(gse, "SummarizedExperiment")
colData(gse)[,10:12]

