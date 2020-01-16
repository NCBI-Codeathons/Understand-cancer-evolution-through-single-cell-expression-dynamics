suppressPackageStartupMessages({
    library(scater)
    library(SingleCellExperiment)
    library(slingshot)
})

snr <- commandArgs()[8]

sample_num = paste0("sample", snr)
data.dir = "/home/ludwig.geistlinger/data"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "slingshot.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

bp <- BiocParallel::registered()[[1]]

## 1. marker genes
cat("Run slingshot ...:\n", file=log.file)

# try: full and filter mat
# try: clusterLabels and celltype labels

reducedDims(sce)$PCA <- reducedDims(sce)$PCA[,1:10]
sce <- slingshot(sce, clusterLabels = "Cluster", reducedDim = "PCA", approx_points=200)

# annotate to sce
ind <- grep("slingPseudotime", colnames(colData(sce)), value=TRUE)

cat(paste0("\n\n Estimated trajectories: ", length(ind)), 
    file=log.file, append=TRUE)

l <- vector("list", length=length(ind))
for(i in seq_along(ind)) l[[i]] <- plotTSNE(sce, colour_by=ind[i])
pdf(file.path(plot.dir, "slingshot.pdf"))
multiplot(plotlist=l, cols=2)
dev.off() 

saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
