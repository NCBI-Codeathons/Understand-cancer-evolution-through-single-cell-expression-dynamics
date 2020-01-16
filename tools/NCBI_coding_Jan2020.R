if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")


### Installing slingshot
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(slingshot)
?read.csv

data1 <- read.csv("/Users/nr267/Desktop/GSM3348303_553_Perit_S1.counts.umiCounts.aboveBackground.table.csv", header = TRUE)
dim (data1)

data1 <- data.matrix(data1)
View (data1)
rm (data1)

### Obtain GSE118828 dataset

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages(("BiocManager"))

library(GEOquery)
library(SummarizedExperiment)
library(SingleCellExperiment)

getwd

setwd("/Users/nr267/Desktop/OVC_data/")

library(GEOquery)
library(SummarizedExperiment)
library(SingleCellExperiment)

gsm.files <- list.files(".", pattern = "^GSM")

#getGEO (GEO = "GSE118828")

gsm.files <- list.files (".", pattern = "^GSM")

readGSM <- function(gsm.file)
{
  gsm <- read.csv(gsm.file, as.is = TRUE, row.names = 1L)
  gsm <- t(gsm)
  SingleCellExperiment(assays = list(counts = gsm))
}

sce <- readGSM(gsm.files[2])
dim(sce)



# information about tumors
gse <- getGEO("GSE118828")[[1]]
gse <- as(gse, "SummarizedExperiment")
colData(gse)[,10:12]

## Install genefilter

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter")

geneFilter <- apply(assays(sce)$counts,1,function(x){sum(x >= 3) >= 10})
sce_genefilter <- sce[geneFilter,]


FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

assays(sce_genefilter)$norm <- FQnorm(assays(sce_genefilter)$counts)

pca <- prcomp(t(log1p(assays(sce_genefilter)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

### Install destiny
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("destiny")

library(destiny)
diffusion_map <- DiffusionMap(t(log1p(assays(sce_genefilter)$norm)))

rd2 <- cbind(DC1 = diffusion_map$DC1, DC2 = diffusion_map$DC2)
plot(rd2, col = rgb(0,0,0,.5), pch = 16,asp = 1)

reducedDims(sce_genefilter) <- SimpleList(PCA = rd1, Diffap = rd2)

### mclust
install.packages("mclust")
library(mclust)

cl1 <- Mclust(rd1)$classification
colData(sce_genefilter)$GMM <- cl1

library (RColorBrewer)

plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce_genefilter)$kmeans <- cl2

plot(rd1, col = brewer.pal(9, "Set1")[cl2], pch = 16, asp = 1)

sce_genefilter_ss <- slingshot (sce_genefilter, clusterLabels = 'GMM', reducedDim = 'PCA')

summary(sce_genefilter_ss$slingPseudotime_1)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_genefilter_ss$slingPseudotime_1, breaks = 100)]

plot(reducedDims(sce_genefilter_ss)$PCA, col = plotcol, pch = 16, asp = 1)
lines(SlingshotDataSet(sce_genefilter_ss), lwd = 2, col = 'black')

plot(reducedDims(sce_genefilter_ss)$PCA, col = brewer.pal(9,'Set1')[sce_genefilter_ss$GMM], pch = 16, asp = 1)
lines(SlingshotDataSet(sce_genefilter_ss), lwd = 2, type = 'lineages', col = 'black')

require(gam)

install.packages(gam)
library(gam)

N <- sce_genefilter_ss$slingPseudotime_1

# for time, only look at the 100 most variable genes
Y <- log1p(assays(sce_genefilter_ss)$norm)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, N=N)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(N), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assays(sce_genefilter_ss)$norm[topgenes, order(N, na.last = NA)]
heatclus <- sce_genefilter_ss$GMM[order(N, na.last = NA)]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

dim(sce_genefilter_ss)

### Install scater
BiocManager::install("scater")
a

### Install SingleR
BiocManager::install("SingleR")

library(scater)
library(SingleR)
library(SingleCellExperiment)

snr <- commandArgs()

sample_num = paste0("sample", snr)
data.dir = "/Users/nr267/Desktop/OVC_data/"
sample.dir = file.path(data.dir, snr)
plot.dir <- file.path(sample.dir, "plots")
log.file <- file.path(sample.dir, "annotate.log")

sce <- readRDS(file.path(sample.dir, paste0("sample", snr, "_sce.rds")))

cat("\n\n Annotate cell types ... \n", file = log.file, append = TRUE)

annotateCellType <- function(sce, ref=c("hpca", "encode"))
{
  ref <- match.arg(ref)
  if (ref == "hpca") se <- HumanPrimaryCellAtlasData()
  else se <- BlueprintEncodeData()
  
  common <- intersect(rownames(sce), rownames(se))
  se <- se[common,]
  sce <- sce[common,]
  
  pred <- SingleR(test = sce, ref = se, labels = se$label.main, assay.type.ref = "logcounts", BPPARAM = BiocParallel::registered()[[1]])


fname <- paste0("sample", snr, "_celltypes_", ref, ".rds")
saveRDS(pred, file = file.path(sample.dir, fname))

return(pred)
}

#annotate and plot

hpc <- annotateCellType(sce, "hpca")
encode <- annotateCellType(sce, "encode")

colData(sce)$hpca.celltype <- hpc$labels
colData(sce)$encode.celltype <- encode$labels

maxScore <- function(res)
  vapply(1:nrow(res), function(i) res$scores[i, res[i, "labels"]], numeric(1))

colData(sce)$hpca.celltype.score <- maxScore(hpc)
colData(sce)$encode.celltype.score <- maxScore(hpc)

plotMainCellTypes <- function(sce, col)
{ 
  main.cats <- sort(table(sce[[col]]), decreasing = TRUE)
  main.cats <- names(main.cats)
  
  sce.sub <- sce[,sce[[col]] %in% main.cats[1:6]]
  sce.sub[[col]] <- factor(as.vector(sce.sub[[col]]), levels=main.cats[1:6])
  plotTSNE(sce.sub, colour_by=col, shape_by=col)
}

pdf(file.path(plot.dir, "cell_types.pdf"))
plotMainCellTypes(sce, "hpca.celltype")
plotMainCellTypes(sce, "encode.celltype")

plotScoreHeatmap(hpc)
plotScoreDistribution(encode)
dev.off()

saveRDS(sce, file = file.path(sample.dir, paste0("sample", snr, "_sce.rds")))
