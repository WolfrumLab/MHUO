# This R script has been used to generate the following figures (separately for subcutaneous and visceral AT depot):
### Figure 1M, 1N, 2L, 2M, 5C

library(Biobase)
library(BisqueRNA)
library(DESeq2)
library(ezRun)
library(limma)
library(Seurat)
library(tidyverse)

# Load bulkRNA raw count matrix as a SummarizedExperiment object
rawData <- readRDS("../rawData.rds")

# Extract sample metadata
param <- metadata(rawData)$param
dataset <- data.frame(colData(rawData), row.names=colnames(rawData),
                      check.names = FALSE)

# Add sample condition info
metadata(rawData)$design <- ezDesignFromDataset(dataset, param)
group <- metadata(rawData)$design

# Remove undetected genes
rawCounts <- assays(rawData)$counts
rawCounts <- rawCounts[,rownames(group)]

rawCounts <- rawCounts[rowSums(rawCounts) > 0, ]

# Extract gene annotations
geneAnnotation <- do.call(cbind.data.frame, elementMetadata(rawData))

# Create DGEList object
dds <- DESeqDataSetFromMatrix(countData = rawCounts2,
                              colData = group,
                              design= ~ 1)

# Remove genes with < 10 reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Apply variance stabilizing transformation
vsd <- vst(dds)

# Adjust for exon mapping rate (EMR)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), covariates=cbind(vsd$EMR))
normCounts <- assay(vsd) 

# Combine normalized counts with gene annotations
normCounts_df <- merge(geneAnnotation, as.data.frame(normCounts), by.x="gene_id", by.y="row.names")

write_tsv(normCounts, "MHUO_normCounts.tsv")

##########################################################################################################################################

# Convert normalized count matrix to eSet
bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(normCounts))

# Load integrated snRNAseq Seurat object
sndata <- readRDS("../snMulti.rds")

# Convert Seurat object to eSet
sn.eset <- BisqueRNA::SeuratToExpressionSet(sndata, delimiter="-", position=2, version="v3")

# Perform snRNA reference based deconvolution of bulkRNA
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sn.eset, markers=NULL, use.overlap = FALSE)

ref.based.estimates <- res$bulk.props
knitr::kable(ref.based.estimates, digits=2)

ref.based.estimates <- data.frame(ref.based.estimates, check.names = FALSE)
ref.based.estimates <- rownames_to_column(ref.based.estimates, "cellType")

# Save information about estimated relative amount of cells 
write_tsv(ref.based.estimates, "bisqueDeconv.tsv")
