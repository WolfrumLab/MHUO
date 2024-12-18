# This R script has been used to generate the following figures (separately for subcutaneous and visceral AT depot):
### Figure S1I, S1J

library(corrplot)
library(ezRun)
library(psych)
library(stringr)

# Load numeric clinical parameters
clinParams <- read.table("../MHUO_numClinParams.txt", header=TRUE, sep = "\t")

# Load estimated relative cell type proportions from Bisque deconvolution of bulkRNA
deconv <- read_tsv("../bisqueDeconv.tsv")

# Check sample overlaps 
idx <- which(colnames(clinParams) %in% rownames(deconv))
clinParams <- clinParams[,idx]

# Transpose count matrix
deconv2 <- t(deconv)
deconv2 <- as.matrix(deconv2)

# Calculate Kendall correlation coefficients
corrRes <- psych::corr.test(clinParams, deconv2, method="kendall", adjust="BH")

corcoeff <- corrRes$r
cor_pVal_adj <- corrRes$p.adj

col<- colorRampPalette(c("blue", "white", "red"))(20)

# Correlation plot
corrplot(t(corcoeff), col = col,is.corr=FALSE, 
         tl.col = "black", tl.srt = 45, tl.cex = 0.8, col.lim = c(-1, 1), p.mat = t(cor_pVal_adj), 
         sig.level = 0.01, pch.cex = 0.9,
         insig = 'blank', pch.col = 'grey20', cl.ratio = 0.08, cl.length = 5)

##########################################################################################################################################

# This R script has been used to generate the following figures:
### Figure 5K, S8A, S8B

# Load normalized counts
ds <- read_tsv("../visAT_MHUO_normCounts.tsv")

# List genes of interest
candGenes <- c("KRT8", "KRT19", "TM4SF1", "CDON", "WT1", "BCHE", "DPP4", "DIAPH3", "CD34", "ITLN1")

# Subset normalized count matrix for genes of interest 
ds_sub <- subset(ds, gene_name %in% candGenes)

rownames(ds_sub) <- NULL
ds_sub <- column_to_rownames(ds_sub, "gene_name")

# Check sample overlaps 
idx <- which(colnames(clinParams) %in% rownames(ds_sub))
clinParams <- clinParams[,idx]

# Transpose count matrix
ds_sub2 <- t(ds_sub)
ds_sub2 <- as.matrix(ds_sub2)

# Calculate Spearman correlation coefficients
corrRes <- psych::corr.test(clinParams, ds_sub2, method="spearman")

corcoeff <- corrRes$r
cor_pVal_adj <- corrRes$p

col<- colorRampPalette(c("blue", "white", "red"))(20)

# Correlation plot
corrplot(t(corcoeff), col = col,is.corr=FALSE,
         tl.col = "black", tl.srt = 45, tl.cex = 0.8, col.lim = c(-0.5, 0.5), p.mat = t(cor_pVal_adj),
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', cl.ratio = 0.08, cl.length = 5)

