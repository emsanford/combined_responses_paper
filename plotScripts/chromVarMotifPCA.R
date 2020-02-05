library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)
library(chromVARmotifs)
library(BiocParallel)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2019)
register(MulticoreParam(4, progressbar = TRUE))

# inputpeaks <- read_tsv(here('extractedData', 'differentialAtacPeaks.tsv'))
fragmentCountsDiffPeaks <- read_rds(here('extractedData', 'atacFragmentCountsAllCondsDifferentialPeaks.rds'))
fragmentCountsDiffPeaks <- addGCBias(fragmentCountsDiffPeaks, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38)
data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper
motif_ix <- matchMotifs(motifSet, fragmentCountsDiffPeaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object = fragmentCountsDiffPeaks, annotations = motif_ix)
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE)
tvar <- as.tibble(variability)

top_n_motifs_to_keep <- 72
TF.names.variabilitycutoff <- arrange(tvar, -variability)$name[1:top_n_motifs_to_keep]
PWM.object.indices.bool <- sapply(strsplit(names(motifSet), "_"), function (x) x[3]) %in% TF.names.variabilitycutoff
PWM.object.indices.num  <- which(PWM.object.indices.bool)
PWM.object.indices.num.1h <- PWM.object.indices.num[1:36]
PWM.object.indices.num.2h <- PWM.object.indices.num[37:72]
# selected.PWM.objects <- motifSet[PWM.object.indices]
# need to filter out the top 72 motifs then do PCA on them, visualize the motifs to categorize as TGFb, RA, or "control/lost" motifs
pcainputmatx1 <- t(as.matrix(deviations(dev[PWM.object.indices.num.1h,])))
pca_res1 <- prcomp(pcainputmatx1)
var_expl1 <- (pca_res1$sdev ^ 2) / (sum(pca_res1$sdev ^ 2))
print(sprintf("Frac var explained PCA set 1: PC1-->%.3f, PC2-->%.3f, rest of PCs-->%.3f", var_expl1[1], var_expl1[2], sum(var_expl1[3:36])))
pcscores1 <- pca_res1$x
pc1_1h <- pcscores1[,1]
pc2_1h <- pcscores1[,2]
textlabels1 <- sapply(strsplit(rownames(pcscores1), "-"), function (x) paste0(x[2], x[3]))
p1 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_1h, y=pc2_1h, label = textlabels1)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# show weights of PC1 and PC2 on a text scatter plot
pcrotation1 <- pca_res1$rotation
pc1_coeffs_1h <- pcrotation1[, 1]
pc2_coeffs_1h <- pcrotation1[, 2]
p2 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_coeffs_1h, y=pc2_coeffs_1h, label = sapply(strsplit(rownames(pcrotation1), "_"), function (x) x[3]))) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)



pcainputmatx2 <- t(as.matrix(deviations(dev[PWM.object.indices.num.2h,])))
pca_res2 <- prcomp(pcainputmatx2)
var_expl2 <- (pca_res2$sdev ^ 2) / (sum(pca_res2$sdev ^ 2))
print(sprintf("Frac var explained PCA set 1: PC1-->%.2f, PC2-->%.2f, rest of PCs-->%.3f", var_expl2[1], var_expl2[2], sum(var_expl2[3:36])))
pcscores2 <- pca_res2$x
pc1_2h <- pcscores2[,1]
pc2_2h <- pcscores2[,2]
textlabels2 <- sapply(strsplit(rownames(pcscores2), "-"), function (x) paste0(x[2], x[3]))
p3 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_2h, y=pc2_2h, label = textlabels2)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# show weights of PC1 and PC2 on a text scatter plot
pcrotation2 <- pca_res2$rotation
pc1_coeffs_2h <- pcrotation2[, 1]
pc2_coeffs_2h <- pcrotation2[, 2]
p4 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_coeffs_2h, y=pc2_coeffs_2h, label = sapply(strsplit(rownames(pcrotation2), "_"), function (x) x[3]))) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)

sort(as.character(TF.names.variabilitycutoff))




