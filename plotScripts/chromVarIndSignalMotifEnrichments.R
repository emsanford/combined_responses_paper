library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)
library(chromVARmotifs)
library(BiocParallel)
library(motifmatchr)
library(SummarizedExperiment)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
library(BSgenome.Hsapiens.UCSC.hg38)
register(MulticoreParam(4, progressbar = TRUE))

# I had to change some plumbing in chromVAR to get access to raw deviation scores. I did the following steps, some of which must have been required:
#    0. clone the chromVAR git repo, edit the compute_deviations source code file, and commit the changes to git. (but do not push)
#    1. uninstall previous chromVAR installation with remove.packages(c("chromVAR"))
#    2. restart R Studio
#    3. install the edited local package using
        # library(devtools)
        # install_local(path = "/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/software/modified_chromVAR/chromVAR", force = TRUE, build = TRUE)
        # (enter an empty line when prompted to update a bunch of package versions)


# inputpeaks <- read_tsv(here('extractedData', 'differentialAtacPeaks.tsv'))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCountsDiffPeaks <- read_rds(here('extractedData', 'atacFragmentCountsAllCondsDifferentialPeaks.rds'))
  outputPlotPrefix <- here("extractedData", "peaks_categorized_by_mode_of_integration", "plots_upreg_alldiff_peaks_")
} else {
  fragmentCountsDiffPeaks <- read_rds(cmdargs[1])
  outputPlotPrefix <- cmdargs[2]
}



# fragmentCountsDiffPeaks <- read_rds(here('extractedData', 'peaks_categorized_by_mode_of_integration', 'upreg_add_peaks.rds'))
# outputPlotPrefix <- here("extractedData", "peaks_categorized_by_mode_of_integration", "plots_upreg_add_peaks_")
# fragmentCountsDiffPeaks <- read_rds(here('extractedData', 'peaks_categorized_by_mode_of_integration', 'upreg_subadd_peaks.rds'))
# outputPlotPrefix <- here("extractedData", "peaks_categorized_by_mode_of_integration", "plots_upreg_subadd_peaks_")

# fragmentCountsDiffPeaks1 <- read_rds(here('extractedData', 'peaks_categorized_by_mode_of_integration', 'upreg_superadd_peaks.rds'))
# outputPlotPrefix <- here("extractedData", "peaks_categorized_by_mode_of_integration", "plots_upreg_superadd_peaks_")
# fragmentCountsDiffPeaks2 <- read_rds(here('extractedData', 'peaks_categorized_by_mode_of_integration', 'upreg_add_peaks.rds'))
# fragmentCountsDiffPeaks <- rbind(fragmentCountsDiffPeaks1, fragmentCountsDiffPeaks2[1:28,])  ## be sure to fix this "dirtying the super-add peaks" short-fix

selected.PWM.objects <- readRDS(here("extractedData", "peaks_categorized_by_mode_of_integration", "top72motifs_alldiffpeaks_robject.rds"))

fragmentCountsDiffPeaks <- addGCBias(fragmentCountsDiffPeaks, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38)
##### select the top 72 most variable motifs from the data set #####
# set.seed(2019)
# data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
# motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper
# motif_ix <- matchMotifs(motifSet, fragmentCountsDiffPeaks, 
#                         genome = BSgenome.Hsapiens.UCSC.hg38)
# dev <- computeDeviations(object = fragmentCountsDiffPeaks, annotations = motif_ix)
# variability <- computeVariability(dev)
# plotVariability(variability, use_plotly = FALSE)
# tvar <- as.tibble(variability)
# 
# top_n_motifs_to_keep <- 72
# TF.names.variabilitycutoff <- arrange(tvar, -variability)$name[1:top_n_motifs_to_keep]
# PWM.object.indices.bool <- sapply(strsplit(names(motifSet), "_"), function (x) x[3]) %in% TF.names.variabilitycutoff
# PWM.object.indices.num  <- which(PWM.object.indices.bool)
# PWM.object.indices.num.1h <- PWM.object.indices.num[1:36]
# PWM.object.indices.num.2h <- PWM.object.indices.num[37:72]
# selected.PWM.objects <- motifSet[PWM.object.indices.num]
################################################################################################################  

# now, use the selected PWM objects to test for motif deviations in a subset of the data: controls and TGFB only, controls and RA only
allSampleNames <- colnames(fragmentCountsDiffPeaks)
etohSampleNames <- allSampleNames[grepl("EtOH", allSampleNames)]
raSampleNames <- allSampleNames[grepl("-RA-", allSampleNames) & !grepl("-TGFb-", allSampleNames)]
tgfbSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & !grepl("-RA-", allSampleNames)]
BothSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & grepl("-RA-", allSampleNames)]

set.seed(2019)
fragCountsRAandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, raSampleNames)]
motif_ix <- matchMotifs(selected.PWM.objects, fragCountsRAandControls, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devRA <- computeDeviations(object = fragCountsRAandControls, annotations = motif_ix)
variabilityRA <- computeVariability(devRA)
plotVariability(variabilityRA, use_plotly = FALSE)
tvarRA <- as.tibble(variabilityRA)

devScoresEtOHconds <- deviations(devRA[, 1:9])
devScoresRAconds   <- deviations(devRA[, 10:18])
avgDevEtOH <- rowMeans(devScoresEtOHconds)
avgDevRA   <- rowMeans(devScoresRAconds)
tfnames <- names(avgDevEtOH)
tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
n_tfs = length(tfnames)
raDevTib <- tibble(tf_name = factor(rep(tfnames, 2)), dev_score = c(avgDevEtOH, avgDevRA), cond = c(rep("EtOH-1", n_tfs), rep("RA", n_tfs)))
ggplot(raDevTib, aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .)
# deviationScores(devRA[, 10:18])

set.seed(2019)
fragCountsTGFbandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, tgfbSampleNames)]
motif_ix <- matchMotifs(selected.PWM.objects, fragCountsTGFbandControls, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devTGFb <- computeDeviations(object = fragCountsTGFbandControls, annotations = motif_ix)
variabilityTGFb <- computeVariability(devTGFb)
plotVariability(variabilityTGFb, use_plotly = FALSE)
tvarTGFb <- as.tibble(variabilityTGFb)

devScoresEtOHconds <- deviations(devTGFb[, 1:9])
devScoresTGFbconds   <- deviations(devTGFb[, 10:18])
avgDevEtOH <- rowMeans(devScoresEtOHconds)
avgDevTGFb   <- rowMeans(devScoresTGFbconds)
tfnames <- names(avgDevEtOH)
tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
n_tfs = length(tfnames)
TGFbDevTib <- tibble(tf_name = factor(rep(tfnames, 2)), dev_score = c(avgDevEtOH, avgDevTGFb), cond = c(rep("EtOH-2", n_tfs), rep("TGFb", n_tfs)))
ggplot(TGFbDevTib, aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .)


set.seed(2019)
fragCountsBothandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, BothSampleNames)]
motif_ix <- matchMotifs(selected.PWM.objects, fragCountsBothandControls, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devBoth <- computeDeviations(object = fragCountsBothandControls, annotations = motif_ix)
variabilityBoth <- computeVariability(devBoth)
plotVariability(variabilityBoth, use_plotly = FALSE)
tvarBoth <- as.tibble(variabilityBoth)

devScoresEtOHconds <- deviations(devBoth[, 1:9])
devScoresBothconds   <- deviations(devBoth[, 10:18])
avgDevEtOH <- rowMeans(devScoresEtOHconds)
avgDevBoth   <- rowMeans(devScoresBothconds)
tfnames <- names(avgDevEtOH)
tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
n_tfs = length(tfnames)
BothDevTib <- tibble(tf_name = factor(rep(tfnames, 2)), dev_score = c(avgDevEtOH, avgDevBoth), cond = c(rep("EtOH-2", n_tfs), rep("Both", n_tfs)))
ggplot(BothDevTib, aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .)


combTib1 <- rbind(raDevTib, TGFbDevTib)
ggplot(filter(combTib1, cond %in% c("TGFb", "RA")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))

combTib2 <- rbind(raDevTib, TGFbDevTib, BothDevTib)
ggplot(filter(combTib2, cond %in% c("TGFb", "RA", "Both")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))

combTib3 <- rbind(raDevTib, TGFbDevTib, BothDevTib, tibble(tf_name = factor(tfnames), dev_score = (avgDevTGFb + avgDevRA), cond = rep("add dev score pred", n_tfs)))
ggplot(filter(combTib3, cond %in% c("TGFb", "RA", "Both", "add dev score pred")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))

combTib4 <- rbind(raDevTib, TGFbDevTib, BothDevTib, 
                  tibble(tf_name = factor(tfnames), dev_score = (avgDevTGFb + avgDevRA), cond = rep("add dev score pred", n_tfs)),
                  tibble(tf_name = factor(tfnames), dev_score = avgDevBoth - (avgDevTGFb + avgDevRA), cond = rep("c-residual", n_tfs)))
p <- ggplot(filter(combTib4, cond %in% c("TGFb", "RA", "Both", "add dev score pred", "c-residual")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))
ggsave(paste0(outputPlotPrefix, "bias_corrected_dev_score_by_tf_name.svg"), width = 12, height = 4.5, plot = p)


############ this is the exact same code as above, but pasted and using the raw deviation scores instead of bias-corrected deviation scores ############
allSampleNames <- colnames(fragmentCountsDiffPeaks)
etohSampleNames <- allSampleNames[grepl("EtOH", allSampleNames)]
raSampleNames <- allSampleNames[grepl("-RA-", allSampleNames) & !grepl("-TGFb-", allSampleNames)]
tgfbSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & !grepl("-RA-", allSampleNames)]
BothSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & grepl("-RA-", allSampleNames)]

set.seed(2019)
fragCountsRAandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, raSampleNames)]
motif_ix <- matchMotifs(selected.PWM.objects, fragCountsRAandControls, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devRA <- computeDeviations(object = fragCountsRAandControls, annotations = motif_ix)
variabilityRA <- computeVariability(devRA)
plotVariability(variabilityRA, use_plotly = FALSE)
tvarRA <- as.tibble(variabilityRA)

devScoresEtOHconds <- (assays(devRA)$raw_deviations[, 1:9])
devScoresRAconds   <- (assays(devRA)$raw_deviations[, 10:18])
avgDevEtOH <- rowMeans(devScoresEtOHconds)
avgDevRA   <- rowMeans(devScoresRAconds)
tfnames <- names(avgDevEtOH)
tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
n_tfs = length(tfnames)
raDevTib <- tibble(tf_name = factor(rep(tfnames, 2)), dev_score = c(avgDevEtOH, avgDevRA), cond = c(rep("EtOH-1", n_tfs), rep("RA", n_tfs)))
ggplot(raDevTib, aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .)
# deviationScores(devRA[, 10:18])

set.seed(2019)
fragCountsTGFbandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, tgfbSampleNames)]
motif_ix <- matchMotifs(selected.PWM.objects, fragCountsTGFbandControls, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devTGFb <- computeDeviations(object = fragCountsTGFbandControls, annotations = motif_ix)
variabilityTGFb <- computeVariability(devTGFb)
plotVariability(variabilityTGFb, use_plotly = FALSE)
tvarTGFb <- as.tibble(variabilityTGFb)

devScoresEtOHconds <- (assays(devTGFb)$raw_deviations[, 1:9])
devScoresTGFbconds   <- (assays(devTGFb)$raw_deviations[, 10:18])
avgDevEtOH <- rowMeans(devScoresEtOHconds)
avgDevTGFb   <- rowMeans(devScoresTGFbconds)
tfnames <- names(avgDevEtOH)
tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
n_tfs = length(tfnames)
TGFbDevTib <- tibble(tf_name = factor(rep(tfnames, 2)), dev_score = c(avgDevEtOH, avgDevTGFb), cond = c(rep("EtOH-2", n_tfs), rep("TGFb", n_tfs)))
ggplot(TGFbDevTib, aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .)


set.seed(2019)
fragCountsBothandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, BothSampleNames)]
motif_ix <- matchMotifs(selected.PWM.objects, fragCountsBothandControls, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devBoth <- computeDeviations(object = fragCountsBothandControls, annotations = motif_ix)
variabilityBoth <- computeVariability(devBoth)
plotVariability(variabilityBoth, use_plotly = FALSE)
tvarBoth <- as.tibble(variabilityBoth)

devScoresEtOHconds <- (assays(devBoth)$raw_deviations[, 1:9])
devScoresBothconds   <- (assays(devBoth)$raw_deviations[, 10:18])
avgDevEtOH <- rowMeans(devScoresEtOHconds)
avgDevBoth   <- rowMeans(devScoresBothconds)
tfnames <- names(avgDevEtOH)
tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
n_tfs = length(tfnames)
BothDevTib <- tibble(tf_name = factor(rep(tfnames, 2)), dev_score = c(avgDevEtOH, avgDevBoth), cond = c(rep("EtOH-2", n_tfs), rep("Both", n_tfs)))
ggplot(BothDevTib, aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .)


combTib1 <- rbind(raDevTib, TGFbDevTib)
ggplot(filter(combTib1, cond %in% c("TGFb", "RA")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))

combTib2 <- rbind(raDevTib, TGFbDevTib, BothDevTib)
ggplot(filter(combTib2, cond %in% c("TGFb", "RA", "Both")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))

combTib3 <- rbind(raDevTib, TGFbDevTib, BothDevTib, tibble(tf_name = factor(tfnames), dev_score = (avgDevTGFb + avgDevRA), cond = rep("add dev score pred", n_tfs)))
ggplot(filter(combTib3, cond %in% c("TGFb", "RA", "Both", "add dev score pred")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))

combTib4 <- rbind(raDevTib, TGFbDevTib, BothDevTib, 
                  tibble(tf_name = factor(tfnames), dev_score = (avgDevTGFb + avgDevRA), cond = rep("add dev score pred", n_tfs)),
                  tibble(tf_name = factor(tfnames), dev_score = avgDevBoth - (avgDevTGFb + avgDevRA), cond = rep("c-residual", n_tfs)))
p <- ggplot(filter(combTib4, cond %in% c("TGFb", "RA", "Both", "add dev score pred", "c-residual")), aes(tf_name, dev_score, color = cond)) + geom_bar(stat="identity") + facet_grid(cond ~ .) + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))
ggsave(paste0(outputPlotPrefix, "raw_dev_score_by_tf_name.svg"), width = 12, height = 4.5, plot = p)

########################################################################################################################################################

# # scatterplot of RA deviation vs TGFb deviation:
# ggplot(tibble(avgDevTGFb, avgDevRA, labels = tfnames), aes(avgDevTGFb, avgDevRA, label = tfnames)) + geom_text() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + coord_fixed()
# 
# # scatterplot of Both treatment vs. prediction
# ggplot(tibble(avgDevAddPred = (avgDevTGFb + avgDevRA), avgDevBoth, labels = tfnames), aes(avgDevAddPred, avgDevBoth, label = tfnames)) + geom_point() + geom_abline(slope = 1) + coord_fixed()

# 
# pcainputmatx1 <- t(as.matrix(deviations(dev[PWM.object.indices.num.1h,])))
# pca_res1 <- prcomp(pcainputmatx1)
# var_expl1 <- (pca_res1$sdev ^ 2) / (sum(pca_res1$sdev ^ 2))
# print(sprintf("Frac var explained PCA set 1: PC1-->%.3f, PC2-->%.3f, rest of PCs-->%.3f", var_expl1[1], var_expl1[2], sum(var_expl1[3:36])))
# pcscores1 <- pca_res1$x
# pc1_1h <- pcscores1[,1]
# pc2_1h <- pcscores1[,2]
# textlabels1 <- sapply(strsplit(rownames(pcscores1), "-"), function (x) paste0(x[2], x[3]))
# p1 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_1h, y=pc2_1h, label = textlabels1)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# # show weights of PC1 and PC2 on a text scatter plot
# pcrotation1 <- pca_res1$rotation
# pc1_coeffs_1h <- pcrotation1[, 1]
# pc2_coeffs_1h <- pcrotation1[, 2]
# p2 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_coeffs_1h, y=pc2_coeffs_1h, label = sapply(strsplit(rownames(pcrotation1), "_"), function (x) x[3]))) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# 
# 
# 
# pcainputmatx2 <- t(as.matrix(deviations(dev[PWM.object.indices.num.2h,])))
# pca_res2 <- prcomp(pcainputmatx2)
# var_expl2 <- (pca_res2$sdev ^ 2) / (sum(pca_res2$sdev ^ 2))
# print(sprintf("Frac var explained PCA set 1: PC1-->%.2f, PC2-->%.2f, rest of PCs-->%.3f", var_expl2[1], var_expl2[2], sum(var_expl2[3:36])))
# pcscores2 <- pca_res2$x
# pc1_2h <- pcscores2[,1]
# pc2_2h <- pcscores2[,2]
# textlabels2 <- sapply(strsplit(rownames(pcscores2), "-"), function (x) paste0(x[2], x[3]))
# p3 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_2h, y=pc2_2h, label = textlabels2)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# # show weights of PC1 and PC2 on a text scatter plot
# pcrotation2 <- pca_res2$rotation
# pc1_coeffs_2h <- pcrotation2[, 1]
# pc2_coeffs_2h <- pcrotation2[, 2]
# p4 <- ggplot(NULL) + geom_text(mapping = aes(x=pc1_coeffs_2h, y=pc2_coeffs_2h, label = sapply(strsplit(rownames(pcrotation2), "_"), function (x) x[3]))) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
# library(gridExtra)
# grid.arrange(p1, p2, p3, p4, ncol = 2)
# 
# sort(as.character(TF.names.variabilitycutoff))
# 



