library(tidyverse)
library(chromVAR)
library(here)
library(chromVARmotifs)
library(BiocParallel)
library(motifmatchr)
library(SummarizedExperiment)
library(patchwork)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
library(BSgenome.Hsapiens.UCSC.hg38)
register(MulticoreParam(4, progressbar = TRUE))

source(here('extractionScripts', 'util.R'))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCountsDiffPeaks <- readRDS("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/final_diffPeaks_fragment_counts_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
  outputPlotPrefix        <- here("plots", "")
} else {
  fragmentCountsDiffPeaks <- readRDS(cmdargs[1])
  outputPlotPrefix        <- cmdargs[2]
}

# tf.factor.order <- c("RARA", "SMAD3", "SMAD4", "SMAD9", 
#                      "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", 
#                      "BACH1",  "BACH2", "BATF", "FOXA1", "FOXA2", "FOXA3", "FOXC2", "FOXD3", 
#                      "HOXA13", "HOXB13", "HOXC10", "HOXC12", "HOXC13", "HOXD13", 
#                      "NFKB1", "REL", "RELA", "SMARCC1", "CDX1",  "CDX2", 
#                      "CTCF", "NFE2", "NFE2L2", "MAFF", "MAFK", "BCL11A", "BCL11B", 
#                      "GRHL1", "SPI1", "SPIB", "SPIC", 
#                      "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5", "ELK4", "ETS2")

n.bootstrap.samples <- 1000

#use the selected PWM objects to test for motif deviations in a subset of the data: controls and TGFB only, controls and RA only
makeDeviationScorePlotTibs <- function(fragmentCountsDiffPeaks, motif_ix, TF.names.of.interest) {
  allSampleNames <- colnames(fragmentCountsDiffPeaks)
  etohSampleNames <- allSampleNames[grepl("EtOH", allSampleNames)]
  raSampleNames <- allSampleNames[grepl("-RA-", allSampleNames) & !grepl("-TGFb-", allSampleNames)]
  tgfbSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & !grepl("-RA-", allSampleNames)]
  BothSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & grepl("-RA-", allSampleNames)]
  
  # calculate deviation scores looking just the RA condition and how it compares to the EtOH condition
  fragCountsRAandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, raSampleNames)]
  devRA <- computeDeviations(object = fragCountsRAandControls, annotations = motif_ix)
  variabilityRA <- computeVariability(devRA)
  tvarRA <- as.tibble(variabilityRA)
  devScoresEtOHconds <- (assays(devRA)$raw_deviations[, 1:9])
  devScoresRAconds   <- (assays(devRA)$raw_deviations[, 10:18])
  avgDevEtOH <- rowMeans(devScoresEtOHconds)
  avgDevRA   <- rowMeans(devScoresRAconds)
  devScoreRA <- (avgDevRA + 1) / (avgDevEtOH + 1) - 1
  tfnames <- names(avgDevEtOH)
  tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
  n.motifs <- length(tfnames)
  n_tfs = length(tfnames)
  raDevTib <- tibble(tf_name = factor(tfnames), dev_score = devScoreRA, chromVARvariability = tvarRA$variability, cond = rep("RA", n_tfs))
  
  print("Retinoic acid conditions motif ranks for motifs of interest:")
  for (tf.name.of.interest in TF.names.of.interest) {
    variabilityRank <- (raDevTib %>% mutate(varRank = rank(chromVARvariability)) %>% filter(tf_name == tf.name.of.interest) %>% pull(varRank))[1] 
    print(sprintf("%s: %d of %d (top %.3f frac)", tf.name.of.interest, variabilityRank, n.motifs, variabilityRank / n.motifs))
  }
  
  # calculate deviation scores looking just the TGFb condition and how it compares to the EtOH condition
  fragCountsTGFbandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, tgfbSampleNames)]
  devTGFb <- computeDeviations(object = fragCountsTGFbandControls, annotations = motif_ix)
  variabilityTGFb <- computeVariability(devTGFb)
  tvarTGFb <- as.tibble(variabilityTGFb)
  devScoresEtOHconds <- (assays(devTGFb)$raw_deviations[, 1:9])
  devScoresTGFbconds   <- (assays(devTGFb)$raw_deviations[, 10:18])
  avgDevEtOH <- rowMeans(devScoresEtOHconds)
  avgDevTGFb   <- rowMeans(devScoresTGFbconds)
  devScoreTGFb <- (avgDevTGFb + 1) / (avgDevEtOH + 1) - 1
  tfnames <- names(avgDevEtOH)
  tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
  n_tfs = length(tfnames)
  TGFbDevTib <- tibble(tf_name = factor(tfnames), dev_score = devScoreTGFb, chromVARvariability = tvarTGFb$variability, cond = rep("TGFb", n_tfs))
  
  print("TGFb conditions motif ranks for motifs of interest:")
  for (tf.name.of.interest in TF.names.of.interest) {
    variabilityRank <- (TGFbDevTib %>% mutate(varRank = rank(chromVARvariability)) %>% filter(tf_name == tf.name.of.interest) %>% pull(varRank))[1] 
    print(sprintf("%s: %d of %d (top %.3f frac)", tf.name.of.interest, variabilityRank, n.motifs, variabilityRank / n.motifs))
  }
  
  
  # calculate deviation scores looking at just the "both" condition and how it compares to the EtOH condition
  fragCountsBothandControls <- fragmentCountsDiffPeaks[, c(etohSampleNames, BothSampleNames)]
  devBoth <- computeDeviations(object = fragCountsBothandControls, annotations = motif_ix)
  variabilityBoth <- computeVariability(devBoth)
  tvarBoth <- as.tibble(variabilityBoth)
  devScoresEtOHconds <- (assays(devBoth)$raw_deviations[, 1:9])
  devScoresBothconds   <- (assays(devBoth)$raw_deviations[, 10:18])
  avgDevEtOH <- rowMeans(devScoresEtOHconds)
  avgDevBoth   <- rowMeans(devScoresBothconds)
  devScoreBoth <- (avgDevBoth + 1) / (avgDevEtOH + 1) - 1
  tfnames <- names(avgDevEtOH)
  tfnames <- sapply( strsplit(tfnames, "_"), function(x) x[3])
  n_tfs = length(tfnames)
  BothDevTib <- tibble(tf_name = factor(tfnames), dev_score = devScoreBoth, chromVARvariability = tvarBoth$variability, cond = rep("Both", n_tfs))

  print("Both signal conditions motif ranks for motifs of interest:")
  for (tf.name.of.interest in TF.names.of.interest) {
    variabilityRank <- (BothDevTib %>% mutate(varRank = rank(chromVARvariability)) %>% filter(tf_name == tf.name.of.interest) %>% pull(varRank))[1] 
    print(sprintf("%s: %d of %d (top %.3f frac)", tf.name.of.interest, variabilityRank, n.motifs, variabilityRank / n.motifs))
  }
  
  
  combTib <- rbind(raDevTib, TGFbDevTib, BothDevTib) 

  return(combTib)  
}

fragmentCountsDiffPeaks <- addGCBias(fragmentCountsDiffPeaks, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38)
data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper

TF.names.of.interest <- c("RARA", "SMAD3", "SMAD4", "SMAD9")

motif_ix_all <- matchMotifs(motifSet, fragmentCountsDiffPeaks, 
                           genome = BSgenome.Hsapiens.UCSC.hg38)


devScoresTib <- makeDeviationScorePlotTibs(fragmentCountsDiffPeaks, motif_ix_all, TF.names.of.interest)

# add bootstrap confidence intervals
set.seed(0)
n.peaks <- nrow(fragmentCountsDiffPeaks)
bootstrap.dev.scores.tib <- NULL
for (ii in 1:n.bootstrap.samples) {
  this.sample.inds <- sample(1:n.peaks, n.peaks, replace = TRUE)
  fragCtsBootstrapSample <- fragmentCountsDiffPeaks[this.sample.inds, ]
  this.motif_ix <- motif_ix_all[this.sample.inds, ]
  this.devScoreTib <- makeDeviationScorePlotTibs(fragCtsBootstrapSample, this.motif_ix, c())
  this.devScores <- this.devScoreTib$dev_score
  bootstrap.dev.scores.tib <- rbind(bootstrap.dev.scores.tib, this.devScores)
}
lower.bootstrap.quantiles <- c()
upper.bootstrap.quantiles <- c()
for (ii in 1:ncol(bootstrap.dev.scores.tib)) {
  lower.bootstrap.quantiles <- c(lower.bootstrap.quantiles, quantile(bootstrap.dev.scores.tib[,ii], .05))
  upper.bootstrap.quantiles <- c(upper.bootstrap.quantiles, quantile(bootstrap.dev.scores.tib[,ii], .95))
}
devScoresTib[["bootstrap_ci_lower"]] <- 2 * devScoresTib$dev_score - upper.bootstrap.quantiles
devScoresTib[["bootstrap_ci_upper"]] <- 2 * devScoresTib$dev_score - lower.bootstrap.quantiles


devscore.plot.list <- list()
colorscheme <- c('dark green', 'blue', 'dark orange')
for (condname in c("RA", "TGFb", "Both")) {
  p <- devScoresTib %>%
    filter(cond == condname,
           tf_name %in% TF.names.of.interest) %>%
    ggplot(aes(tf_name, dev_score, ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper)) + 
    geom_bar(stat="identity", fill = colorscheme[length(devscore.plot.list) + 1]) + 
    geom_errorbar(width = 0) +
    geom_hline(yintercept = 0) +
    ylab(paste0("dev score ", condname)) + 
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5)) 
  devscore.plot.list[[length(devscore.plot.list) + 1]] <- p
}
# brute force, standardize y axes
y_limits <- ggplot_build(devscore.plot.list[[3]])$layout$panel_scales_y[[1]]$range$range
devscore.plot.list[[1]] <- devscore.plot.list[[1]] + ylim(y_limits)
devscore.plot.list[[2]] <- devscore.plot.list[[2]] + ylim(y_limits)
devscore.plot.list[[3]] <- devscore.plot.list[[3]] + ylim(y_limits)
patchplot1 <- devscore.plot.list[[1]] + devscore.plot.list[[2]] + devscore.plot.list[[3]]

ggsave(paste0(outputPlotPrefix, "supp_motif_analysis_canonical_signal_TF_activity.svg"), plot = patchplot1, width = 12, height = 6)

