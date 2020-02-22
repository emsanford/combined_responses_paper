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
  fragmentCountsDiffPeaks <- read_rds("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/final_diffPeaks_fragment_counts_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
  outputPlotPrefix        <- here("extractedData", "peaks_categorized_by_mode_of_integration", "plots_upreg_alldiff_peaks_")
  siUpregPeaks            <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  selected.PWM.objects    <- readRDS("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/mostVariableMotifs_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
} else {
  fragmentCountsDiffPeaks <- read_rds(cmdargs[1])
  siUpregPeaks <- read_rds(cmdargs[2])
  selected.PWM.objects <- readRDS(cmdargs[3])
  outputPlotPrefix <- cmdargs[4]
}

fragmentCountsDiffPeaks <- addGCBias(fragmentCountsDiffPeaks, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38)

tf.factor.order <- c("RARA", "SMAD3", "SMAD4", "SMAD9", 
                     "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", 
                     "BACH1",  "BACH2", "BATF", "FOXA1", "FOXA2", "FOXA3", "FOXC2", "FOXD3", 
                     "HOXA13", "HOXB13", "HOXC10", "HOXC12", "HOXC13", "HOXD13", 
                     "NFKB1", "REL", "RELA", "SMARCC1", "CDX1",  "CDX2", 
                     "CTCF", "NFE2", "NFE2L2", "MAFF", "MAFK", "BCL11A", "BCL11B", 
                     "GRHL1", "SPI1", "SPIB", "SPIC", 
                     "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5", "ELK4", "ETS2")

#use the selected PWM objects to test for motif deviations in a subset of the data: controls and TGFB only, controls and RA only
allSampleNames <- colnames(fragmentCountsDiffPeaks)
etohSampleNames <- allSampleNames[grepl("EtOH", allSampleNames)]
raSampleNames <- allSampleNames[grepl("-RA-", allSampleNames) & !grepl("-TGFb-", allSampleNames)]
tgfbSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & !grepl("-RA-", allSampleNames)]
BothSampleNames <- allSampleNames[grepl("-TGFb-", allSampleNames) & grepl("-RA-", allSampleNames)]

# calculate deviation scores looking just the RA condition and how it compares to the EtOH condition
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

# calculate deviation scores looking just the TGFb condition and how it compares to the EtOH condition
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

# calculate deviation scores looking at just the "both" condition and how it compares to the EtOH condition
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


combTib <- rbind(raDevTib, TGFbDevTib, BothDevTib)
combTib$tf_name <- factor(combTib$tf_name , levels = tf.factor.order)
devscore.plot.list <- list()
colorscheme <- c('dark green', 'blue', 'dark orange')
for (condname in c("RA", "TGFb", "Both")) {
  p <- combTib %>%
    filter(cond == condname) %>%
    ggplot(aes(tf_name, dev_score)) + 
    geom_bar(stat="identity", fill = colorscheme[length(devscore.plot.list) + 1]) + 
    geom_hline(yintercept = 0) +
    ylab(paste0("dev score ", condname)) + 
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5)) 
  devscore.plot.list[[length(devscore.plot.list) + 1]] <- p
}

patchplot1 <- devscore.plot.list[[1]] / devscore.plot.list[[2]] / devscore.plot.list[[3]]

# make plot of motif matches sub, add, superadd
motif.col.names     <- colnames(siUpregPeaks)[which(grepl("_motifMatchScore", colnames(siUpregPeaks)))]
motif.names         <- sapply(strsplit(motif.col.names, "_"), function(x) x[[1]])
corresponding.avg.d.scores    <- c()
corresponding.median.d.scores <- c()
corresponding.perc.d.scores   <- c()
peak.d.values <- pull(siUpregPeaks, "peakAdditivePredFcResidual-med")
n.peaks       <- length(peak.d.values)

indmotifplots <- list()
counter <- 1
for (motif.col.name in motif.col.names) {
  motif.match.indices <- which(siUpregPeaks[, motif.col.name] > 0)
  frac.peaks.with.motif <- length(motif.match.indices) / n.peaks
  print(sprintf("%s, in %0.3f percent of %d motifs", motif.col.name, frac.peaks.with.motif, n.peaks))
  this.avg.d          <- mean(peak.d.values[motif.match.indices])
  this.median.d       <- median(peak.d.values[motif.match.indices])
  this.80perc.d       <- quantile(peak.d.values[motif.match.indices], .8)
  corresponding.avg.d.scores    <- c(corresponding.avg.d.scores, this.avg.d)
  corresponding.median.d.scores <- c(corresponding.median.d.scores, this.median.d)
  corresponding.perc.d.scores   <- c(corresponding.perc.d.scores, this.80perc.d)
  p <- qplot(peak.d.values[motif.match.indices]) + ggtitle(strsplit(motif.col.name, "_")[[1]][1]) + ylab("counts")
  indmotifplots[[counter]] <- p
  counter <- counter + 1
}
motif.d.scores.tib <- tibble(motif.name = factor(motif.names, levels = tf.factor.order), 
                             median.d.score = corresponding.median.d.scores, 
                             avg.d.score = corresponding.avg.d.scores,
                             perc.d.score = corresponding.perc.d.scores)

motif.d.scores.plot <- motif.d.scores.tib %>%
  ggplot(aes(x = motif.name, y = median.d.score)) +
  geom_bar(stat = "identity") + 
  geom_hline(yintercept = median(peak.d.values)) + 
  xlab("") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5))


# look at frequency of each motif at sub-additive, additive, and super-additive peaks
subadditive.peak.categories   <- c("sub-additive")
additive.peak.categories      <- c("additive", "ambiguous")
superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")

frac.motif.matches.by.category.tib    <- NULL
for (motif.col.name in motif.col.names) {
  subadd.peaks <- siUpregPeaks %>% filter(`peak_integrationCategory-med-dose` %in% subadditive.peak.categories)
  n.subadd.peaks <- nrow(subadd.peaks)
  subadd.motif.match.indices <- which(subadd.peaks[, motif.col.name] > 0)
  subadd.frac.peaks.with.motif <- length(subadd.motif.match.indices) / n.subadd.peaks
  
  add.peaks <- siUpregPeaks %>% filter(`peak_integrationCategory-med-dose` %in% additive.peak.categories)
  n.add.peaks <- nrow(add.peaks)
  add.motif.match.indices <- which(add.peaks[, motif.col.name] > 0)
  add.frac.peaks.with.motif <- length(add.motif.match.indices) / n.add.peaks
  
  superadd.peaks <- siUpregPeaks %>% filter(`peak_integrationCategory-med-dose` %in% superadditive.peak.categories)
  n.superadd.peaks <- nrow(superadd.peaks)
  superadd.motif.match.indices <- which(superadd.peaks[, motif.col.name] > 0)
  superadd.frac.peaks.with.motif <- length(superadd.motif.match.indices) / n.superadd.peaks
  
  motif.name   <- factor(rep(strsplit(motif.col.name, "_")[[1]][1], 3), levels = tf.factor.order)
  frac.peaks.with.motif.matches <- c(subadd.frac.peaks.with.motif, add.frac.peaks.with.motif, superadd.frac.peaks.with.motif)
  peak.category <- factor(c('subadditive', 'additive', 'superadditive'), levels = c('subadditive', 'additive', 'superadditive'))
  peak.category.total.num.peaks <- c(n.subadd.peaks, n.add.peaks, n.superadd.peaks)
  
  frac.motif.matches.by.category.tib <- rbind(frac.motif.matches.by.category.tib, tibble(motif.name, frac.peaks.with.motif.matches, peak.category, peak.category.total.num.peaks))
}


frac.motif.matches.by.category.plot <- frac.motif.matches.by.category.tib %>% ggplot(aes(x = motif.name, y = frac.peaks.with.motif.matches, fill = peak.category)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("") +
  guides(fill=FALSE) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5))


composite.plot <- devscore.plot.list[[1]] / devscore.plot.list[[2]] / devscore.plot.list[[3]] / motif.d.scores.plot / frac.motif.matches.by.category.plot
ggsave(paste0(outputPlotPrefix, "motif_analysis_composite_plot.svg"), plot = composite.plot, width = 24, height = 18)

