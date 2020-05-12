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
  outputPlotPrefix        <- here("extractedData", "peaks_categorized_by_mode_of_integration", "plots_upreg_alldiff_peaks_")
  siUpregPeaks            <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  selected.PWM.objects    <- readRDS("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/mostVariableMotifs_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
} else {
  fragmentCountsDiffPeaks <- readRDS(cmdargs[1])
  siUpregPeaks            <- read_tsv(cmdargs[2])
  selected.PWM.objects    <- readRDS(cmdargs[3])
  outputPlotPrefix        <- cmdargs[4]
}

tf.factor.order <- c("RARA", "SMAD3", "SMAD4", "SMAD9", 
                     "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", 
                     "BACH1",  "BACH2", "BATF", "FOXA1", "FOXA2", "FOXA3", "FOXC2", "FOXD3", 
                     "HOXA13", "HOXB13", "HOXC10", "HOXC12", "HOXC13", "HOXD13", 
                     "NFKB1", "REL", "RELA", "SMARCC1", "CDX1",  "CDX2", 
                     "CTCF", "NFE2", "NFE2L2", "MAFF", "MAFK", "BCL11A", "BCL11B", 
                     "GRHL1", "SPI1", "SPIB", "SPIC", 
                     "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5", "ELK4", "ETS2")

n.bootstrap.samples <- 1000

#use the selected PWM objects to test for motif deviations in a subset of the data: controls and TGFB only, controls and RA only
makeDeviationScorePlotTibs <- function(fragmentCountsDiffPeaks, motif_ix) {
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
  n_tfs = length(tfnames)
  raDevTib <- tibble(tf_name = factor(tfnames, levels = tf.factor.order), dev_score = devScoreRA, cond = rep("RA", n_tfs))
  
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
  TGFbDevTib <- tibble(tf_name = factor(tfnames, levels = tf.factor.order), dev_score = devScoreTGFb, cond = rep("TGFb", n_tfs))
  
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
  BothDevTib <- tibble(tf_name = factor(tfnames, levels = tf.factor.order), dev_score = devScoreBoth, cond = rep("Both", n_tfs))

  combTib <- rbind(raDevTib, TGFbDevTib, BothDevTib) 

  return(combTib)  
}

fragmentCountsDiffPeaks <- addGCBias(fragmentCountsDiffPeaks, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38)
motif_ix <- matchMotifs(selected.PWM.objects, fragmentCountsDiffPeaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
devScoresTib <- makeDeviationScorePlotTibs(fragmentCountsDiffPeaks, motif_ix)

# add bootstrap confidence intervals
set.seed(0)
n.peaks <- nrow(fragmentCountsDiffPeaks)
bootstrap.dev.scores.tib <- NULL
for (ii in 1:n.bootstrap.samples) {
  this.sample.inds <- sample(1:n.peaks, n.peaks, replace = TRUE)
  fragCtsBootstrapSample <- fragmentCountsDiffPeaks[this.sample.inds, ]
  this.motif_ix <- motif_ix[this.sample.inds, ]
  this.devScoreTib <- makeDeviationScorePlotTibs(fragCtsBootstrapSample, this.motif_ix)
  this.devScores <- this.devScoreTib$dev_score
  bootstrap.dev.scores.tib <- rbind(bootstrap.dev.scores.tib, this.devScores)
  
  print(sprintf("%d of %d bootstrap samples complete for motif enrichment by TF plot", ii, n.bootstrap.samples))
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
    filter(cond == condname) %>%
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
patchplot1 <- devscore.plot.list[[1]] / devscore.plot.list[[2]] / devscore.plot.list[[3]]


##### make median D score by motif plot
makeMotifDScoresTib <- function(anno.peaks, motif.column.names) {
  motif.names         <- sapply(strsplit(motif.column.names, "_"), function(x) x[[1]])
  peak.d.values <- pull(anno.peaks, "peakAdditivePredFcResidual-med")
  n.peaks       <- length(peak.d.values)
  corresponding.avg.d.scores    <- c()
  corresponding.median.d.scores <- c()
  for (motif.col.name in motif.column.names) {
    motif.match.indices <- which(anno.peaks[, motif.col.name] > 0)
    frac.peaks.with.motif <- length(motif.match.indices) / n.peaks
    this.avg.d          <- mean(peak.d.values[motif.match.indices])
    this.median.d       <- median(peak.d.values[motif.match.indices])
    corresponding.avg.d.scores    <- c(corresponding.avg.d.scores, this.avg.d)
    corresponding.median.d.scores <- c(corresponding.median.d.scores, this.median.d)
    p <- qplot(peak.d.values[motif.match.indices]) + ggtitle(strsplit(motif.col.name, "_")[[1]][1]) + ylab("counts")
  }
  motif.d.scores.tib <- tibble(motif.name = factor(motif.names, levels = tf.factor.order), 
                               median.d.score = corresponding.median.d.scores, 
                               avg.d.score = corresponding.avg.d.scores)
  return(motif.d.scores.tib)
}

motif.col.names     <- colnames(siUpregPeaks)[which(grepl("_motifMatchScore", colnames(siUpregPeaks)))]
motif.names         <- sapply(strsplit(motif.col.names, "_"), function(x) x[[1]])
corresponding.avg.d.scores    <- c()
corresponding.median.d.scores <- c()
corresponding.perc.d.scores   <- c()

motif.d.scores.tib <- makeMotifDScoresTib(siUpregPeaks, motif.col.names)

# now do bootstrap error bars
set.seed(0)
median.dscores.tib <- NULL
for (ii in 1:n.bootstrap.samples) {
  siUpregPeaksBootstrapSample   <- sample_n(siUpregPeaks, nrow(siUpregPeaks), replace = TRUE)
  this.sample.motif.d.score.tib <- makeMotifDScoresTib(siUpregPeaksBootstrapSample, motif.col.names)
  this.median.dscores <- this.sample.motif.d.score.tib$median.d.score
  median.dscores.tib <- rbind(median.dscores.tib, this.median.dscores)
}
colnames(median.dscores.tib) <- this.sample.motif.d.score.tib$motif.name
rownames(median.dscores.tib) <- NULL

lower.bootstrap.quantiles <- c()
upper.bootstrap.quantiles <- c()
for (ii in 1:ncol(median.dscores.tib)) {
  lower.bootstrap.quantiles <- c(lower.bootstrap.quantiles, quantile(median.dscores.tib[,ii], .05))
  upper.bootstrap.quantiles <- c(upper.bootstrap.quantiles, quantile(median.dscores.tib[,ii], .95))
}
motif.d.scores.tib[["bootstrap_ci_lower"]] <- 2 * motif.d.scores.tib$median.d.score - upper.bootstrap.quantiles
motif.d.scores.tib[["bootstrap_ci_upper"]] <- 2 * motif.d.scores.tib$median.d.score - lower.bootstrap.quantiles

motif.d.scores.plot <- motif.d.scores.tib %>%
  ggplot(aes(x = motif.name, y = median.d.score, ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper)) +
  geom_bar(stat = "identity") + 
  geom_errorbar(width = 0) +
  geom_hline(yintercept = median(pull(siUpregPeaks, "peakAdditivePredFcResidual-med")), color = "grey") + 
  geom_hline(yintercept = 0, color = "black") + 
  xlab("") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5))

##### make number of motif matches per 150bp at different classes of peaks: sub-additive, additive, and super-additive
makeFracMotifMatchesByCategoryTib <- function(anno.peaks, motif.col.names, subadditive.peak.categories, additive.peak.categories, superadditive.peak.categories) {
  frac.motif.matches.by.category.tib    <- NULL
  for (motif.col.name in motif.col.names) {
    subadd.peaks <- anno.peaks %>% filter(`peak_integrationCategory-med-dose` %in% subadditive.peak.categories)
    n.subadd.peaks <- nrow(subadd.peaks)
    n.subadd.motif.matches <- sum(subadd.peaks[, motif.col.name])
    total.subadd.peak.sequence <- sum(subadd.peaks$endLocs - subadd.peaks$startLocs + 1)
    subadd.num.motifs.per.150bp <- 150 * n.subadd.motif.matches / total.subadd.peak.sequence
    
    add.peaks <- anno.peaks %>% filter(`peak_integrationCategory-med-dose` %in% additive.peak.categories)
    n.add.peaks <- nrow(add.peaks)
    n.add.motif.matches <- sum(add.peaks[, motif.col.name])
    total.add.peak.sequence <- sum(add.peaks$endLocs - add.peaks$startLocs + 1)
    add.num.motifs.per.150bp <- 150 * n.add.motif.matches / total.add.peak.sequence
    
    superadd.peaks <- anno.peaks %>% filter(`peak_integrationCategory-med-dose` %in% superadditive.peak.categories)
    n.superadd.peaks <- nrow(superadd.peaks)
    n.superadd.motif.matches <- sum(superadd.peaks[, motif.col.name])
    total.superadd.peak.sequence <- sum(superadd.peaks$endLocs - superadd.peaks$startLocs + 1)
    superadd.num.motifs.per.150bp <- 150 * n.superadd.motif.matches / total.superadd.peak.sequence
    
    motif.name   <- factor(rep(strsplit(motif.col.name, "_")[[1]][1], 3), levels = tf.factor.order)
    num.motifs.per.150bp.matches <- c(subadd.num.motifs.per.150bp, add.num.motifs.per.150bp, superadd.num.motifs.per.150bp)
    peak.category <- factor(c('subadditive', 'additive', 'superadditive'), levels = c('subadditive', 'additive', 'superadditive'))
    peak.category.total.num.peaks <- c(n.subadd.peaks, n.add.peaks, n.superadd.peaks)
    
    frac.motif.matches.by.category.tib <- rbind(frac.motif.matches.by.category.tib, tibble(motif.name, num.motifs.per.150bp.matches, peak.category, peak.category.total.num.peaks))
  }
  return(frac.motif.matches.by.category.tib)
}

subadditive.peak.categories   <- c("sub-additive")
additive.peak.categories      <- c("additive", "ambiguous")
superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")

motif.col.names     <- colnames(siUpregPeaks)[which(grepl("_numMotifMatches", colnames(siUpregPeaks)))]
frac.motif.matches.by.category.tib <- makeFracMotifMatchesByCategoryTib(siUpregPeaks, motif.col.names, subadditive.peak.categories, additive.peak.categories, superadditive.peak.categories)

# now do bootstrap error bars
set.seed(0)
num.motifs.per.150bp.match.tib <- NULL
for (ii in 1:n.bootstrap.samples) {
  siUpregPeaksBootstrapSample   <- sample_n(siUpregPeaks, nrow(siUpregPeaks), replace = TRUE)
  this.frac.motif.matches.by.category <- makeFracMotifMatchesByCategoryTib(siUpregPeaksBootstrapSample, motif.col.names, subadditive.peak.categories, additive.peak.categories, superadditive.peak.categories)
  this.num.motifs.per.150bp.match <- this.frac.motif.matches.by.category$num.motifs.per.150bp.matches
  num.motifs.per.150bp.match.tib <- rbind(num.motifs.per.150bp.match.tib, this.num.motifs.per.150bp.match)
}
colnames(num.motifs.per.150bp.match.tib) <- paste0(this.frac.motif.matches.by.category$motif.name, "_", this.frac.motif.matches.by.category$peak.category)
rownames(num.motifs.per.150bp.match.tib) <- NULL

lower.bootstrap.quantiles <- c()
upper.bootstrap.quantiles <- c()
for (ii in 1:ncol(num.motifs.per.150bp.match.tib)) {
  lower.bootstrap.quantiles <- c(lower.bootstrap.quantiles, quantile(num.motifs.per.150bp.match.tib[,ii], .05))
  upper.bootstrap.quantiles <- c(upper.bootstrap.quantiles, quantile(num.motifs.per.150bp.match.tib[,ii], .95))
}
frac.motif.matches.by.category.tib[["bootstrap_ci_lower"]] <- 2 * frac.motif.matches.by.category.tib$num.motifs.per.150bp.matches - upper.bootstrap.quantiles
frac.motif.matches.by.category.tib[["bootstrap_ci_upper"]] <- 2 * frac.motif.matches.by.category.tib$num.motifs.per.150bp.matches - lower.bootstrap.quantiles

frac.motif.matches.by.category.plot <- frac.motif.matches.by.category.tib %>% 
  ggplot(aes(x = motif.name, y = num.motifs.per.150bp.matches, fill = peak.category,
             ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar( position = position_dodge(width=0.9), colour="black", width=0.0) +
  xlab("") +
  guides(fill=FALSE) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,  hjust = 1, vjust=0.5))


composite.plot <- devscore.plot.list[[1]] / devscore.plot.list[[2]] / devscore.plot.list[[3]] / motif.d.scores.plot / frac.motif.matches.by.category.plot
ggsave(paste0(outputPlotPrefix, "motif_analysis_composite_plot.svg"), plot = composite.plot, width = 24, height = 18)

### are super-additive peaks more likely to contain dual motifs?
getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks <- function(peak.tib.anno, peak.integration.category, peak.dose) {
  int.categories <- sapply(pull(peak.tib.anno, paste0("peak_integrationCategory-", peak.dose, "-dose")), convertUpregCvalCatToDvalCat)
  int.cat.inds   <- int.categories == peak.integration.category
  filt.peak.tib.anno <- peak.tib.anno[int.cat.inds, ]
  
  dual.motif.tib <- filt.peak.tib.anno %>%
    mutate(hasTGFbMatch = `group-TGFbdominant_maxMotifMatchScore` > 0,
           hasRAMatch   = `group-RAdominant_maxMotifMatchScore` > 0) %>%
    mutate(hasDualMotifMatch = hasTGFbMatch & hasRAMatch)
  
  n.peaks.this.cat <- nrow(dual.motif.tib)
  measured_frac_dual_motif = sum(dual.motif.tib$hasDualMotifMatch) / n.peaks.this.cat
  expected_frac_dual_motif = (sum(dual.motif.tib$hasRAMatch) / n.peaks.this.cat) * (sum(dual.motif.tib$hasTGFbMatch) / n.peaks.this.cat)

  return(list(measured_frac_dual_motif, expected_frac_dual_motif))
}

frac.dual.motif.matches.measured.subadditive   <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeaks, "sub-additive", "med")[[1]]
frac.dual.motif.matches.expected.subadditive   <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeaks, "sub-additive", "med")[[2]]
frac.dual.motif.matches.measured.additive      <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeaks, "additive", "med")[[1]]
frac.dual.motif.matches.expected.additive      <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeaks, "additive", "med")[[2]]
frac.dual.motif.matches.measured.superadditive <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeaks, "super-additive", "med")[[1]]
frac.dual.motif.matches.expected.superadditive <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeaks, "super-additive", "med")[[2]]

# do bootstrap right here, add CI's
set.seed(0)
frac.dual.motif.matches.measured.subadditive.bootstrap.values <- c()
frac.dual.motif.matches.expected.subadditive.bootstrap.values   <- c()
frac.dual.motif.matches.measured.additive.bootstrap.values      <- c()
frac.dual.motif.matches.expected.additive.bootstrap.values      <- c()
frac.dual.motif.matches.measured.superadditive.bootstrap.values <- c()
frac.dual.motif.matches.expected.superadditive.bootstrap.values <- c()
for (ii in 1:n.bootstrap.samples) {
  siUpregPeakBootstrapSample <- sample_n(siUpregPeaks, nrow(siUpregPeaks), replace = TRUE)
  
  bootstrap.frac.dual.motif.matches.measured.subadditive   <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeakBootstrapSample, "sub-additive", "med")[[1]] - frac.dual.motif.matches.measured.subadditive
  bootstrap.frac.dual.motif.matches.expected.subadditive   <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeakBootstrapSample, "sub-additive", "med")[[2]] - frac.dual.motif.matches.expected.subadditive
  bootstrap.frac.dual.motif.matches.measured.additive      <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeakBootstrapSample, "additive", "med")[[1]] - frac.dual.motif.matches.measured.additive
  bootstrap.frac.dual.motif.matches.expected.additive      <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeakBootstrapSample, "additive", "med")[[2]] - frac.dual.motif.matches.expected.additive
  bootstrap.frac.dual.motif.matches.measured.superadditive <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeakBootstrapSample, "super-additive", "med")[[1]] - frac.dual.motif.matches.measured.superadditive
  bootstrap.frac.dual.motif.matches.expected.superadditive <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(siUpregPeakBootstrapSample, "super-additive", "med")[[2]] - frac.dual.motif.matches.expected.superadditive
  
  frac.dual.motif.matches.measured.subadditive.bootstrap.values   <- c(frac.dual.motif.matches.measured.subadditive.bootstrap.values, bootstrap.frac.dual.motif.matches.measured.subadditive)
  frac.dual.motif.matches.expected.subadditive.bootstrap.values   <- c(frac.dual.motif.matches.expected.subadditive.bootstrap.values, bootstrap.frac.dual.motif.matches.expected.subadditive)
  frac.dual.motif.matches.measured.additive.bootstrap.values      <- c(frac.dual.motif.matches.measured.additive.bootstrap.values, bootstrap.frac.dual.motif.matches.measured.additive)
  frac.dual.motif.matches.expected.additive.bootstrap.values      <- c(frac.dual.motif.matches.expected.additive.bootstrap.values, bootstrap.frac.dual.motif.matches.expected.additive)
  frac.dual.motif.matches.measured.superadditive.bootstrap.values <- c(frac.dual.motif.matches.measured.superadditive.bootstrap.values, bootstrap.frac.dual.motif.matches.measured.superadditive)
  frac.dual.motif.matches.expected.superadditive.bootstrap.values <- c(frac.dual.motif.matches.expected.superadditive.bootstrap.values, bootstrap.frac.dual.motif.matches.expected.superadditive)
}
# now build tibble for plot, row-by-row
reduced.dual.motif.analysis.tib <- tibble(expected.vs.measured = c("measured","expected", "measured","expected", "measured", "expected"),
                                          peak.category        = factor(c('sub-additive', 'sub-additive', 'additive', 'additive', 'super-additive', 'super-additive'), levels = c("sub-additive", "additive", "super-additive")),
                                          frac.dual.motifs = c(frac.dual.motif.matches.measured.subadditive, frac.dual.motif.matches.expected.subadditive,
                                                               frac.dual.motif.matches.measured.additive, frac.dual.motif.matches.expected.additive,
                                                               frac.dual.motif.matches.measured.superadditive, frac.dual.motif.matches.expected.superadditive),
                                          upper.ci = frac.dual.motifs - c(quantile(frac.dual.motif.matches.measured.subadditive.bootstrap.values, .05), quantile(frac.dual.motif.matches.expected.subadditive.bootstrap.values, .05),
                                                                          quantile(frac.dual.motif.matches.measured.additive.bootstrap.values, .05), quantile(frac.dual.motif.matches.expected.additive.bootstrap.values, .05),
                                                                          quantile(frac.dual.motif.matches.measured.superadditive.bootstrap.values, .05), quantile(frac.dual.motif.matches.expected.superadditive.bootstrap.values, .05)),
                                          lower.ci = frac.dual.motifs - c(quantile(frac.dual.motif.matches.measured.subadditive.bootstrap.values, .95), quantile(frac.dual.motif.matches.expected.subadditive.bootstrap.values, .95),
                                                                          quantile(frac.dual.motif.matches.measured.additive.bootstrap.values, .95), quantile(frac.dual.motif.matches.expected.additive.bootstrap.values, .95),
                                                                          quantile(frac.dual.motif.matches.measured.superadditive.bootstrap.values, .95), quantile(frac.dual.motif.matches.expected.superadditive.bootstrap.values, .95)),
                                          )

dual.motif.analysis.plot <- reduced.dual.motif.analysis.tib %>%
  ggplot(aes(x = peak.category, y = frac.dual.motifs, fill = expected.vs.measured, ymin = lower.ci, ymax = upper.ci)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(position = position_dodge(width=0.9), width = 0) +
  theme_classic()

ggsave(paste0(outputPlotPrefix, "motif_analysis_freq_dual_motif_matches_by_peakIntCategory.svg"), plot = dual.motif.analysis.plot, width = 12, height = 12)
