library(tidyverse)
# 
# allpeaks   <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.tsv")
# joinedupregpeaks <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/upregJoinedPeakGeneTib_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.tsv")
# upreggenes       <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
# 

upregpeaks <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")

use.d.scores.within.this.quantile.range <- 0.00
lower_d_filter <- quantile(upregpeaks$`peakAdditivePredFcResidual-med`, use.d.scores.within.this.quantile.range)
upper_d_filter <- quantile(upregpeaks$`peakAdditivePredFcResidual-med`, 1 - use.d.scores.within.this.quantile.range)
filtupregpeaks <- upregpeaks %>% filter(`peakAdditivePredFcResidual-med` > lower_d_filter, `peakAdditivePredFcResidual-med` < upper_d_filter,
                                        `peak_integrationCategory-med-dose` %in% superadditive.peak.categories)

motif.col.names     <- colnames(filtupregpeaks)[which(grepl("_motifMatchScore", colnames(filtupregpeaks)))]
motif.names <- sapply(strsplit(motif.col.names, "_"), function(x) x[[1]])
corresponding.avg.d.scores <- c()
corresponding.median.d.scores <- c()
peak.d.values <- pull(filtupregpeaks, "peakAdditivePredFcResidual-med")
n.peaks <- length(peak.d.values)

# make histograms of d values of motifs across all upregulated peaks
indmotifplots <- list()
counter <- 1
for (motif.col.name in motif.col.names) {
  motif.match.indices <- which(filtupregpeaks[, motif.col.name] > 0)
  frac.peaks.with.motif <- length(motif.match.indices) / n.peaks
  print(sprintf("%s, in %0.3f percent of %d motifs", motif.col.name, frac.peaks.with.motif, n.peaks))
  this.avg.d          <- mean(peak.d.values[motif.match.indices])
  this.median.d       <- median(peak.d.values[motif.match.indices])
  corresponding.avg.d.scores    <- c(corresponding.avg.d.scores, this.avg.d)
  corresponding.median.d.scores <- c(corresponding.median.d.scores, this.median.d)
  p <- qplot(peak.d.values[motif.match.indices]) + ggtitle(strsplit(motif.col.name, "_")[[1]][1]) + ylab("counts")
  indmotifplots[[counter]] <- p
  counter <- counter + 1
  # print(p)
}

tfp <- tibble(motif_name = motif.names, avg_d_val = corresponding.avg.d.scores, median_d_val = corresponding.median.d.scores)

qplot(upregpeaks$`peakAdditivePredFcResidual-med`, bins = 100) + xlim(-5, 5) + ggtitle("d-value distribution, upregulated peaks, medium dose") + xlab("d-value")

ggplot(tfp, aes(x = motif.names, y = corresponding.avg.d.scores)) + 
  theme_minimal(base_size = 10) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5)) + 
  geom_hline(yintercept = mean(filtupregpeaks$`peakAdditivePredFcResidual-med`))

ggplot(tfp, aes(x = motif.names, y = median_d_val)) + 
  theme_minimal(base_size = 10) +
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5)) + 
  geom_hline(yintercept = median(filtupregpeaks$`peakAdditivePredFcResidual-med`))


# look at frequency of motif for sub-additive, additive, and super-additive peaks

superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")
additive.peak.categories      <- c("additive", "ambiguous")
subadditive.peak.categories   <- c("sub-additive")


npeaks.all <- nrow(upregpeaks)
res.tib    <- NULL
for (motif.col.name in motif.col.names) {
  subadd.peaks <- upregpeaks %>% filter(`peak_integrationCategory-med-dose` %in% subadditive.peak.categories)
  n.subadd.peaks <- nrow(subadd.peaks)
  subadd.motif.match.indices <- which(subadd.peaks[, motif.col.name] > 0)
  subadd.frac.peaks.with.motif <- length(subadd.motif.match.indices) / n.subadd.peaks
  
  add.peaks <- upregpeaks %>% filter(`peak_integrationCategory-med-dose` %in% additive.peak.categories)
  n.add.peaks <- nrow(add.peaks)
  add.motif.match.indices <- which(add.peaks[, motif.col.name] > 0)
  add.frac.peaks.with.motif <- length(add.motif.match.indices) / n.add.peaks
  
  superadd.peaks <- upregpeaks %>% filter(`peak_integrationCategory-med-dose` %in% superadditive.peak.categories)
  n.superadd.peaks <- nrow(superadd.peaks)
  superadd.motif.match.indices <- which(superadd.peaks[, motif.col.name] > 0)
  superadd.frac.peaks.with.motif <- length(superadd.motif.match.indices) / n.superadd.peaks
  
  motif.name   <- rep(strsplit(motif.col.name, "_")[[1]][1], 3)
  frac.peaks.with.motif.matches <- c(subadd.frac.peaks.with.motif, add.frac.peaks.with.motif, superadd.frac.peaks.with.motif)
  peak.category <- factor(c('subadditive', 'additive', 'superadditive'), levels = c('subadditive', 'additive', 'superadditive'))
  peak.category.total.num.peaks <- c(n.subadd.peaks, n.add.peaks, n.superadd.peaks)
  
  res.tib <- rbind(res.tib, tibble(motif.name, frac.peaks.with.motif.matches, peak.category, peak.category.total.num.peaks))
}


res.tib %>% ggplot(aes(x = motif.name, y = frac.peaks.with.motif.matches, fill = peak.category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust=0.5))













