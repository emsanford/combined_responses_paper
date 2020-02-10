library(tidyverse)
# 
# allpeaks   <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.tsv")
# joinedupregpeaks <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/upregJoinedPeakGeneTib_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.tsv")
# upreggenes       <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
# 

upregpeaks <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")

use.d.scores.within.this.quantile.range <- 0.01
lower_d_filter <- quantile(upregpeaks$`peakAdditivePredFcResidual-med`, use.d.scores.within.this.quantile.range)
upper_d_filter <- quantile(upregpeaks$`peakAdditivePredFcResidual-med`, 1 - use.d.scores.within.this.quantile.range)
filtupregpeaks <- upregpeaks %>% filter(`peakAdditivePredFcResidual-med` > lower_d_filter, `peakAdditivePredFcResidual-med` < upper_d_filter)

motif.col.names     <- colnames(filtupregpeaks)[which(grepl("_motifMatchScore", colnames(filtupregpeaks)))]
motif.names <- sapply(strsplit(motif.col.names, "_"), function(x) x[[1]])
corresponding.avg.d.scores <- c()
corresponding.median.d.scores <- c()
peak.d.values <- pull(filtupregpeaks, "peakAdditivePredFcResidual-med")
n.peaks <- length(peak.d.values)
indmotifplots <- list()
counter <- 1
for (motif.col.name in motif.col.names) {
  motif.match.indices <- which(filtupregpeaks[, motif.col.name] > 0)
  print(sprintf("%s, in %0.3f percent of %d motifs", motif.col.name, length(motif.match.indices) / n.peaks, n.peaks))
  this.avg.d          <- mean(peak.d.values[motif.match.indices])
  this.median.d       <- median(peak.d.values[motif.match.indices])
  corresponding.avg.d.scores    <- c(corresponding.avg.d.scores, this.avg.d)
  corresponding.median.d.scores <- c(corresponding.median.d.scores, this.median.d)
  p <- qplot(peak.d.values[motif.match.indices]) + ggtitle(strsplit(motif.col.name, "_")[[1]][1]) + ylab("counts")
  indmotifplots[[counter]] <- p
  counter <- counter + 1
  print(p)
}

print(indmotifplots[1])

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


