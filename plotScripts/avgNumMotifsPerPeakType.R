library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)
library(chromVARmotifs)
library(BiocParallel)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
source(here('extractionScripts', 'util.R'))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  upreg.peaks          <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  mostVariableMotifSet <- read_rds("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/mostVariableMotifs_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
  outputPlotPrefix     <- here("plots", "")
} else {
  upreg.peaks           <- read_tsv(cmdargs[1])
  mostVariableMotifSet <- read_rds(cmdargs[2])
  outputPlotPrefix     <- cmdargs[3]
}

n.bootstrap.samples <- 1000

#0 define motif set (top 50 or all in cisBP data set)
data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
cisbp_motifs <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper

calcAvgNumMotifsPerPeakType <- function(upreg.peaks.granges, motif_ix, motifsetname) {
  restib <- NULL
  
  # 3: count number of motif matches 
  motif.count.matx <- as.matrix(motifCounts(motif_ix))
  motif.sum.each.peak <- rowSums(motif.count.matx)
  peak.types <- names(motif.sum.each.peak)
  
  inds.subadditive   <- peak.types == "sub-additive"
  inds.additive      <- peak.types == "additive"
  inds.superadditive <- peak.types == "super-additive"
  
  peak_sizes <- width(ranges(upreg.peaks.granges))
  subadd.total.peak.size     <- sum(peak_sizes[inds.subadditive])  
  add.total.peak.size        <- sum(peak_sizes[inds.additive])  
  superadd.total.peak.size   <- sum(peak_sizes[inds.superadditive]) 
  
  subadd.avg.peak.size     <- sum(peak_sizes[inds.subadditive]) / sum(inds.subadditive)
  add.avg.peak.size        <- sum(peak_sizes[inds.additive]) / sum(inds.additive)  
  superadd.avg.peak.size   <- sum(peak_sizes[inds.superadditive]) / sum(inds.superadditive)  
  
  subadd.num.matches.per.150.peak.bp   <- 150 * sum(motif.sum.each.peak[inds.subadditive])   / subadd.total.peak.size
  add.num.matches.per.150.peak.bp      <- 150 * sum(motif.sum.each.peak[inds.additive])      / add.total.peak.size
  superadd.num.matches.per.150.peak.bp <- 150 * sum(motif.sum.each.peak[inds.superadditive]) / superadd.total.peak.size
  
  # 4: store results in result tibble
  restib <- tibble(avg_motifs_per_peak = c(subadd.num.matches.per.150.peak.bp, add.num.matches.per.150.peak.bp, superadd.num.matches.per.150.peak.bp), 
                   avg_peak_size       = c(subadd.avg.peak.size, add.avg.peak.size, superadd.avg.peak.size),
                   peak_type = factor(c("sub-additive", "additive", "super-additive"), levels = c("sub-additive", "additive", "super-additive")), 
                   motif_set = rep(motifsetname, 3))

  return(restib)
}

# 1: make genomic ranges for upreg peaks (medium dose)

peaktypes        <- sapply(upreg.peaks$`peak_integrationCategory-med-dose`, convertUpregCvalCatToDvalCat)
names(peaktypes) <- peaktypes

upreg.peaks.granges <- GRanges(seqnames = upreg.peaks$chrom,
                               ranges = IRanges(start = upreg.peaks$startLocs,
                                                end   = upreg.peaks$endLocs),
                               peaktype = peaktypes)

# 2: match motifs to these genomic ranges (may need hg38 reference seq). pick full motif set or top50 motifs
motif_ix_mostVariable <- matchMotifs(mostVariableMotifSet, upreg.peaks.granges, genome = BSgenome.Hsapiens.UCSC.hg38, out = "scores")
motif_ix_allCisBP     <- matchMotifs(cisbp_motifs, upreg.peaks.granges, genome = BSgenome.Hsapiens.UCSC.hg38, out = "scores")

result.tib.topmotifs <- calcAvgNumMotifsPerPeakType(upreg.peaks.granges, motif_ix_mostVariable, "most variable motifs")
result.tib.allmotifs <- calcAvgNumMotifsPerPeakType(upreg.peaks.granges, motif_ix_allCisBP, "all cisBP motifs")

 
# add bootstrap confidence intervals
set.seed(0)
n.peaks <- nrow(upreg.peaks)
bootstrap.res.tib <- NULL
for (ii in 1:n.bootstrap.samples) {
  this.sample.inds <- sample(1:n.peaks, n.peaks, replace = TRUE)

  this.result.tib.topmotifs <- calcAvgNumMotifsPerPeakType(upreg.peaks.granges[this.sample.inds], motif_ix_mostVariable[this.sample.inds, ], "most variable motifs")
  this.result.tib.allmotifs <- calcAvgNumMotifsPerPeakType(upreg.peaks.granges[this.sample.inds], motif_ix_allCisBP[this.sample.inds, ], "all cisBP motifs")

  bootstrap.res.tib <- rbind(bootstrap.res.tib, this.result.tib.topmotifs, this.result.tib.allmotifs)
}



bootstrap.res.topmotifs <- bootstrap.res.tib %>% filter(motif_set == "most variable motifs")
bootstrap.res.allmotifs <- bootstrap.res.tib %>% filter(motif_set == "all cisBP motifs")

for (col.of.interest in c("avg_motifs_per_peak", "avg_peak_size")) {
  counter <- 1
  for (peak.type.of.interest in c("sub-additive", "additive", "super-additive")) {
    avg.value        <- result.tib.topmotifs    %>% filter(peak_type == peak.type.of.interest) %>% pull(col.of.interest)
    bootstrap.values <- bootstrap.res.topmotifs %>% filter(peak_type == peak.type.of.interest) %>% pull(col.of.interest)
    result.tib.topmotifs[[paste0(col.of.interest, "_bootstrap_ci_lower")]][counter] <- 2 * avg.value - quantile(bootstrap.values, 0.95)
    result.tib.topmotifs[[paste0(col.of.interest, "_bootstrap_ci_upper")]][counter] <- 2 * avg.value - quantile(bootstrap.values, 0.05)
    
    avg.value        <- result.tib.allmotifs    %>% filter(peak_type == peak.type.of.interest) %>% pull(col.of.interest)
    bootstrap.values <- bootstrap.res.allmotifs %>% filter(peak_type == peak.type.of.interest) %>% pull(col.of.interest)
    result.tib.allmotifs[[paste0(col.of.interest, "_bootstrap_ci_lower")]][counter] <- 2 * avg.value - quantile(bootstrap.values, 0.95)
    result.tib.allmotifs[[paste0(col.of.interest, "_bootstrap_ci_upper")]][counter] <- 2 * avg.value - quantile(bootstrap.values, 0.05)
    
    counter <- counter + 1
  }
}

result.tib.topmotifs
result.tib.allmotifs

p1 <- ggplot(result.tib.topmotifs, aes(y = avg_motifs_per_peak, x = peak_type, ymin = avg_motifs_per_peak_bootstrap_ci_lower, ymax = avg_motifs_per_peak_bootstrap_ci_upper)) +
  geom_bar(stat = "identity") + ggtitle( "Enriched Motif Set") + 
  geom_errorbar(position = position_dodge(width=0.9), width = 0) +
  ylab("average number of motif matches per 150 bp of sequence") +
  theme_classic()

p2 <- ggplot(result.tib.allmotifs, aes(y = avg_motifs_per_peak, x = peak_type,  ymin = avg_motifs_per_peak_bootstrap_ci_lower, ymax = avg_motifs_per_peak_bootstrap_ci_upper)) +
  geom_bar(stat = "identity") + ggtitle("All cisBP Motifs") + 
  geom_errorbar(position = position_dodge(width=0.9), width = 0) +
  ylab("average number of motif matches per 150 bp of sequence") +
  theme_classic() 

p3 <- ggplot(result.tib.allmotifs, aes(y = avg_peak_size, x = peak_type,  ymin = avg_peak_size_bootstrap_ci_lower, ymax = avg_peak_size_bootstrap_ci_upper)) +
   geom_bar(stat = "identity") + ggtitle("peak size") + 
   geom_errorbar(position = position_dodge(width=0.9), width = 0) +
   ylab("average peak width (bp)") +
   theme_classic()
  
# library(patchwork)
# p3 + p1 + p2

ggsave(paste0(outputPlotPrefix, "motifByPeakTypeDensityEnrichedMotifSet.svg"), plot = p1, width = 12, height = 12)
ggsave(paste0(outputPlotPrefix, "motifByPeakTypeDensityFullMotifSet.svg"), plot = p2, width = 12, height = 12)
ggsave(paste0(outputPlotPrefix, "avgPeakWidthByPeakType.svg"), plot = p3, width = 12, height = 12)
