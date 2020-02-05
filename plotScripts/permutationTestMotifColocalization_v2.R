library(tidyverse)
library(here)

# joinedTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv'))
# geneTib   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTib   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))

n.permutations.to.make.null.distributions <- 10000
superadditive.peak.threshold.orig.expression.units.above.additive.prediction <- 1.5

peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > superadditive.peak.threshold.orig.expression.units.above.additive.prediction) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()

permTestMotifEnrichment <- function(peaksUp, motif_column_name, n.permutations.to.make.null.distributions) {
  peaksUp.nonSuperAdditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
  peaksUp.superAdditive     <- filter(peaksUp,   `isSuperAdditivePeak-med`)
  
  motifMatchesNotsuperadditive <- pull(peaksUp.nonSuperAdditive, motif_column_name) > 1
  motifMatchesSuperadditive    <- pull(peaksUp.superAdditive, motif_column_name) > 1
  motifMatchesAll              <- pull(peaksUp, motif_column_name) > 1
  
  nonsuperadditive.average <- sum(motifMatchesNotsuperadditive) / length(motifMatchesNotsuperadditive)
  superadditive.average    <- sum(motifMatchesSuperadditive) / length(motifMatchesSuperadditive)
  print(sprintf("fraction single match, non-superadditive: %0.3f", nonsuperadditive.average))
  print(sprintf("fraction single match, superadditive: %0.3f", superadditive.average))
  
  samplesize <- length(motifMatchesSuperadditive)
  resultsvec <- permutationTest(samplesize, motifMatchesAll, n.permutations.to.make.null.distributions)
  f <- ecdf(resultsvec)
  percentile_superadditive_peaks <- f(superadditive.average)

  motif_name <- motif_column_name
  p <- qplot(resultsvec, bins = 500) +
    geom_vline(xintercept = superadditive.average, color = "red") +
    theme_minimal(base_size = 16) +
    ggtitle(sprintf("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have %s motif matches", motif_name)) +
    xlab(sprintf("Fraction of peaks with a %s motif match \n (red line is %.3f quantile)", motif_name, percentile_superadditive_peaks))

  return(list(p, percentile_superadditive_peaks))
}

permTestDualMotifEnrichment <- function(peaksUp, motif_column_name1, motif_column_name2, n.permutations.to.make.null.distributions) {
  peaksUp.nonSuperAdditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
  peaksUp.superAdditive     <- filter(peaksUp,   `isSuperAdditivePeak-med`)
  
  motifMatchesNotsuperadditive <- (pull(peaksUp.nonSuperAdditive, motif_column_name1) > 1) & (pull(peaksUp.nonSuperAdditive, motif_column_name2) > 1)
  motifMatchesSuperadditive    <- (pull(peaksUp.superAdditive, motif_column_name1) > 1) & (pull(peaksUp.superAdditive, motif_column_name2) > 1)
  motifMatchesAll              <- (pull(peaksUp, motif_column_name1) > 1) & (pull(peaksUp, motif_column_name2) > 1)
  
  nonsuperadditive.average <- sum(motifMatchesNotsuperadditive) / length(motifMatchesNotsuperadditive)
  superadditive.average    <- sum(motifMatchesSuperadditive) / length(motifMatchesSuperadditive)
  print(sprintf("fraction dual match, non-superadditive: %0.3f", nonsuperadditive.average))
  print(sprintf("fraction dual match, superadditive: %0.3f", superadditive.average))
  
  samplesize <- length(motifMatchesSuperadditive)
  resultsvec <- permutationTest(samplesize, motifMatchesAll, n.permutations.to.make.null.distributions)
  f <- ecdf(resultsvec)
  percentile_superadditive_peaks <- f(superadditive.average)
  
  motif_name1 <- motif_column_name1
  motif_name2 <- motif_column_name2
  
  p <- qplot(resultsvec, bins = 500) +
    geom_vline(xintercept = superadditive.average, color = "red") +
    theme_minimal(base_size = 12) +
    ggtitle(sprintf("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have %s and %s motif matches", motif_name1, motif_name2)) +
    xlab(sprintf("Fraction of peaks with a dual motif group match \n (red line is %.3f quantile)", percentile_superadditive_peaks))
  
  return(list(p, percentile_superadditive_peaks))
}

permutationTest <- function(samplesize, motifMatchLogicalVector, n.permutations.to.make.null.distributions) {
  #samplesize <- length(dualMotifMatchesSuperadditive)
  set.seed(0)
  resultsvec <- c()
  for (ii in 1:n.permutations.to.make.null.distributions) {
    this.sample <- sample(motifMatchLogicalVector, samplesize)
    fracDual <- sum(this.sample) / samplesize
    resultsvec <- c(resultsvec, fracDual)
  }
  return(resultsvec)
}

loginds_indMotifs <- grepl("motifMatchScore", colnames(peaksUp)) & !(grepl("group", colnames(peaksUp)))
indMotifNames <- colnames(peaksUp)[loginds_indMotifs]
for (motif_name in indMotifNames) {
  res <- permTestMotifEnrichment(peaksUp, motif_name, n.permutations.to.make.null.distributions)
  output_fn <- here("plots", "permutationTestMotifPlots_SuperAdditivePeaks", sprintf("%s_singleMotifPermTestUpregPeaks.pdf", motif_name))
  ggsave(output_fn, res[[1]])
}


loginds_motifGroups <- grepl("maxMotifMatchScore", colnames(peaksUp)) & (grepl("group", colnames(peaksUp)))
motifGroupNames <- colnames(peaksUp)[loginds_motifGroups]
for (motif_name in motifGroupNames) {
  res <- permTestMotifEnrichment(peaksUp, motif_name, n.permutations.to.make.null.distributions)
  output_fn <- here("plots", "permutationTestGroupMotifPlots_SuperAdditivePeaks", sprintf("%s_singleMotifPermTestUpregPeaks.pdf", motif_name))
  ggsave(output_fn, res[[1]])
}


loginds_motifGroups <- grepl("maxMotifMatchScore", colnames(peaksUp)) & (grepl("group", colnames(peaksUp)))
motifGroupNames <- colnames(peaksUp)[loginds_motifGroups]
nGroups <- length(motifGroupNames)
# for (ii in 1:(nGroups-1)) {
#   for (jj in (ii+1):nGroups) {
for (ii in 1:nGroups) {
  for (jj in 1:nGroups) {
    motifGroupName1 <- motifGroupNames[ii]
    motifGroupName2 <- motifGroupNames[jj]
    res <- permTestDualMotifEnrichment(peaksUp, motifGroupName1, motifGroupName2, n.permutations.to.make.null.distributions)
    output_fn <- here("plots", "permutationTestDualMotifPlots_SuperAdditivePeaks", sprintf("%s_%s_dualMotifPermTestUpregPeaks.pdf", motifGroupName1, motifGroupName2))
    ggsave(output_fn, res[[1]])
  }
}


