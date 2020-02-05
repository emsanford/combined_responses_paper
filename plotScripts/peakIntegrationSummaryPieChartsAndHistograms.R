library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  # siUpregPeaks   <- read_tsv(here('plots', 'sensitivity_analysis',  'differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.tsv'))
  outputPrefix  <- here("plots", paste0("stackedBarHist_upregPeakIntegrationConstants_"))
} else {
  siUpregPeaks  <- read_tsv(cmdargs[1])
  outputPrefix <- cmdargs[2]
}



# svg(filename=here("plots", "rTpmBeeSwarmsAll", sprintf("piechart_%s_%s.svg", signal.effect.direction, "0-low")),width=8,height=8)
pie(table(siUpregPeaks$`peak_integrationCategory-low-dose`))
pie(table(siUpregPeaks$`peak_integrationCategory-med-dose`))
pie(table(siUpregPeaks$`peak_integrationCategory-high-dose`))
siUpregPeaks$`peak_integrationConstant-low` %>% qplot() + xlim(-2, 2) + ylim(0, 155)
siUpregPeaks$`peak_integrationConstant-med` %>% qplot() + xlim(-2, 2) + ylim(0, 155)
siUpregPeaks$`peak_integrationConstant-high` %>% qplot() + xlim(-2, 2) + ylim(0, 155)
# dev.off()


####### make stacked histogram plot showing signal integration constant & colored by frequency of each category in the plot

# function: reassign weird categories that may come up when a different dose integrates in a different direction to "uncategorized" 
mapCatsToReducesCatSet <- function(cat.values) {
  allowedCategories <- c("ambiguous", "sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative")
  cat.values[! cat.values %in% allowedCategories] <- "uncategorized"
  return(factor(cat.values, levels = rev(c("uncategorized", allowedCategories))))
}

# function: assign values to a bin
# input -- a vector of integration contstants
# output -- bin values for the vector. order of values in the vector doesn't change
#   note: values outside the range get assigned to the lowest or highest bin available
assignValuesToHistBin <- function(values, bin.midpoints, bin.radius) {
  n.vals <- length(values)
  outputVec <- c()
  for (ii in 1:n.vals) {
    this.val <- values[ii]
    distvec <- abs(this.val - bin.midpoints)
    lowest.bin.distance <- min(distvec)
    bin.index <- which(distvec == lowest.bin.distance)[1]
    outputVec <- c(outputVec, bin.midpoints[bin.index])
  }
  return(outputVec)
}


for (dosage in c("low", "med", "high")) {
  categorical.values <- pull(siUpregPeaks, paste0("peak_integrationCategory-", dosage ,"-dose"))
  hist.values        <- pull(siUpregPeaks, paste0("peak_integrationConstant-", dosage))
  
  bin.midpoints <- seq(-4, 8, by = 0.1)
  bin.radius    <- (bin.midpoints[2] - bin.midpoints[1]) / 2
  
  mapped.categorical.values <- mapCatsToReducesCatSet(categorical.values)
  bin.values <- assignValuesToHistBin(hist.values, bin.midpoints, bin.radius)
  
  stackedBarHistTib <- tibble(intConstantHhistBin = bin.values, intCategory = mapped.categorical.values)
  
  stackedBarHist <- ggplot(stackedBarHistTib, aes(x = intConstantHhistBin, fill = intCategory)) +
    geom_bar(stat="count", width = bin.radius * 2 * .90) + 
    theme_minimal(base_size = 12) + 
    xlab("integration constant value for a peak") +
    ylab("number of peaks") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 1) +
    ggtitle(paste0("Distribution of integration constants for upregulated peaks\n", dosage, " dose"))
  
  ggsave(paste0(outputPrefix, "_", dosage, "_dose.svg"), width = 12, height = 4, plot = stackedBarHist)
}

# "d" value histograms for peaks
d1 <- siUpregPeaks$`peakAdditivePredFcResidual-low` %>% qplot(bins = 50) + xlim(-5, 5) + ylab("number of upregulated peaks") + xlab("fold-change difference from additive prediction") + theme_minimal(base_size = 16)
d2 <- siUpregPeaks$`peakAdditivePredFcResidual-med` %>% qplot(bins = 50) + xlim(-5, 5) + ylab("number of upregulated peaks") + xlab("fold-change difference from additive prediction") + theme_minimal(base_size = 16)
d3 <- siUpregPeaks$`peakAdditivePredFcResidual-high` %>% qplot(bins = 50) + xlim(-5, 5) + ylab("number of upregulated peaks") + xlab("fold-change difference from additive prediction") + theme_minimal(base_size = 16)

d4 <- siUpregPeaks$`peakMultiplicativePredFcResidual-low` %>% qplot(bins = 50) + xlim(-5, 5) + ylab("number of upregulated peaks") + xlab("fold-change difference from multiplicative prediction") + theme_minimal(base_size = 16)
d5 <- siUpregPeaks$`peakMultiplicativePredFcResidual-med` %>% qplot(bins = 50) + xlim(-5, 5) + ylab("number of upregulated peaks") + xlab("fold-change difference from multiplicative prediction") + theme_minimal(base_size = 16)
d6 <- siUpregPeaks$`peakMultiplicativePredFcResidual-high` %>% qplot(bins = 50) + xlim(-5, 5) + ylab("number of upregulated peaks") + xlab("fold-change difference from multiplicative prediction") + theme_minimal(base_size = 16)

