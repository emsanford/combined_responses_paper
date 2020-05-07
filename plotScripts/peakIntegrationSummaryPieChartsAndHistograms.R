library(tidyverse)
library(here)
library(patchwork)

source(here('extractionScripts', 'util.R'))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregPeaks   <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  output.folder  <- here("plots", "peak_integration_summary_plots")
} else {
  siUpregPeaks  <- read_tsv(cmdargs[1])
  output.folder <- cmdargs[2]
}

piechart.location.prefix <- paste0(output.folder, '/peak_integration_category_pie_chart_')
factor.order.peak.categories <- c("sub-additive", "additive", "super-additive")
n.upreg.peaks <- nrow(siUpregPeaks)

svg(filename=paste0(piechart.location.prefix, "low_dose_", "n", n.upreg.peaks, ".svg"),width=8,height=8)
input.vector <- factor(sapply(siUpregPeaks$`peak_integrationCategory-low-dose`, convertUpregCvalCatToDvalCat), levels = factor.order.peak.categories)
pieplot.low  <- pie(table(input.vector), main = paste0("low dose, N =", n.upreg.peaks), init.angle = 150, clockwise = TRUE)
print(pieplot.low)
dev.off()

svg(filename=paste0(piechart.location.prefix, "med_dose_", "n", n.upreg.peaks, ".svg"),width=8,height=8)
input.vector <- factor(sapply(siUpregPeaks$`peak_integrationCategory-med-dose`, convertUpregCvalCatToDvalCat), levels = factor.order.peak.categories)
pieplot.med  <- pie(table(input.vector), main = paste0("med dose, N =", n.upreg.peaks), init.angle = 150, clockwise = TRUE)
print(pieplot.med)
dev.off()

svg(filename=paste0(piechart.location.prefix, "high_dose_", "n", n.upreg.peaks, ".svg"),width=8,height=8)
input.vector <- factor(sapply(siUpregPeaks$`peak_integrationCategory-high-dose`, convertUpregCvalCatToDvalCat), levels = factor.order.peak.categories)
pieplot.high  <- pie(table(input.vector), main = paste0("high dose, N =", n.upreg.peaks), init.angle = 150, clockwise = TRUE)
print(pieplot.high)
dev.off()


####### make stacked histogram plot showing signal integration constant & colored by frequency of each category in the plot

# function: reassign weird categories that may come up when a different dose integrates in a different direction to "uncategorized" 
mapCatsToReducesCatSet <- function(cat.values) {
  allowedCategories <- c("ambiguous", "sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative")
  cat.values[! cat.values %in% allowedCategories] <- "uncategorized"
  return(factor(cat.values, levels = rev(c("uncategorized", allowedCategories))))
}


# make stacked bar histograms for D value
bin.step.size   <-  0.05
bin.leftmost    <- -2
bin.rightmost   <-  2
plot.width      <-  16
plot.height     <-  8
bin.radius      <- bin.step.size / 2
bin.midpoints   <- seq(bin.leftmost + bin.step.size, bin.rightmost, by = bin.step.size) - bin.radius
use.common.scale <- T
standard.ylim <- 500
threshold.for.reporting.upper.end.of.histogram <- 1.5

output.fig.list <- list()
counter <- 1
for (dosage in c("low", "med", "high")) {
  categorical.values <- sapply(pull(siUpregPeaks, paste0("peak_integrationCategory-", dosage ,"-dose")), convertUpregCvalCatToDvalCat)
  hist.values        <- pull(siUpregPeaks, paste0("peakAdditivePredFcResidual-", dosage))
  
  stackedBarRes <- makeHistogramOfValues(hist.values, categorical.values, bin.leftmost, bin.rightmost,
                                         bin.step.size, paste0(dosage, " dose, d-values"), 
                                         xlabel = "d-value", ylabel = "count", color.by.category = T)
  
  stackedBarHist <- stackedBarRes[[1]]
  
  if (use.common.scale) {
    stackedBarHist <- stackedBarHist + ylim(0, standard.ylim)
  }
  
  stackedBarTib  <- stackedBarRes[[2]]
  stackedBarTib[["dose"]] <- dosage
  reduced.tib <- stackedBarTib %>% 
    group_by(intConstantHhistBin, intCategory, dose) %>% 
    mutate(n_this_bin = n(), freq_this_bin = n_this_bin / nrow(stackedBarTib)) %>%
    ungroup() %>%
    unique()
  
  freq.above.c.2 <- reduced.tib %>%
    filter(intConstantHhistBin >= threshold.for.reporting.upper.end.of.histogram) %>%
    pull("freq_this_bin") %>%
    sum()
  
  print(sprintf("freq d.val above %0.2f: (%s dose) %0.3f", threshold.for.reporting.upper.end.of.histogram, dosage, freq.above.c.2))
  
  ggsave(paste0(output.folder, "/dval_addPredDiff_histogram_", dosage, "_dose.svg"), width = plot.width, height = plot.height, plot = stackedBarHist)
  
  output.fig.list[[counter]] <- stackedBarHist
  counter <- counter + 1
}

grand.plot.dvals <- output.fig.list[[1]] + output.fig.list[[2]] + output.fig.list[[3]] 
ggsave(paste0(output.folder, "/dval_addPredDiff_composed_histogram.svg"), width = plot.width * 3, height = plot.height, plot = grand.plot.dvals) 

## uncomment this block of code to see similar plots but for the fold-change difference from the multiplicative prediction
for (dosage in c("low", "med", "high")) {
  categorical.values <- sapply(pull(siUpregPeaks, paste0("peak_integrationCategory-", dosage ,"-dose")), convertUpregCvalCatToDvalCat)
  hist.values        <- pull(siUpregPeaks, paste0("peakMultiplicativePredFcResidual-", dosage))

  stackedBarRes <- makeHistogramOfValues(hist.values, categorical.values, bin.leftmost, bin.rightmost,
                                         bin.step.size, paste0(dosage, " dose, d-values"),
                                         xlabel = "d-value", ylabel = "count", color.by.category = T)

  stackedBarHist <- stackedBarRes[[1]]

  if (use.common.scale) {
    stackedBarHist <- stackedBarHist + ylim(0, standard.ylim)
  }

  stackedBarTib  <- stackedBarRes[[2]]
  stackedBarTib[["dose"]] <- dosage
  reduced.tib <- stackedBarTib %>%
    group_by(intConstantHhistBin, intCategory, dose) %>%
    mutate(n_this_bin = n(), freq_this_bin = n_this_bin / nrow(stackedBarTib)) %>%
    ungroup() %>%
    unique()

  freq.above.c.2 <- reduced.tib %>%
    filter(intConstantHhistBin >= threshold.for.reporting.upper.end.of.histogram) %>%
    pull("freq_this_bin") %>%
    sum()

  print(sprintf("%s %0.3f", dosage, freq.above.c.2))

  ggsave(paste0(output.folder, "/dval_multPredDiff_histogram_", dosage, "_dose.svg"), width = plot.width, height = plot.height, plot = stackedBarHist)

  output.fig.list[[counter]] <- stackedBarHist
  counter <- counter + 1
}

grand.plot.dvals <- output.fig.list[[1]] + output.fig.list[[2]] + output.fig.list[[3]] 
ggsave(paste0(output.folder, "/dval_multPredDiff_composed_histogram.svg"), width = plot.width * 3, height = plot.height, plot = grand.plot.dvals) 
