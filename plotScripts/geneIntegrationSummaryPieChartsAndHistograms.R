library(tidyverse)
library(here)

source(here('extractionScripts', 'util.R'))

use.common.scale <- F

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes     <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
  addPredFcDiffMin <- 0
  minTpmDiff       <- 0
  output.folder    <- here('plots', 'gene_integration_summary_plots')
} else {
  siUpregGenes     <- read_tsv(cmdargs[1])
  addPredFcDiffMin <- as.numeric(cmdargs[2])
  minTpmDiff       <- as.numeric(cmdargs[3])
  output.folder    <- cmdargs[4]
}

filtSiUpregGenes <- siUpregGenes %>% 
  filter(`addMultPredFcDiff-low` >= addPredFcDiffMin,
         (`multPred-low` - `addPred-low`) >= minTpmDiff,
         `addMultPredFcDiff-med` >= addPredFcDiffMin,
         (`multPred-med` - `addPred-med`) >= minTpmDiff,
         `addMultPredFcDiff-high` >= addPredFcDiffMin,
         (`multPred-high` - `addPred-high`) >= minTpmDiff)

piechart.location.prefix <- paste0(output.folder, '/gene_integration_mode_pie_chart_')
stackedBarHistogram.location.prefix <- paste0(output.folder, '/gene_integration_mode_stackedBarHistogram_')

n_upreg_genes <- nrow(siUpregGenes)
####### first, make pie charts for the categorical description for each upregulated gene
factor.order.gene.categories <- c("sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative", "ambiguous")

svg(filename=paste0(piechart.location.prefix, "low_dose_", "n", n_upreg_genes, ".svg"),width=8,height=8)
pieplot.low  <- pie(table(factor(siUpregGenes$`integrationCategory-low-dose`, levels = factor.order.gene.categories)), main = paste0("low dose, N =", n_upreg_genes))
print(pieplot.low)
dev.off()

svg(filename=paste0(piechart.location.prefix, "med_dose_", "n", n_upreg_genes, ".svg"),width=8,height=8)
pieplot.med  <- pie(table(factor(siUpregGenes$`integrationCategory-med-dose`, levels = factor.order.gene.categories)), main = paste0("med dose, N =", n_upreg_genes))
print(pieplot.med)
dev.off()

svg(filename=paste0(piechart.location.prefix, "high_dose_", "n", n_upreg_genes, ".svg"),width=8,height=8)
pieplot.high <- pie(table(factor(siUpregGenes$`integrationCategory-high-dose`, levels = factor.order.gene.categories)), main = paste0("high dose, N =", n_upreg_genes))
print(pieplot.high)
dev.off()


####### make stacked histogram plot showing signal integration constant & colored by frequency of each category in the plot

# function: reassign weird categories that may come up when a specific dose integrates in a different direction to "uncategorized" 
mapCatsToReducesCatSet <- function(cat.values) {
  allowedCategories <- c("ambiguous", "sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative")
  cat.values[! cat.values %in% allowedCategories] <- "uncategorized"
  return(factor(cat.values, levels = rev(c("uncategorized", allowedCategories))))
}

# function: assign values to a bin
# input -- a vector of integration contstants
# output -- bin values for the vector. order of values in the vector doesn't change
#   note: values outside the range get assigned to the lowest or highest bin available


# first loop: get the maximum bin y value to standardize the y axis limits when making plots
bin.step.size   <-  0.125
bin.leftmost    <- -3
bin.rightmost   <-  5
plot.width    <-   18
plot.height   <-    8.15
bin.radius    <- bin.step.size / 2
bin.midpoints <- seq(bin.leftmost + bin.step.size, bin.rightmost, by = bin.step.size) - bin.radius


max.bin.vals <- c()
for (dosage in c("low", "med", "high")) {
  hist.values <- pull(filtSiUpregGenes, paste0("integrationConstant-", dosage))
  bin.values  <- assignValuesToHistBin(hist.values, bin.midpoints, bin.radius)
  max.bin.val <- max(table(bin.values)[2:(length(table(bin.values)) - 2)]) # do not select the edges for y limits, we will use "broken bars" to illustrate their N
  max.bin.vals <- c(max.bin.vals, max.bin.val)
}

# second loop: make the stacked bar histograms
for (dosage in c("low", "med", "high")) {
  categorical.values <- pull(filtSiUpregGenes, paste0("integrationCategory-", dosage ,"-dose"))
  hist.values        <- pull(filtSiUpregGenes, paste0("integrationConstant-", dosage))
  mapped.categorical.values <- mapCatsToReducesCatSet(categorical.values)

  stackedBarRes <- makeHistogramOfValues(hist.values, mapped.categorical.values, bin.leftmost, bin.rightmost,
                                             bin.step.size, paste0(dosage, " dose, c-values"), 
                                             xlabel = "c-value", ylabel = "count", color.by.category = T)
  
  stackedBarHist <- stackedBarRes[[1]]
  
  if (use.common.scale) {
    stackedBarHist <- stackedBarHist + ylim(0, max(max.bin.vals))
  }

  stackedBarTib  <- stackedBarRes[[2]]
  stackedBarTib[["dose"]] <- dosage
  threshold.for.reporting.upper.end.of.histogram <- 2
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

  ggsave(paste0(stackedBarHistogram.location.prefix, "upregGeneIntegrationConstants_", dosage, "_dose.svg"), plot = stackedBarHist, width = plot.width, height = plot.height)
}
  
