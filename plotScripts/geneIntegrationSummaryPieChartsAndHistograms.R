library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv'))
  addPredFcDiffMin <- 0.20
  minTpmDiff       <- 1
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

# first loop: get the maximum bin y value to standardize the y axis limits when making plots
bin.step.size   <- 0.10
mid.point.shift <- bin.step.size / 2
bin.midpoints <- seq(-2.5 + bin.step.size, 5, by = bin.step.size) - mid.point.shift
bin.radius    <- (bin.midpoints[2] - bin.midpoints[1]) / 2

max.bin.vals <- c()
for (dosage in c("low", "med", "high")) {
  hist.values <- pull(filtSiUpregGenes, paste0("integrationConstant-", dosage))
  bin.values  <- assignValuesToHistBin(hist.values, bin.midpoints, bin.radius)
  max.bin.val <- max(table(bin.values))
  max.bin.vals <- c(max.bin.vals, max.bin.val)
}

# second loop: make the stacked bar histograms
for (dosage in c("low", "med", "high")) {
  categorical.values <- pull(filtSiUpregGenes, paste0("integrationCategory-", dosage ,"-dose"))
  hist.values        <- pull(filtSiUpregGenes, paste0("integrationConstant-", dosage))
  
  mapped.categorical.values <- mapCatsToReducesCatSet(categorical.values)
  bin.values <- assignValuesToHistBin(hist.values, bin.midpoints, bin.radius)
  
  stackedBarHistTib <- tibble(intConstantHhistBin = bin.values, intCategory = mapped.categorical.values)
  
  stackedBarHist <- ggplot(stackedBarHistTib, aes(x = intConstantHhistBin, fill = intCategory)) +
    geom_bar(stat="count", width = bin.radius * 2 * .90) + 
    theme_minimal(base_size = 12) + 
    # ylim(0, 45) +
    xlab("integration constant value for a gene") +
    ylab("number of genes") +
    ggtitle(paste0("Distribution of integration constants for upregulated genes\n", dosage, " dose")) +
    geom_vline(xintercept = 0) + geom_vline(xintercept = 1) +
    ylim(0, max(max.bin.vals) * 1.05) + 
    xlim(min(bin.midpoints) - bin.radius, max(bin.midpoints) + bin.radius)
  
  ggsave(paste0(stackedBarHistogram.location.prefix, "upregGeneIntegrationConstants_", dosage, "_dose.svg"), width = 8, height = 4)
}
  
