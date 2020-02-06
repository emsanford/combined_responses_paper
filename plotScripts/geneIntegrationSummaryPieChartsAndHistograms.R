library(tidyverse)
library(here)

siUpregGenes   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv'))
addPredFcDiffMin <- 0.20
minTpmDiff <- 1
output.location.prefix = here('plots', 'gene_summary_plots_')

####### first, make pie charts for the categorical description for each upregulated gene
# svg(filename=here("plots", "rTpmBeeSwarmsAll", sprintf("piechart_%s_%s.svg", signal.effect.direction, "0-low")),width=8,height=8)
pie(table(siUpregGenes$`integrationCategory-low-dose`), main = "low dose")
pie(table(siUpregGenes$`integrationCategory-med-dose`), main = "med dose")
pie(table(siUpregGenes$`integrationCategory-high-dose`), main = "high dose")
# dev.off()


####### make stacked histogram plot showing signal integration constant & colored by frequency of each category in the plot
# function: reassign weird categories that may come up when a different dose integrates in a different direction to "uncategorized" 
filtSiUpregGenes <- siUpregGenes %>% 
  filter(`addMultPredFcDiff-low` >= addPredFcDiffMin,
         (`multPred-low` - `addPred-low`) >= minTpmDiff,
         `addMultPredFcDiff-med` >= addPredFcDiffMin,
         (`multPred-med` - `addPred-med`) >= minTpmDiff,
         `addMultPredFcDiff-high` >= addPredFcDiffMin,
         (`multPred-high` - `addPred-high`) >= minTpmDiff)
# nrow(filtSiUpregGenes)

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
  categorical.values <- pull(filtSiUpregGenes, paste0("integrationCategory-", dosage ,"-dose"))
  hist.values        <- pull(filtSiUpregGenes, paste0("integrationConstant-", dosage))
  
  bin.midpoints <- seq(-3, 5, by = 0.25)
  bin.radius    <- (bin.midpoints[2] - bin.midpoints[1]) / 2
  
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
    geom_vline(xintercept = 0) + geom_vline(xintercept = 1) 
  print(stackedBarHist)
  ggsave(here("plots", paste0("stackedBarHist_upregGeneIntegrationConstants_", dosage, "_dose.svg")), width = 8, height = 4)
}
  
