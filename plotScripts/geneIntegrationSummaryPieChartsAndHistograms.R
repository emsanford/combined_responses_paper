library(tidyverse)
library(here)

geneTibAnno   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))

# filter out diff peaks for which to do signal integration analysis on
siUpregGenes1 <- geneTibAnno %>% filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`,
                                        `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`,
                                        `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
                                        `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
                                        `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0)

# this definition is "upregulated by both individually OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals
siUpregGenes2 <- geneTibAnno %>% filter(`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene` | ((`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
                                       `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
                                       `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
                                       `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0)

# this definition is "upregulated by both individually OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals AND not upregulated in a single individual condition
siUpregGenes3 <- geneTibAnno %>% filter(`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene` | ((`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
                                       !((`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`) &  (`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & !(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
                                       !((`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`) & !(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) &  (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
                                       `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
                                       `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
                                       `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0,
                                       `TGFb-and-RA-low_log2fc` > 0, `TGFb-and-RA-med_log2fc` > 0, `TGFb-and-RA-high_log2fc` > 0)
table(siUpregGenes1$`integrationCategory-med-dose`)
table(siUpregGenes2$`integrationCategory-med-dose`)
table(siUpregGenes3$`integrationCategory-med-dose`)

genes_of_interest <- setdiff(siUpregGenes2$gene_name, siUpregGenes3$gene_name)
tib.genes.of.interest <- siUpregGenes2 %>% filter(gene_name %in% genes_of_interest)
new.supermultiplicative.geneIDs <- tib.genes.of.interest %>% filter(`integrationCategory-med-dose` == "super-multiplicative") %>% pull("gene_name")
table(tib.genes.of.interest$`integrationCategory-med-dose`)

siUpregGenes <- siUpregGenes3



# svg(filename=here("plots", "rTpmBeeSwarmsAll", sprintf("piechart_%s_%s.svg", signal.effect.direction, "0-low")),width=8,height=8)
pie(table(siUpregGenes$`integrationCategory-low-dose`))
pie(table(siUpregGenes$`integrationCategory-med-dose`))
pie(table(siUpregGenes$`integrationCategory-high-dose`))
siUpregGenes$`integrationConstant-low` %>% qplot() + xlim(-2, 3)
siUpregGenes$`integrationConstant-med` %>% qplot() + xlim(-2, 3)
siUpregGenes$`integrationConstant-high` %>% qplot() + xlim(-2, 3) 
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
  categorical.values <- pull(siUpregGenes, paste0("integrationCategory-", dosage ,"-dose"))
  hist.values        <- pull(siUpregGenes, paste0("integrationConstant-", dosage))
  
  bin.midpoints <- seq(-10, 10, by = 0.25)
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
  ggsave(here("plots", paste0("stackedBarHist_upregGeneIntegrationConstants_", dosage, "_dose.svg")))
}
  



addPredFcDiffMin <- 1
minTpmDiff <- 2

tfp <- siUpregGenes2 %>% 
  filter(siUpregGenes2$`addMultPredFcDiff-low` >= addPredFcDiffMin,
         (siUpregGenes2$`multPred-low` - siUpregGenes2$`addPred-low`) >= minTpmDiff,
         siUpregGenes2$`addMultPredFcDiff-med` >= addPredFcDiffMin,
         (siUpregGenes2$`multPred-med` - siUpregGenes2$`addPred-med`) >= minTpmDiff,
         siUpregGenes2$`addMultPredFcDiff-high` >= addPredFcDiffMin,
         (siUpregGenes2$`multPred-high` - siUpregGenes2$`addPred-high`) >= minTpmDiff) %>%
  
  dplyr::select(gene_name, `integrationConstant-low`, `integrationConstant-med`, `integrationConstant-high`) %>% 
  mutate(order = rank(`integrationConstant-high`))
nrow(tfp)
genes_of_interest <- tfp$gene_name
# filter some kind of cutoff based on difference...

tfp2 <- gather(tfp, `integrationConstant-low`, `integrationConstant-med`, `integrationConstant-high`, value = "c_value", key = "origColName")
tfp2[["order"]][which(tfp2$origColName == "integrationConstant-low")] <- tfp2[["order"]][which(tfp2$origColName == "integrationConstant-low")] - 0.3
tfp2[["order"]][which(tfp2$origColName == "integrationConstant-high")] <- tfp2[["order"]][which(tfp2$origColName == "integrationConstant-high")] + 0.3

p <- ggplot(data = NULL)
  
p <- p +
  geom_point(data = filter(tfp2, origColName == "integrationConstant-high"), aes(order, c_value)) + 
  geom_hline(yintercept = 0) + geom_hline(yintercept = 1)

for (gn in unique(tfp2$gene_name)) {
  gn.only.tib <- tfp2 %>% filter(gene_name == gn)
  p <- p + geom_line(data = gn.only.tib, aes(order, c_value))
}
xtickLocations <- tfp2 %>% filter(origColName == "integrationConstant-med") %>% pull("order")
xtickLabels    <- tfp2 %>% filter(origColName == "integrationConstant-med") %>% pull("gene_name") 
p <- p + scale_x_continuous(breaks = xtickLocations, labels = xtickLabels)
p <- p + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle=90, hjust=1)) + ggtitle(sprintf("integration is ~stable but also tends to saturate with increasing dose? All upreg genes, such that:\nmin FC diff between add and mult pred = %0.2f.\nmin TPM diff between add and mult pred = %0.2f.", addPredFcDiffMin, minTpmDiff))
p 
# p + ylim(-5, 5)
