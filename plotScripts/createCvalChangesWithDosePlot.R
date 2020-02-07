library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv'))
  addPredFcDiffMin <- 1
  minTpmDiff       <- 2
  output_location <- here('plots', 'cvalsEachDoseForSetOfGenes.svg')
} else {
  siUpregGenes     <- read_tsv(cmdargs[1])
  addPredFcDiffMin <- as.numeric(cmdargs[2])
  minTpmDiff       <- as.numeric(cmdargs[3])
  output_location <- cmdargs[4]
}


tfp <- siUpregGenes %>% 
  filter(`addMultPredFcDiff-low` >= addPredFcDiffMin,
         (`multPred-low` - `addPred-low`) >= minTpmDiff,
         `addMultPredFcDiff-med` >= addPredFcDiffMin,
         (`multPred-med` - `addPred-med`) >= minTpmDiff,
         `addMultPredFcDiff-high` >= addPredFcDiffMin,
         (`multPred-high` - `addPred-high`) >= minTpmDiff) %>%
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
p <- p + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle=90, hjust=1)) + ggtitle(sprintf("integration is ~stable but also tends to saturate with increasing dose? All upreg genes, such that:\nmin FC diff between add and mult pred = %0.2f.\nmin TPM diff between add and mult pred = %0.2f.", addPredFcDiffMin, minTpmDiff))
ggsave(output_location, plot = p, width = 24, height = 6)
# p + ylim(-5, 5)
