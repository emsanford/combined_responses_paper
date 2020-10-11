library(tidyverse)
library(here)
library(patchwork)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes     <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
  addPredFcDiffMin <- 1
  minTpmDiff       <- 2
  output_prefix <- here('plots', 'cvals')
} else {
  siUpregGenes     <- read_tsv(cmdargs[1])
  addPredFcDiffMin <- as.numeric(cmdargs[2])
  minTpmDiff       <- as.numeric(cmdargs[3])
  output_prefix <- cmdargs[4]
}

control.tpm.zero.rows <- which(siUpregGenes$`EtOH-nlDensity_avgTPM` == 0)
siUpregGenes <- siUpregGenes[-control.tpm.zero.rows, ] # remove rows with control TPM = 0 because the erroneously estimate c to be zero

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
  geom_point(data = filter(tfp2, origColName == "integrationConstant-high"), aes(order, c_value), color = "#05106E") + 
  geom_point(data = filter(tfp2, origColName == "integrationConstant-med"), aes(order, c_value), color = "#007DFF") + 
  geom_point(data = filter(tfp2, origColName == "integrationConstant-low"), aes(order, c_value), color = "#0DE0FF") + 
  geom_hline(yintercept = 0) + geom_hline(yintercept = 1)

for (gn in unique(tfp2$gene_name)) {
  gn.only.tib <- tfp2 %>% filter(gene_name == gn)
  p <- p + geom_line(data = gn.only.tib, aes(order, c_value))
}
xtickLocations <- tfp2 %>% filter(origColName == "integrationConstant-med") %>% pull("order")
xtickLabels    <- tfp2 %>% filter(origColName == "integrationConstant-med") %>% pull("gene_name") 
p <- p + scale_x_continuous(breaks = xtickLocations, labels = xtickLabels)
p <- p + theme_classic(base_size = 20) + theme(axis.text.x = element_text(angle=90, hjust=1)) + ggtitle(sprintf("integration is ~stable but also tends to saturate with increasing dose? All upreg genes, such that:\nmin FC diff between add and mult pred = %0.2f.\nmin TPM diff between add and mult pred = %0.2f.", addPredFcDiffMin, minTpmDiff))
ggsave(paste0(output_prefix, "EachDoseForSetOfGenes.svg"), plot = p, width = 24, height = 6)
# p + ylim(-5, 5)

### make histograms of c value changes

# all c values
ld.cvals <- siUpregGenes$`integrationConstant-low`
md.cvals <- siUpregGenes$`integrationConstant-med`
hd.cvals <- siUpregGenes$`integrationConstant-high`

low.to.mid.cval.diff  <- md.cvals - ld.cvals
mid.to.high.cval.diff <- hd.cvals - md.cvals
low.to.high.cval.diff <- hd.cvals - ld.cvals

cval.hist.lowerbound <- -4
cval.hist.upperbound <-  4
bin.step.size <-  0.125
max.yval <- 80

clipped.low.to.mid.cval.diff <- ifelse(low.to.mid.cval.diff > cval.hist.upperbound, cval.hist.upperbound, low.to.mid.cval.diff)
clipped.low.to.mid.cval.diff <- ifelse(clipped.low.to.mid.cval.diff < cval.hist.lowerbound, cval.hist.lowerbound, clipped.low.to.mid.cval.diff)
n.gt.upperbound <- sum(low.to.mid.cval.diff > cval.hist.upperbound)
n.lt.lowerbound <- sum(low.to.mid.cval.diff < cval.hist.lowerbound)
phist1 <- qplot(clipped.low.to.mid.cval.diff, binwidth = bin.step.size) + theme_classic() + ggtitle(sprintf('low to med dose\n%d < lower bound\n%d > upper bound', n.lt.lowerbound, n.gt.upperbound)) + ylab("number of genes") + xlab("change in c value") + ylim(0, max.yval) + xlim(cval.hist.lowerbound, cval.hist.upperbound)

clipped.mid.to.high.cval.diff <- ifelse(mid.to.high.cval.diff > cval.hist.upperbound, cval.hist.upperbound, mid.to.high.cval.diff)
clipped.mid.to.high.cval.diff <- ifelse(clipped.mid.to.high.cval.diff < cval.hist.lowerbound, cval.hist.lowerbound, clipped.mid.to.high.cval.diff)
n.gt.upperbound <- sum(mid.to.high.cval.diff > cval.hist.upperbound)
n.lt.lowerbound <- sum(mid.to.high.cval.diff < cval.hist.lowerbound)
phist2 <- qplot(clipped.mid.to.high.cval.diff, binwidth = bin.step.size) + theme_classic() + ggtitle(sprintf('med to high dose\n%d < lower bound\n%d > upper bound', n.lt.lowerbound, n.gt.upperbound)) + ylab("number of genes") + xlab("change in c value") + ylim(0, max.yval) + xlim(cval.hist.lowerbound, cval.hist.upperbound)

clipped.low.to.high.cval.diff <- ifelse(low.to.high.cval.diff > cval.hist.upperbound, cval.hist.upperbound, low.to.high.cval.diff)
clipped.low.to.high.cval.diff <- ifelse(clipped.low.to.high.cval.diff < cval.hist.lowerbound, cval.hist.lowerbound, clipped.low.to.high.cval.diff)
n.gt.upperbound <- sum(low.to.high.cval.diff > cval.hist.upperbound)
n.lt.lowerbound <- sum(low.to.high.cval.diff < cval.hist.lowerbound)
phist3 <- qplot(clipped.low.to.high.cval.diff, binwidth = bin.step.size) + theme_classic() + ggtitle(sprintf('low to high dose\n%d < lower bound\n%d > upper bound', n.lt.lowerbound, n.gt.upperbound)) + ylab("number of genes") + xlab("change in c value") + ylim(0, max.yval) + xlim(cval.hist.lowerbound, cval.hist.upperbound)

patchworkplot <- phist1 + phist2 + phist3
ggsave(paste0(output_prefix, "histogramsCvalChangesIncreasingDose.svg"), plot = patchworkplot, width = 18, height = 6)


### c values that are filtered such that min tpm diff and min fold change diff between ind. conditions is sufficiently high to stabilize the c val estimate
ld.cvals <- tfp$`integrationConstant-low`
md.cvals <- tfp$`integrationConstant-med`
hd.cvals <- tfp$`integrationConstant-high`

low.to.mid.cval.diff  <- md.cvals - ld.cvals
mid.to.high.cval.diff <- hd.cvals - md.cvals
low.to.high.cval.diff <- hd.cvals - ld.cvals

cval.hist.lowerbound <- -4
cval.hist.upperbound <-  4
bin.step.size <-  0.125
max.yval <- 20

clipped.low.to.mid.cval.diff <- ifelse(low.to.mid.cval.diff > cval.hist.upperbound, cval.hist.upperbound, low.to.mid.cval.diff)
clipped.low.to.mid.cval.diff <- ifelse(clipped.low.to.mid.cval.diff < cval.hist.lowerbound, cval.hist.lowerbound, clipped.low.to.mid.cval.diff)
n.gt.upperbound <- sum(low.to.mid.cval.diff > cval.hist.upperbound)
n.lt.lowerbound <- sum(low.to.mid.cval.diff < cval.hist.lowerbound)
phist1 <- qplot(clipped.low.to.mid.cval.diff, binwidth = bin.step.size) + theme_classic() + ggtitle(sprintf('low to med dose\n%d < lower bound\n%d > upper bound', n.lt.lowerbound, n.gt.upperbound)) + ylab("number of genes") + xlab("change in c value") + ylim(0, max.yval) + xlim(cval.hist.lowerbound, cval.hist.upperbound)

clipped.mid.to.high.cval.diff <- ifelse(mid.to.high.cval.diff > cval.hist.upperbound, cval.hist.upperbound, mid.to.high.cval.diff)
clipped.mid.to.high.cval.diff <- ifelse(clipped.mid.to.high.cval.diff < cval.hist.lowerbound, cval.hist.lowerbound, clipped.mid.to.high.cval.diff)
n.gt.upperbound <- sum(mid.to.high.cval.diff > cval.hist.upperbound)
n.lt.lowerbound <- sum(mid.to.high.cval.diff < cval.hist.lowerbound)
phist2 <- qplot(clipped.mid.to.high.cval.diff, binwidth = bin.step.size) + theme_classic() + ggtitle(sprintf('med to high dose\n%d < lower bound\n%d > upper bound', n.lt.lowerbound, n.gt.upperbound)) + ylab("number of genes") + xlab("change in c value") + ylim(0, max.yval) + xlim(cval.hist.lowerbound, cval.hist.upperbound)

clipped.low.to.high.cval.diff <- ifelse(low.to.high.cval.diff > cval.hist.upperbound, cval.hist.upperbound, low.to.high.cval.diff)
clipped.low.to.high.cval.diff <- ifelse(clipped.low.to.high.cval.diff < cval.hist.lowerbound, cval.hist.lowerbound, clipped.low.to.high.cval.diff)
n.gt.upperbound <- sum(low.to.high.cval.diff > cval.hist.upperbound)
n.lt.lowerbound <- sum(low.to.high.cval.diff < cval.hist.lowerbound)
phist3 <- qplot(clipped.low.to.high.cval.diff, binwidth = bin.step.size) + theme_classic() + ggtitle(sprintf('low to high dose\n%d < lower bound\n%d > upper bound', n.lt.lowerbound, n.gt.upperbound)) + ylab("number of genes") + xlab("change in c value") + ylim(0, max.yval) + xlim(cval.hist.lowerbound, cval.hist.upperbound)

patchworkplot <- phist1 + phist2 + phist3
ggsave(paste0(output_prefix, "histogramsCvalChangesIncreasingDose_Filtered.svg"), plot = patchworkplot, width = 18, height = 6)






