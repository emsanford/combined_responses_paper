library(tidyverse)
library(here)

superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")
additive.peak.categories      <- c("additive", "ambiguous")
subadditive.peak.categories   <- c("sub-additive")
other.possible.peak.categories <- c("super-additive", "sub-multiplicative", "uncategorized") # these categories only show up when the individual effects downregulate the peak or have opposing effects

factor.order.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous", "between-add-and-mult")
plot.these.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous")

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes        <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv'))
  # siUpregJoinedPeaks  <- read_tsv(here('extractedData', 'joinedTableAllDiffPeaksNearUpregGenes.tsv'))
  siUpregJoinedPeaks  <- read_tsv(here('extractedData', 'joinedTableAllDiffPeaksNearUpregGenes.tsv'))
  # siUpregJoinedPeaks  <- siUpregJoinedPeaks %>% filter(peak_integrationEffectDir == "both up")
  selected.peak.category.arg <- "superadditive"
  outputloc.prefix <- here('plots', paste0('selected_', selected.peak.category.arg,'_peak_stats_near_upreg_genesets_'))
} else {
  siUpregGenes       <- read_tsv(cmdargs[1])
  siUpregJoinedPeaks <- read_tsv(cmdargs[2])
  outputloc.prefix   <- cmdargs[3]
  selected.peak.category.arg <- cmdargs[4] 
}

if (selected.peak.category.arg == "superadditive") {
  selected.peak.category <- superadditive.peak.categories
  plot.category.string   <- "super-additive"
} else if (selected.peak.category.arg == "additive") {
  selected.peak.category <- additive.peak.categories
  plot.category.string   <- "additive"
} else if (selected.peak.category.arg == "subadditive") {
  selected.peak.category <- subadditive.peak.categories
  plot.category.string   <- "sub-additive"
} else if (selected.peak.category.arg == "all") {
  selected.peak.category <- c(subadditive.peak.categories, additive.peak.categories, superadditive.peak.categories, other.possible.peak.categories)
  plot.category.string   <- "any upregulated"
}

geneSet <- siUpregGenes$ensg

# Plot frequency of selected (e.g., super-additive) peaks near genes of different integration categories.

peaksNearGeneSet       <- siUpregJoinedPeaks %>% filter(ensg %in% geneSet)
genesWithPeaksNearby   <- peaksNearGeneSet$ensg %>% unique()
genesWithNoPeaksNearby <- setdiff(geneSet, genesWithPeaksNearby)

longTibAllGenes <- NULL # can add the ones with zeros here...
for (dosage in c("low", "med", "high")) {
  numGenesWithNoPeaksNearby <- length(genesWithNoPeaksNearby)
  integrationModesThisDose  <- siUpregGenes %>% filter(ensg %in% genesWithNoPeaksNearby) %>% pull(paste0("integrationCategory-", dosage ,"-dose"))
  geneIDsThisDose           <- siUpregGenes %>% filter(ensg %in% genesWithNoPeaksNearby) %>% pull("ensg")
  longTibAllGenes <- rbind(longTibAllGenes, tibble(ensg = geneIDsThisDose, 
                                                   dose = dosage, 
                                                   modeOfIntegration = integrationModesThisDose, 
                                                   numNearbyPeaksThisType = rep(as.integer(0), numGenesWithNoPeaksNearby)))
}
# fields: gene ID, dosage, number of nearby superadditive peaks, gene's mode of integration
for (dosage in c("low", "med", "high")) {
  longTibGenesWithPeaksNearbyThisDose <- peaksNearGeneSet %>%
    group_by(ensg) %>%
    mutate(is_selected_peak_type = UQ(as.symbol(paste0("peak_integrationCategory-", dosage, "-dose"))) %in% selected.peak.category) %>%
    mutate(numNearbyPeaksThisType = sum(is_selected_peak_type)) %>%
    mutate(modeOfIntegration = UQ(as.symbol(paste0("integrationCategory-", dosage, "-dose")))) %>%
    mutate(dose = dosage) %>%
    dplyr::select(ensg, dose, modeOfIntegration, numNearbyPeaksThisType) %>%
    ungroup() %>%
    unique() 
    
  longTibAllGenes <- rbind(longTibAllGenes, longTibGenesWithPeaksNearbyThisDose)

}

longTibAllGenes <- filter(longTibAllGenes, modeOfIntegration %in% plot.these.gene.categories)
longTibAllGenes[["dose"]] <- factor(longTibAllGenes[["dose"]], levels = c("low", "med", "high"))
longTibAllGenes[["modeOfIntegration"]] <- factor(longTibAllGenes[["modeOfIntegration"]], levels = factor.order.gene.categories)

tfp <- longTibAllGenes %>%
  group_by(dose, modeOfIntegration) %>%
  mutate(freqNearbySuperaddPeak   = sum(numNearbyPeaksThisType >= 1) / n()) %>%
  mutate(avgnumNearbyPeaksThisType = mean(numNearbyPeaksThisType)) %>%
  mutate(numGenesThisCategory = n()) %>%
  dplyr::select(dose, modeOfIntegration, freqNearbySuperaddPeak, avgnumNearbyPeaksThisType, numGenesThisCategory) %>%
  ungroup() %>%
  unique()

p1 <- ggplot(tfp, aes(x= modeOfIntegration, y= freqNearbySuperaddPeak)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~dose) + xlab("integration category") + ylab(paste0("Fraction of genes with >=1 ", plot.category.string, " peak")) + 
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
  ggtitle(paste0(plot.category.string, " peak frequency near genes by integration category"))
ggsave(paste0(outputloc.prefix, "freq_near_genecats.svg"), width = 8, height = 5, plot = p1)

p2 <- ggplot(tfp, aes(x= modeOfIntegration, y= avgnumNearbyPeaksThisType)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~dose) + xlab("integration category") + ylab(paste0("Average number of ", plot.category.string, " peaks nearby")) + 
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
  ggtitle(paste0("Number of ", plot.category.string, " peaks near genes by integration category"))
ggsave(paste0(outputloc.prefix, "avgNum_near_genecats.svg"), width = 8, height = 5, plot = p2)

# filtLongTibAllGenes <- longTibAllGenes
p3 <- ggplot(longTibAllGenes, aes(x= modeOfIntegration, y= numNearbyPeaksThisType)) + 
  geom_boxplot() + 
  facet_wrap(~dose) + xlab("integration category") + ylab(paste0("number of ", plot.category.string, " peaks nearby")) + 
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
  ggtitle(paste0("Number of ", plot.category.string, " peaks near genes by integration category"))
ggsave(paste0(outputloc.prefix, "boxplots_numPeaks_near_genecats.svg"), width = 8, height = 5, plot = p3)

