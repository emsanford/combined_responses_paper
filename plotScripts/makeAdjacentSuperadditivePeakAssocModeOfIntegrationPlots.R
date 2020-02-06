library(tidyverse)
library(here)

superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")
additive.peak.categories      <- c("additive", "ambiguous")
subadditive.peak.categories   <- c("sub-additive")

plot.these.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous")

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes        <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv'))
  siUpregJoinedPeaks  <- read_tsv(here('extractedData', 'joinedTableUpregPeaksNearUpregGenes.tsv'))
  min.ControlTPM      <- 0
  selected.peak.category.arg <- "subadditive"
  outputloc.prefix <- here('plots', paste0('selected_', selected.peak.category.arg,'_peak_stats_near_upreg_genesets_'))
} else {
  siUpregGenes       <- read_tsv(cmdargs[1])
  siUpregJoinedPeaks <- read_tsv(cmdargs[2])
  min.ControlTPM     <- as.numeric(cmdargs[3])
  outputloc.prefix   <- cmdargs[4]
  selected.peak.category.arg <- cmdargs[5] 
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
}

geneSet <- siUpregGenes$ensg

# Plot frequency of selected (e.g., super-additive) peaks near genes of different integration categories.

peaksNearGeneSet       <- siUpregJoinedPeaks %>% filter(ensg %in% geneSet)
genesWithPeaksNearby   <- peaksNearGeneSet$ensg %>% unique()
genesWithNoPeaksNearby <- setdiff(geneSet, genesWithPeaksNearby)
dosages <- c()
intcats <- c()
freqNearbySuperaddPeak <- c()
avgNumNearbySuperaddPeak <- c()
numGenes <- c()
longTibAllGenes <- NULL # can add the ones with zeros here...
for (dosage in c("low", "med", "high")) {
  numGenesWithNoPeaksNearby <- length(genesWithNoPeaksNearby)
  integrationModesThisDose  <- siUpregGenes %>% filter(ensg %in% genesWithNoPeaksNearby) %>% pull(paste0("integrationCategory-", dosage ,"-dose"))
  integrationModesThisDoseFactor <- factor(integrationModesThisDose, levels = c("sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative", "ambiguous"))
  geneIDsThisDose           <- siUpregGenes %>% filter(ensg %in% genesWithNoPeaksNearby) %>% pull("ensg")
  longTibAllGenes <- rbind(longTibAllGenes, tibble(ensg = geneIDsThisDose, 
                                                   dose = dosage, 
                                                   modeOfIntegration = integrationModesThisDoseFactor, 
                                                   numNearbySuperaddPeak = rep(0, numGenesWithNoPeaksNearby)))
}
# fields: gene ID, dosage, number of nearby superadditive peaks, gene's mode of integration
for (dosage in c("low", "med", "high")) {
genesWithNearbyPeakSummaryInfo <- peaksNearGeneSet %>%
  group_by(ensg) %>%
  mutate(is_selected_peak_type = UQ(as.symbol(paste0("peak_integrationCategory-", dosage, "-dose"))) %in% selected.peak.category) %>%
  mutate(numNearbySelectedPeaks = sum(is_selected_peak_type)) %>%
  mutate(hasNearbySelectedPeak = sum(is_selected_peak_type) >= 1) %>%
  dplyr::select(ensg, hasNearbySelectedPeak, numNearbySelectedPeaks, paste0("integrationCategory-", dosage, "-dose")) %>%
  unique()
# to do: add genes that have zero peaks nearby
  for (intmode in c( "sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative",  "ambiguous")) {
    genesInThisDoseAndCategory <- filter(genesWithNearbyPeakSummaryInfo, UQ(as.symbol(paste0("integrationCategory-", dosage, "-dose")))  == intmode)
    freqGenesNearbySuperAddPeakThisCat <- sum(genesInThisDoseAndCategory$hasNearbySelectedPeak)/nrow(genesInThisDoseAndCategory)
    avgNumNearbySuperaddPeakThisCat <- sum(genesInThisDoseAndCategory$numNearbySelectedPeaks)/nrow(genesInThisDoseAndCategory)
    nGenesThisCat <- nrow(genesInThisDoseAndCategory)
    print(sprintf("dosage: %s, int-mode: %s, n-genes-this-cat: %d, fracNearbySuperAdd: %.3f", dosage, intmode, nGenesThisCat, freqGenesNearbySuperAddPeakThisCat))
    dosages <- c(dosages, dosage)
    intcats <- c(intcats, intmode)
    freqNearbySuperaddPeak <- c(freqNearbySuperaddPeak, freqGenesNearbySuperAddPeakThisCat)
    avgNumNearbySuperaddPeak <- c(avgNumNearbySuperaddPeak, avgNumNearbySuperaddPeakThisCat)
    numGenes <- c(numGenes, nGenesThisCat)
    
    longTibAllGenes <- rbind(longTibAllGenes, tibble(ensg = genesInThisDoseAndCategory$ensg, 
                                                     dose = rep(dosage, nGenesThisCat), 
                                                     modeOfIntegration = rep(intmode, nGenesThisCat), 
                                                     numNearbySuperaddPeak = genesInThisDoseAndCategory$numNearbySelectedPeaks))
  }
}
intcats <- as_factor(intcats, levels = c("sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative", "ambiguous"))
dosages <- as_factor(dosages, levels = c("low", "med", "high"))
tfp <- tibble(dosages, intcats, freqNearbySuperaddPeak, avgNumNearbySuperaddPeak, numGenes)
tfp <- filter(tfp, intcats %in% plot.these.gene.categories)
p1 <- ggplot(tfp, aes(x= intcats, y= freqNearbySuperaddPeak)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~dosages) + xlab("integration category") + ylab(paste0("Fraction of genes with >=1 ", plot.category.string, " peak")) + 
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
  ggtitle(paste0(plot.category.string, " peak frequency near genes by integration category"))
ggsave(paste0(outputloc.prefix, "freq_near_genecats.svg"), width = 8, height = 5, plot = p1)

p2 <- ggplot(tfp, aes(x= intcats, y= avgNumNearbySuperaddPeak)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~dosages) + xlab("integration category") + ylab(paste0("Average number of ", plot.category.string, " peaks nearby")) + 
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
  ggtitle(paste0("Number of ", plot.category.string, " peaks near genes by integration category"))
ggsave(paste0(outputloc.prefix, "avgNum_near_genecats.svg"), width = 8, height = 5, plot = p2)

longTibAllGenes[["dose"]] <- factor(longTibAllGenes[["dose"]], levels = c("low", "med", "high"))
filtLongTibAllGenes <- longTibAllGenes %>% filter(modeOfIntegration %in% plot.these.gene.categories)
# filtLongTibAllGenes <- longTibAllGenes
p3 <- ggplot(filtLongTibAllGenes, aes(x= modeOfIntegration, y= numNearbySuperaddPeak)) + 
  geom_boxplot() + 
  facet_wrap(~dose) + xlab("integration category") + ylab(paste0("number of ", plot.category.string, " peaks nearby")) + 
  theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
  ggtitle(paste0("Number of ", plot.category.string, " peaks near genes by integration category"))
print(p3)
ggsave(paste0(outputloc.prefix, "boxplots_numPeaks_near_genecats.svg"), width = 8, height = 5, plot = p3)

