library(tidyverse)
library(here)
library(Matching)
source(here('extractionScripts', 'util.R'))

# modified code from makeAdjacentPeakTypeAssocToGeneTypeModeOfIntegrationPlots.R

superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")
additive.peak.categories      <- c("additive", "ambiguous")
subadditive.peak.categories   <- c("sub-additive")
other.possible.peak.categories <- c("super-additive", "sub-multiplicative", "uncategorized") # these categories only show up when the individual effects downregulate the peak or have opposing effects

factor.order.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous", "between-add-and-mult")
plot.these.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous")

mutual.exclusivity.threshold <- 0.90

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  # siUpregGenes        <- read_tsv(here('extractedData', 'differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated'))
  siUpregGenes             <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
  # siUpregJoinedUpregPeaks  <- read_tsv(here('extractedData', 'upregJoinedPeakGeneTib_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.tsv'))
  siUpregJoinedUpregPeaks  <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/upregGenesUpregPeaksJoinedTib_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.tsv")
  siUpregJoinedAllPeaks    <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/upregGenesAllPeaksJoinedTib_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.tsv") 
  outputloc.prefix <- here('plots', 'peaks_near_gene_types_plots', "")
} else {
  siUpregGenes            <- read_tsv(cmdargs[1])
  siUpregJoinedUpregPeaks <- read_tsv(cmdargs[2])
  siUpregJoinedAllPeaks   <- read_tsv(cmdargs[3])
  outputloc.prefix        <- cmdargs[4]
}

grand.list.of.list.of.plots <- list()
grand.counter <- 1
for (selected.peak.category.arg in c("subadditive", "additive", "superadditive", "mutualExclusivePairs")) {
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
  } else if (selected.peak.category.arg == "mutualExclusivePairs") {
    plot.category.string   <- "mutualExclusivePairs"
  }
  
  geneSet <- siUpregGenes$ensg
  upregPeaksNearGeneSet  <- siUpregJoinedUpregPeaks %>% filter(ensg %in% geneSet)
  allPeaksNearGeneSet    <- siUpregJoinedAllPeaks %>% filter(ensg %in% geneSet)
  
  genesWithUpregPeaksNearby   <- upregPeaksNearGeneSet$ensg %>% unique()
  genesWithNoUpregPeaksNearby <- setdiff(geneSet, genesWithUpregPeaksNearby)
  
  genesWithPeaksNearby   <- allPeaksNearGeneSet$ensg %>% unique()
  genesWithNoPeaksNearby <- setdiff(geneSet, genesWithPeaksNearby)
  
  longTibAllGenes <- NULL # can add the ones with zeros here...
  for (dosage in c("low", "med", "high")) {
    if (selected.peak.category.arg == "mutualExclusivePairs") {
      genesWithPeaksNearby   <- allPeaksNearGeneSet$ensg %>% unique()
      genesWithNoPeaksNearby <- setdiff(geneSet, genesWithPeaksNearby)
    } else {
      genesWithUpregPeaksNearby <- upregPeaksNearGeneSet$ensg %>% unique()
      genesWithNoPeaksNearby <- setdiff(geneSet, genesWithUpregPeaksNearby)
    }
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
    if (selected.peak.category.arg == "mutualExclusivePairs") {
      longTibGenesWithPeaksNearbyThisDose <- allPeaksNearGeneSet %>%
        group_by(ensg) %>%
        mutate(num_ME_peaks_RA   = sum(PeakMutualExclusivityScoreAsymmetricAdditive > mutual.exclusivity.threshold)) %>%
        mutate(num_ME_peaks_TGFb = sum(PeakMutualExclusivityScoreAsymmetricAdditive < (1 - mutual.exclusivity.threshold))) %>%
        mutate(num_ME_peaks_RA_up     = sum((PeakMutualExclusivityScoreAsymmetricAdditive > mutual.exclusivity.threshold)        & (`RA-med-avgFoldchange`   > 1))) %>%
        mutate(num_ME_peaks_TGFb_up   = sum((PeakMutualExclusivityScoreAsymmetricAdditive < (1 - mutual.exclusivity.threshold))  & (`TGFb-med-avgFoldchange` > 1))) %>%
        mutate(num_ME_peaks_RA_down   = sum((PeakMutualExclusivityScoreAsymmetricAdditive > mutual.exclusivity.threshold)        & (`RA-med-avgFoldchange`   < 1))) %>%
        mutate(num_ME_peaks_TGFb_down = sum((PeakMutualExclusivityScoreAsymmetricAdditive < (1 - mutual.exclusivity.threshold))  & (`TGFb-med-avgFoldchange` < 1))) %>%
        mutate(num_ME_peak_pairs = min(num_ME_peaks_RA_up, num_ME_peaks_TGFb_up)) %>%
        mutate(numNearbyPeaksThisType = num_ME_peak_pairs >= 1) %>%
        mutate(modeOfIntegration = UQ(as.symbol(paste0("integrationCategory-", dosage, "-dose")))) %>%
        mutate(dose = dosage) %>%
        dplyr::select(ensg, dose, modeOfIntegration, numNearbyPeaksThisType) %>%
        ungroup() %>%
        unique() 
    } else {
      longTibGenesWithPeaksNearbyThisDose <- upregPeaksNearGeneSet %>%
        group_by(ensg) %>%
        mutate(is_selected_peak_type = UQ(as.symbol(paste0("peak_integrationCategory-", dosage, "-dose"))) %in% selected.peak.category) %>%
        mutate(numNearbyPeaksThisType = sum(is_selected_peak_type)) %>%
        mutate(modeOfIntegration = UQ(as.symbol(paste0("integrationCategory-", dosage, "-dose")))) %>%
        mutate(dose = dosage) %>%
        dplyr::select(ensg, dose, modeOfIntegration, numNearbyPeaksThisType) %>%
        ungroup() %>%
        unique() 
    }
    longTibAllGenes <- rbind(longTibAllGenes, longTibGenesWithPeaksNearbyThisDose)
  }
  
  longTibAllGenes <- filter(longTibAllGenes, modeOfIntegration %in% plot.these.gene.categories)
  longTibAllGenes[["dose"]] <- factor(longTibAllGenes[["dose"]], levels = c("low", "med", "high"))
  longTibAllGenes[["modeOfIntegration"]] <- factor(longTibAllGenes[["modeOfIntegration"]], levels = factor.order.gene.categories)
  
  #### add bootstrap confidence intervals here, using longTibAllGenes
  longTibAllGenes[["n_genes_this_dose_and_intmode"]] <- NA
  set.seed(0)
  for (dosage in c("low", "med", "high")) {
    for (intmode in c(plot.these.gene.categories)) {
      sample.of.n.peak.types.near.gene <- longTibAllGenes %>%
        filter(dose == dosage, modeOfIntegration == intmode) %>%
        pull(numNearbyPeaksThisType)
      n.genes.this.cat <- length(sample.of.n.peak.types.near.gene)
      sample.mean <- mean(sample.of.n.peak.types.near.gene)
     
      longTibAllGenes[["n_genes_this_dose_and_intmode"]][(longTibAllGenes$dose == dosage) & (longTibAllGenes$modeOfIntegration == intmode)] <- n.genes.this.cat
    }
  }
  
  # tibble used for plotting the bar graphs for each dose
  tfp <- longTibAllGenes %>%
    group_by(dose, modeOfIntegration) %>%
    mutate(freqNearbySuperaddPeak   = sum(numNearbyPeaksThisType >= 1) / n()) %>%
    mutate(avgnumNearbyPeaksThisType = mean(numNearbyPeaksThisType)) %>%
    dplyr::select(dose, modeOfIntegration, freqNearbySuperaddPeak, avgnumNearbyPeaksThisType, n_genes_this_dose_and_intmode) %>%
    ungroup() %>%
    unique()
  
  for (dosage in c("low", "med", "high")) {
    comparison1 <- "additive"
    comparison2 <- "multiplicative"
    addNPeaksNearbyMedDose <- longTibAllGenes %>% filter(dose == dosage, modeOfIntegration == comparison1) %>% pull(numNearbyPeaksThisType)
    multNPeaksNearbyMedDose <- longTibAllGenes %>% filter(dose == dosage, modeOfIntegration == comparison2) %>% pull(numNearbyPeaksThisType)
    
    print(sprintf("peak type = %s, dose = %s, comparison = %s vs. %s genes", selected.peak.category.arg, dosage, comparison1, comparison2))
    print(t.test(addNPeaksNearbyMedDose, multNPeaksNearbyMedDose))
  }
}






