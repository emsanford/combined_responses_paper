library(tidyverse)
library(here)
library(patchwork)

source(here('extractionScripts', 'util.R'))


superadditive.peak.categories <- c("between-add-and-mult", "multiplicative", "super-multiplicative")
additive.peak.categories      <- c("additive", "ambiguous")
subadditive.peak.categories   <- c("sub-additive")
other.possible.peak.categories <- c("super-additive", "sub-multiplicative", "uncategorized") # these categories only show up when the individual effects downregulate the peak or have opposing effects

factor.order.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous", "between-add-and-mult")
plot.these.gene.categories <- c("sub-additive", "additive", "multiplicative", "super-multiplicative", "ambiguous")

mutual.exclusivity.threshold <- 0.90
n.bootstrap.samples <- 10000

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
        mutate(num_ME_peak_pairs = min(num_ME_peaks_RA, num_ME_peaks_TGFb)) %>%
        mutate(numNearbyPeaksThisType = num_ME_peak_pairs) %>%
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
  longTibAllGenes[["bootstrap_ci_0.05"]] <- NA
  longTibAllGenes[["bootstrap_ci_0.95"]] <- NA
  longTibAllGenes[["n_genes_this_dose_and_intmode"]] <- NA
  set.seed(0)
  for (dosage in c("low", "med", "high")) {
    for (intmode in c(plot.these.gene.categories)) {
      sample.of.n.peak.types.near.gene <- longTibAllGenes %>%
        filter(dose == dosage, modeOfIntegration == intmode) %>%
        pull(numNearbyPeaksThisType)
      n.genes.this.cat <- length(sample.of.n.peak.types.near.gene)
      sample.mean <- mean(sample.of.n.peak.types.near.gene)
      bootstrap.deltas <- c()
      for (ii in 1:n.bootstrap.samples) {
        bootstrap.sample <- sample(sample.of.n.peak.types.near.gene, n.genes.this.cat, replace = TRUE)
        this.bootstrap.mean  <- mean(bootstrap.sample)
        this.bootstrap.delta <- this.bootstrap.mean - sample.mean
        bootstrap.deltas <- c(bootstrap.deltas, this.bootstrap.delta)
      }
      lower.delta <- quantile(bootstrap.deltas, 0.05)
      upper.delta <- quantile(bootstrap.deltas, 0.95)
      ci.lower    <- sample.mean - upper.delta
      ci.upper    <- sample.mean - lower.delta
      longTibAllGenes[["bootstrap_ci_0.05"]][(longTibAllGenes$dose == dosage) & (longTibAllGenes$modeOfIntegration == intmode)] <- ci.lower
      longTibAllGenes[["bootstrap_ci_0.95"]][(longTibAllGenes$dose == dosage) & (longTibAllGenes$modeOfIntegration == intmode)] <- ci.upper
      longTibAllGenes[["n_genes_this_dose_and_intmode"]][(longTibAllGenes$dose == dosage) & (longTibAllGenes$modeOfIntegration == intmode)] <- n.genes.this.cat
    }
  }
  
  tfp <- longTibAllGenes %>%
    group_by(dose, modeOfIntegration) %>%
    mutate(freqNearbySuperaddPeak   = sum(numNearbyPeaksThisType >= 1) / n()) %>%
    mutate(avgnumNearbyPeaksThisType = mean(numNearbyPeaksThisType)) %>%
    dplyr::select(dose, modeOfIntegration, freqNearbySuperaddPeak, avgnumNearbyPeaksThisType, bootstrap_ci_0.05, bootstrap_ci_0.95, n_genes_this_dose_and_intmode) %>%
    ungroup() %>%
    unique()
  
  p1 <- ggplot(tfp, aes(x= modeOfIntegration, y= freqNearbySuperaddPeak)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(~dose) + xlab("integration category") + ylab(paste0("Fraction of genes with >=1 ", plot.category.string, " peak")) + 
    theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
    ggtitle(paste0(plot.category.string, " peak frequency near genes by integration category"))
  ggsave(paste0(outputloc.prefix, plot.category.string, "_freq_near_genecats.svg"), width = 8, height = 5, plot = p1)
  
  p2 <- ggplot(tfp, aes(x= modeOfIntegration, y = avgnumNearbyPeaksThisType, 
                        ymin = bootstrap_ci_0.05, ymax = bootstrap_ci_0.95,
                        fill = modeOfIntegration)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    geom_errorbar(width = 0) + 
    facet_wrap(~dose) + xlab("integration category") + ylab(paste0("Average number of ", plot.category.string, " peaks nearby")) + 
    theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
    ggtitle(paste0("Number of ", plot.category.string, " peaks near genes by integration category")) + guides(fill=FALSE)
  ggsave(paste0(outputloc.prefix, plot.category.string, "_avgNum_near_genecats.svg"), width = 8, height = 5, plot = p2)
  
  # now let's make individual plots so we can customize x-tick labels and compose them later with patchwork
  single.avg.num.peaks.nearby.plot.list <- list()
  counter <- 1
  for (dosage in c('low', 'med', 'high')) {
    filt.tfp <- tfp %>% 
      filter(dose == dosage) %>%
      arrange( modeOfIntegration)
    p <- filt.tfp %>% 
      ggplot(aes(x= modeOfIntegration, y = avgnumNearbyPeaksThisType, 
                      ymin = bootstrap_ci_0.05, ymax = bootstrap_ci_0.95,
                      fill = modeOfIntegration)) + 
        geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
        geom_errorbar(width = 0) + 
        xlab("") + ylab(paste0("Average number of ", plot.category.string, " peaks nearby"))  +
        scale_x_discrete(labels= paste0(filt.tfp$modeOfIntegration, ", N = ", filt.tfp$n_genes_this_dose_and_intmode)) + 
        theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1)) + 
        ggtitle(paste0(dosage, " dose")) + guides(fill=FALSE)
    single.avg.num.peaks.nearby.plot.list[[counter]] <- p
    counter <- counter + 1
  }
  # standardize y axes
  std.ymax <- max(ggplot_build(single.avg.num.peaks.nearby.plot.list[[1]])$layout$panel_scales_y[[1]]$range$range[2],
                  ggplot_build(single.avg.num.peaks.nearby.plot.list[[2]])$layout$panel_scales_y[[1]]$range$range[2],
                  ggplot_build(single.avg.num.peaks.nearby.plot.list[[3]])$layout$panel_scales_y[[1]]$range$range[2])
  for (ii in 1:3) {
    single.avg.num.peaks.nearby.plot.list[[ii]] <- single.avg.num.peaks.nearby.plot.list[[ii]] + ylim(0, std.ymax)
  }
  grand.list.of.list.of.plots[[grand.counter]] <- single.avg.num.peaks.nearby.plot.list
  grand.counter <- grand.counter + 1
  
  # filtLongTibAllGenes <- longTibAllGenes
  p3 <- ggplot(longTibAllGenes, aes(x= modeOfIntegration, y= numNearbyPeaksThisType)) + 
    geom_boxplot() + 
    facet_wrap(~dose) + xlab("integration category") + ylab(paste0("number of ", plot.category.string, " peaks nearby")) + 
    theme_classic(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4)) + 
    ggtitle(paste0("Number of ", plot.category.string, " peaks near genes by integration category"))
  ggsave(paste0(outputloc.prefix, plot.category.string, "_boxplots_numPeaks_near_genecats.svg"), width = 8, height = 5, plot = p3)
}

# make 3 x 3 grid with patchwork
patchwork.plot <- (grand.list.of.list.of.plots[[1]][[1]] + grand.list.of.list.of.plots[[1]][[2]] + grand.list.of.list.of.plots[[1]][[3]]) /
                  (grand.list.of.list.of.plots[[2]][[1]] + grand.list.of.list.of.plots[[2]][[2]] + grand.list.of.list.of.plots[[2]][[3]]) /
                  (grand.list.of.list.of.plots[[3]][[1]] + grand.list.of.list.of.plots[[3]][[2]] + grand.list.of.list.of.plots[[3]][[3]])
ggsave(paste0(outputloc.prefix, "patchwork_plot_avgNumPeakTypesNearGeneTypes.svg"), plot = patchwork.plot, width = 24, height = 16)

patchwork.plot2 <- (grand.list.of.list.of.plots[[4]][[1]] + grand.list.of.list.of.plots[[4]][[2]] + grand.list.of.list.of.plots[[4]][[3]])
ggsave(paste0(outputloc.prefix, "patchwork_plot_mePairsNearNearGeneTypes.svg"), plot = patchwork.plot2, width = 24, height = 16 / 3)







# supermultiplicative gene analysis, looking at genes where one signal seems to require the other signal to have an effect: 
     # this analysis may need to be thrown out
mutual.exclusivity.threshold <- 0.90
siUpregJoinedAllPeaks
siUpregJoinedUpregPeaks

gene.mutual.exclusivity.threshold <- 0.75
peak.mutual.exclusivity.threshold <- 0.75

# for now, we are ignoring genes w/ no peaks nearby
ra.dominant.supermult.genes <- siUpregJoinedAllPeaks %>%
  filter(`integrationCategory-med-dose` == "super-multiplicative") %>%
  filter(gene_asymmmetricMutualExclusivityScore > gene.mutual.exclusivity.threshold) %>%
  # filter(PeakMutualExclusivityScoreAsymmetricAdditive > peak.mutual.exclusivity.threshold) %>%
  group_by(gene_name) %>%
  mutate(max.d.val.of.exclusive.peak = max(`peakAdditivePredFcResidual-med`)) %>%
  ungroup() %>%
  dplyr::select(gene_name, max.d.val.of.exclusive.peak, gene_asymmmetricMutualExclusivityScore) %>%
  unique() %>%
  pull(max.d.val.of.exclusive.peak) %>% qplot() + xlim(-3, 3)

tgfb.dominant.supermult.genes <- siUpregJoinedAllPeaks %>%
  filter(`integrationCategory-med-dose` == "super-multiplicative") %>%
  filter(gene_asymmmetricMutualExclusivityScore < (1 - gene.mutual.exclusivity.threshold)) %>%
  # filter(PeakMutualExclusivityScoreAsymmetricAdditive < (1 - peak.mutual.exclusivity.threshold)) %>%
  group_by(gene_name) %>%
  mutate(max.d.val.of.exclusive.peak = max(`peakAdditivePredFcResidual-med`)) %>%
  ungroup() %>%
  dplyr::select(gene_name, max.d.val.of.exclusive.peak, gene_asymmmetricMutualExclusivityScore) %>%
  unique() %>%
  pull(max.d.val.of.exclusive.peak) %>%  qplot() + xlim(-3, 3)


ra.dominant.add.or.mult.genes <- siUpregJoinedAllPeaks %>%
  filter(`integrationCategory-med-dose` %in% c('additive', 'ambiguous', 'multiplicative')) %>%
  filter(gene_asymmmetricMutualExclusivityScore > gene.mutual.exclusivity.threshold) %>%
  # filter(PeakMutualExclusivityScoreAsymmetricAdditive > peak.mutual.exclusivity.threshold) %>%
  group_by(gene_name) %>%
  mutate(max.d.val.of.exclusive.peak = max(`peakAdditivePredFcResidual-med`)) %>%
  ungroup() %>%
  dplyr::select(gene_name, max.d.val.of.exclusive.peak, gene_asymmmetricMutualExclusivityScore) %>%
  unique() %>%
  pull(max.d.val.of.exclusive.peak) %>%  qplot() + xlim(-3, 3)
  
tgfb.dominant.add.or.mult.genes <- siUpregJoinedAllPeaks %>%
  filter(`integrationCategory-med-dose` %in% c('additive', 'ambiguous', 'multiplicative')) %>%
  filter(gene_asymmmetricMutualExclusivityScore < (1 - gene.mutual.exclusivity.threshold)) %>%
  # filter(PeakMutualExclusivityScoreAsymmetricAdditive < (1 - peak.mutual.exclusivity.threshold)) %>%
  group_by(gene_name) %>%
  mutate(max.d.val.of.exclusive.peak = max(`peakAdditivePredFcResidual-med`)) %>%
  ungroup() %>%
  dplyr::select(gene_name, max.d.val.of.exclusive.peak, gene_asymmmetricMutualExclusivityScore) %>%
  unique() %>%
  pull(max.d.val.of.exclusive.peak) %>%  qplot() + xlim(-3, 3)





