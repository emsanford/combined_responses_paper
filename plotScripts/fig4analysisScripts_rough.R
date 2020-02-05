library(tidyverse)
library(here)

joinedTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv'))
geneTib   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTib   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))


genesUp <- geneTib %>%
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0)

genesDown <- geneTib %>%
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` < 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` < 0)


qplot(genesUp$`integrationConstant-low`, bins = 28) + xlim(-4, 9) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + ggtitle("Signal integration constant value histogram\nUpregulated genes (N=293)\nconstant value reference lines: 0=additive, 1=multiplicative\nlow dose") + theme_minimal(base_size = 14)
qplot(genesUp$`integrationConstant-med`, bins = 100) + xlim(-4, 9) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + ggtitle("Signal integration constant value histogram\nUpregulated genes (N=293)\nconstant value reference lines: 0=additive, 1=multiplicative\nmedium dose") + theme_minimal(base_size = 14)
qplot(genesUp$`integrationConstant-high`, bins = 28) + xlim(-4, 9) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + ggtitle("Signal integration constant value histogram\nUpregulated genes (N=293)\nconstant value reference lines: 0=additive, 1=multiplicative\nhigh dose") + theme_minimal(base_size = 14)

jtUpGenes <- joinedTib %>%
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0) 


mutualexclusivityThresholdsOrdered <- c(.25, .33)
numberOfExceptionsGrantedOrdered   <- c(2, 1)
paramSetIndex <- 2

mutualexclusivityThresh   <- mutualexclusivityThresholdsOrdered[paramSetIndex]
meThreshLow               <- mutualexclusivityThresh
meThreshHigh              <- 1 - mutualexclusivityThresh
numberOfExceptionsGranted <- numberOfExceptionsGrantedOrdered[paramSetIndex]


jtUpGenesMEpeaksOnly <- jtUpGenes %>%
  group_by(ensg) %>% 
  mutate(numPeaksNearby = n()) %>%
  mutate(numNonMutualExclusivePeakNearby = sum(!(  (PeakMutualExclusivityScoreAsymmetricAdditive < meThreshLow)  | 
                                                   (PeakMutualExclusivityScoreAsymmetricAdditive > meThreshHigh)   ))) %>%
  filter(numNonMutualExclusivePeakNearby <= numberOfExceptionsGranted) %>%
  mutate(bothTypesOfMutualExclusivePeaksNearby = any(PeakMutualExclusivityScoreAsymmetricAdditive < meThreshLow) & any(PeakMutualExclusivityScoreAsymmetricAdditive > meThreshHigh)) %>%
  filter(bothTypesOfMutualExclusivePeaksNearby)

jtUpGenesMEpeaksOnly %>% ungroup %>% dplyr::select(ensg, peak_chrom, peak_startLoc, peak_endLoc, PeakMutualExclusivityScoreAsymmetricAdditive, numNonMutualExclusivePeakNearby, numPeaksNearby, bothTypesOfMutualExclusivePeaksNearby) %>% print
genesWithMEPeaks <- unique(jtUpGenesMEpeaksOnly$ensg)
geneTibMEPeaks <- filter(geneTib, ensg %in% genesWithMEPeaks)
geneTibMEPeaks %>% dplyr::select(ensg, gene_name, integrationEffectDir, `integrationCategory-low-dose`, `integrationCategory-med-dose`, `integrationCategory-high-dose`)

freqMultiplicativeThisDataset <- sum(geneTibMEPeaks$`integrationCategory-med-dose` == "multiplicative") / nrow(geneTibMEPeaks)
freqMultiplicativeFullDataset <- sum(genesUp$`integrationCategory-med-dose` == "multiplicative") / nrow(genesUp)
print(sprintf("freq multiplicative in nearby-peak-mutual-exclusive-signal-response-enriched genes: %.3f --- freq multiplicative in all upregulated genes: %.3f", freqMultiplicativeThisDataset, freqMultiplicativeFullDataset))
# peaksUp <- peakTib %>%
#   mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
#   mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
#   filter(`TGFb-med-avgFoldchange` > 1) %>%
#   filter(`RA-med-avgFoldchange` > 1) %>%
#   filter(`isUpRegulatedPeak-med`) %>%
#   dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
#   unique()
# # filter(`peakAdditivePredFcResidual-med` > 1.5)
# 
# peaksUp.nonSuperAdditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
# peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)
