library(tidyverse)
library(here)


joinedPeakTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.tsv'))
threshold.superadditivepeak.fc.difference <- 1.5
outputFileLoc <- here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv')

# annotations to add:

# gene-based annotation: superadditive peak nearby for each dose
joinedPeakTibAnno <- joinedPeakTib
for (dosage in c("low", "med", "high")) {
  temp <- joinedPeakTib %>%
    group_by(ensg) %>%
    mutate(temp = UQ(as.symbol(paste0("peakAdditivePredFcResidual-", dosage))) > threshold.superadditivepeak.fc.difference) %>%
    mutate(numNearbySuperAddPeaks = sum(temp)) %>%
    mutate(hasNearbySuperAddPeak = sum(temp) >= 1) 
  joinedPeakTibAnno[[paste0("geneHasNearbySuperadditivePeak-", dosage)]] <- temp$hasNearbySuperAddPeak
  joinedPeakTibAnno[[paste0("numNearbySuperAdditivePeaks-", dosage)]] <- temp$numNearbySuperAddPeaks
}

# gene-based annotation: number of peaks nearby and for each condition-dosage
joinedPeakTibAnno <- joinedPeakTibAnno %>%
  group_by(ensg) %>%
  mutate(n_peaks_nearby_gene = n()) %>%
  mutate("nPeaksNearby_RA-low" = sum(`RA-low-isDiffPeak`)) %>%
  mutate("nPeaksNearby_RA-med" = sum(`RA-med-isDiffPeak`)) %>%
  mutate("nPeaksNearby_RA-high" = sum(`RA-high-isDiffPeak`)) %>%
  mutate("nPeaksNearby_TGFb-low" = sum(`TGFb-low-isDiffPeak`)) %>%
  mutate("nPeaksNearby_TGFb-med" = sum(`TGFb-med-isDiffPeak`)) %>%
  mutate("nPeaksNearby_TGFb-high" = sum(`TGFb-high-isDiffPeak`)) %>%
  mutate("nPeaksNearby_TGFb-and-RA-low" = sum(`TGFb-and-RA-low-isDiffPeak`)) %>%
  mutate("nPeaksNearby_TGFb-and-RA-med" = sum(`TGFb-and-RA-med-isDiffPeak`)) %>%
  mutate("nPeaksNearby_TGFb-and-RA-high" = sum(`TGFb-and-RA-high-isDiffPeak`)) %>%
  ungroup()


write_tsv(joinedPeakTibAnno, outputFileLoc, col_names = T)