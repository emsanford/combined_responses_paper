library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  diffpeak.tib <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.upregulatedPeakSet.tsv'))
  outputloc    <- here('extractedData', 'differentialAtacPeaks.upregulated.bed')
} else {
  diffpeak.tib <- read_tsv(cmdargs[1])
  outputloc   <- cmdargs[2]
}

# optional: add diffpeak summary info / fold-change string?
bed_file_fields <- diffpeak.tib %>% 
  mutate(fold_change_string = sprintf("%0.2f_%0.2f_%0.2f", `RA-med-avgFoldchange`, `TGFb-med-avgFoldchange`, `TGFb-and-RA-med-avgFoldchange`)) %>%
  dplyr::select(chrom, startLocs, endLocs, fold_change_string) 
write_tsv(bed_file_fields, outputloc, col_names = F)
