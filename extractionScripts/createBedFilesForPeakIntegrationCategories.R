library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  upreg.peaks      <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.upregulatedPeakSet.tsv'))
  outputloc.prefix <- here('extractedData', 'peaks_categorized_by_mode_of_integration', "")
} else {
  upreg.peaks      <- read_tsv(cmdargs[1])
  outputloc.prefix <- (cmdargs[2])
}

subadd.peaks <- upreg.peaks %>% 
  filter(`peak_integrationCategory-med-dose` %in% c("sub-additive")) %>%
  dplyr::select(chrom, startLocs, endLocs)
write_tsv(subadd.peaks, paste0(outputloc.prefix, "_upreg_subadd_peaks.bed"), col_names = F)

add.peaks <- upreg.peaks %>% 
  filter(`peak_integrationCategory-med-dose` %in% c("additive", "ambiguous")) %>%  # note: including ambiguous here at the moment since we're looking at peaks categorically as additive or not additive 
  dplyr::select(chrom, startLocs, endLocs)
write_tsv(add.peaks, paste0(outputloc.prefix, "_upreg_add_peaks.bed"), col_names = F)

superadd.peaks <- upreg.peaks %>% 
  filter(`peak_integrationCategory-med-dose` %in% c("between-add-and-mult", "multiplicative", "super-multiplicative")) %>%
  dplyr::select(chrom, startLocs, endLocs)
write_tsv(superadd.peaks, paste0(outputloc.prefix, "_upreg_superadd_peaks.bed"), col_names = F)



