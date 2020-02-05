library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  peakTibAnno <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))
  outputLoc   <- here("extractedData", "differentialAtacPeaks.annotated.upregulatedPeakSet.tsv")
} else {
  peakTibAnno  <- read_tsv(cmdargs[1])
  outputLoc <- cmdargs[2]
}




# this definition is "upregulated by combined treatment OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals
siUpregPeaks <- peakTibAnno %>% filter(`TGFb-and-RA-low-isDiffPeak` | `TGFb-and-RA-med-isDiffPeak` | `TGFb-and-RA-high-isDiffPeak` | ((`RA-low-isDiffPeak` | `RA-med-isDiffPeak` | `RA-high-isDiffPeak`) & (`TGFb-low-isDiffPeak` | `TGFb-med-isDiffPeak` | `TGFb-high-isDiffPeak`)),
                                       `RA-low-avgFoldchange` > 1, `TGFb-low-avgFoldchange` > 1,
                                       `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1,
                                       `RA-high-avgFoldchange` > 1, `TGFb-high-avgFoldchange` > 1)

write_tsv(siUpregPeaks, outputLoc, col_names = T)

### filter out diff peaks for which to do signal integration analysis on

# # this definition requires the peaks to be individually upregulated by RA and TGFb
# siUpregPeaks <- peakTibAnno %>% filter(`RA-low-isDiffPeak` | `RA-med-isDiffPeak` | `RA-high-isDiffPeak`,
#                                        `TGFb-low-isDiffPeak` | `TGFb-med-isDiffPeak` | `TGFb-high-isDiffPeak`,
#                                        `RA-low-avgFoldchange` > 1, `TGFb-low-avgFoldchange` > 1,
#                                        `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1,
#                                        `RA-high-avgFoldchange` > 1, `TGFb-high-avgFoldchange` > 1)

# # this definition just requires the peaks to have a measured (but not necessarily "significant") increase in
# # accessibility from each individual signal, as well as be upregulated in the "both" condition
# siUpregPeaks <- peakTibAnno %>% filter(`TGFb-and-RA-low-isDiffPeak` | `TGFb-and-RA-med-isDiffPeak` | `TGFb-and-RA-high-isDiffPeak`,
#                                        `RA-low-avgFoldchange` > 1, `TGFb-low-avgFoldchange` > 1,
#                                        `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1,
#                                        `RA-high-avgFoldchange` > 1, `TGFb-high-avgFoldchange` > 1)  

# # this definition lets us look at peaks at are upregulated by both, but not individually by either
# siUpregPeaks <- peakTibAnno %>% filter(`TGFb-and-RA-low-isDiffPeak` | `TGFb-and-RA-med-isDiffPeak` | `TGFb-and-RA-high-isDiffPeak`,
#                                        ! (`RA-low-isDiffPeak` | `RA-med-isDiffPeak` | `RA-high-isDiffPeak`),
#                                        ! (`TGFb-low-isDiffPeak` | `TGFb-med-isDiffPeak` | `TGFb-high-isDiffPeak`),
#                                        `RA-low-avgFoldchange` > 1, `TGFb-low-avgFoldchange` > 1,
#                                        `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1,
#                                        `RA-high-avgFoldchange` > 1, `TGFb-high-avgFoldchange` > 1)  

# this definition is "upregulated by both individually OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals
# siUpregPeaks <- peakTibAnno %>% filter(`TGFb-and-RA-low-isDiffPeak` | `TGFb-and-RA-med-isDiffPeak` | `TGFb-and-RA-high-isDiffPeak` | ((`RA-low-isDiffPeak` | `RA-med-isDiffPeak` | `RA-high-isDiffPeak`) & (`TGFb-low-isDiffPeak` | `TGFb-med-isDiffPeak` | `TGFb-high-isDiffPeak`)),
#                                        `RA-low-avgFoldchange` > 1, `TGFb-low-avgFoldchange` > 1,
#                                        `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1,
#                                        `RA-high-avgFoldchange` > 1, `TGFb-high-avgFoldchange` > 1)

