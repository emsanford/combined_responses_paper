library(tidyverse)
library(GenomicFeatures)
library(here)

##########
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  # deseqtib     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
  # diffpeakstib <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))
  # output.file  <- here('extractedData', 'joinedTablePeaksNearGenes.tsv')
  deseqtib     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv'))
  diffpeakstib <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.upregulatedPeakSet.tsv'))
  output.file  <- here('extractedData', 'joinedTableUpregPeaksNearUpregGenes.tsv')
} else {
  deseqtib     <- read_tsv(cmdargs[1])
  diffpeakstib <- read_tsv(cmdargs[2])
  output.file  <- cmdargs[3]
}



diffpeakstib <- diffpeakstib %>%
                  mutate(peak_chrom = chrom,
                         peak_startLoc = startLocs,
                         peak_endLoc   = endLocs) %>%
                  dplyr::select(-chrom, -startLocs, -endLocs)
TSS.window.radius <- 100000
##########

grTssWindows <- GRanges(seqnames  = c(deseqtib$chrom),
                        ranges    = IRanges(start = deseqtib$TSS_loc - TSS.window.radius,
                                            end   = deseqtib$TSS_loc + TSS.window.radius),
                        strand    = deseqtib$strand,
                        gene_id   = deseqtib$ensg)

grDiffPeaks  <- GRanges(seqnames = diffpeakstib$peak_chrom,
                        ranges = IRanges(start = diffpeakstib$peak_startLoc,
                                         end   = diffpeakstib$peak_endLoc))

TssWindowPeakOverlaps <- findOverlaps(grTssWindows, grDiffPeaks, ignore.strand=TRUE)
gene_indices <- queryHits(TssWindowPeakOverlaps)
peak_indices <- subjectHits(TssWindowPeakOverlaps)

joinedTibble <- as_tibble(cbind(deseqtib[gene_indices,], diffpeakstib[peak_indices,]))

write_tsv(joinedTibble, output.file, col_names = T)