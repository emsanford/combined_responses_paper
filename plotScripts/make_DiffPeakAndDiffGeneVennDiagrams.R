library(tidyverse)
library(here)
library(VennDiagram)
library(venneuler)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  deseqTibAnno     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
  peakTibAnno      <- read_tsv(here('extractedData', 'differentialAtacPeaks_final_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5'))
  output.folder    <- here("plots", "venn_diagrams", '')
} else {
  deseqTibAnno     <- read_tsv(cmdargs[1])
  peakTibAnno      <- read_tsv(cmdargs[2])
  output.folder    <- cmdargs[3]
}

ra.color    <- "#39CF29"
tgfb.color  <- "#2396D6"
both.color  <- "#F29727"
alpha.level <- 1

ra.diff.gene.ids <- filter(deseqTibAnno, `RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) %>% pull("ensg")
ra.upreg.gene.ids <- filter(deseqTibAnno, `RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>% pull("ensg") 

tgfb.diff.gene.ids <- filter(deseqTibAnno, `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`) %>% pull("ensg")
tgfb.upreg.gene.ids <- filter(deseqTibAnno, `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0) %>% pull("ensg") 

both.diff.gene.ids <- filter(deseqTibAnno, `TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`) %>% pull("ensg")
both.upreg.gene.ids <- filter(deseqTibAnno, `TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`, `TGFb-and-RA-high_log2fc` > 0) %>% pull("ensg") 

# install.packages('venneuler')

svg(paste0(output.folder, "allDiffGenes_RA_TGFb_forSetNumbers.svg"))
venn.plot <- venn.diagram(list(ra.diff.gene.ids, tgfb.diff.gene.ids, both.diff.gene.ids), NULL, fill=c(ra.color, tgfb.color, both.color), alpha=c(alpha.level, alpha.level, alpha.level), cex = 2, cat.fontface=4, category.names=c("ra", "tgfb", "both"))
grid.draw(venn.plot)
dev.off()

svg(paste0(output.folder, "allDiffGenes_RA_TGFb.svg"))
df.for.venn.plot1 <- data.frame(elements=c(ra.diff.gene.ids, tgfb.diff.gene.ids, both.diff.gene.ids), 
                        sets=c(rep('RA', length(ra.diff.gene.ids)), rep('TGFb', length(tgfb.diff.gene.ids)), rep('Both', length(both.diff.gene.ids))))
venn.plot1 <- venneuler(df.for.venn.plot1)
plot(venn.plot1)
dev.off()


svg(paste0(output.folder, "upregGenes_RA_TGFb_forSetNumbers.svg"))
venn.plot <- venn.diagram(list(ra.upreg.gene.ids, tgfb.upreg.gene.ids, both.upreg.gene.ids), NULL, fill=c(ra.color, tgfb.color, both.color), alpha=c(alpha.level, alpha.level, alpha.level), cex = 2, cat.fontface=4, category.names=c("ra", "tgfb", "both"))
grid.draw(venn.plot)
dev.off()

svg(paste0(output.folder, "upregGenes_RA_TGFb.svg"))
df.for.venn.plot2 <- data.frame(elements=c(ra.upreg.gene.ids, tgfb.upreg.gene.ids, both.upreg.gene.ids), 
                               sets=c(rep('RA', length(ra.upreg.gene.ids)), rep('TGFb', length(tgfb.upreg.gene.ids)), rep('Both', length(both.upreg.gene.ids))))
venn.plot2 <- venneuler(df.for.venn.plot2)
plot(venn.plot2)
dev.off()


ra.diff.peaks <- peakTibAnno %>%
  filter(peakTibAnno$`RA-low-isDiffPeak` | peakTibAnno$`RA-med-isDiffPeak` | peakTibAnno$`RA-high-isDiffPeak`) %>%
  mutate(loc_string = sprintf("%s:%s", chrom, startLocs)) %>%
  pull("loc_string")

ra.upreg.peaks <- peakTibAnno %>%
  filter(peakTibAnno$`RA-low-isDiffPeak` | peakTibAnno$`RA-med-isDiffPeak` | peakTibAnno$`RA-high-isDiffPeak`, `RA-high-avgFoldchange` > 1) %>%
  mutate(loc_string = sprintf("%s:%s", chrom, startLocs)) %>%
  pull("loc_string")

tgfb.diff.peaks <- peakTibAnno %>%
  filter(peakTibAnno$`TGFb-low-isDiffPeak` | peakTibAnno$`TGFb-med-isDiffPeak` | peakTibAnno$`TGFb-high-isDiffPeak`) %>%
  mutate(loc_string = sprintf("%s:%s", chrom, startLocs)) %>%
  pull("loc_string")

tgfb.upreg.peaks <- peakTibAnno %>%
  filter(peakTibAnno$`TGFb-low-isDiffPeak` | peakTibAnno$`TGFb-med-isDiffPeak` | peakTibAnno$`TGFb-high-isDiffPeak`, `TGFb-high-avgFoldchange` > 1) %>%
  mutate(loc_string = sprintf("%s:%s", chrom, startLocs)) %>%
  pull("loc_string")

both.diff.peaks <- peakTibAnno %>%
  filter(peakTibAnno$`TGFb-and-RA-low-isDiffPeak` | peakTibAnno$`TGFb-and-RA-med-isDiffPeak` | peakTibAnno$`TGFb-and-RA-high-isDiffPeak`) %>%
  mutate(loc_string = sprintf("%s:%s", chrom, startLocs)) %>%
  pull("loc_string")

both.upreg.peaks <- peakTibAnno %>%
  filter(peakTibAnno$`TGFb-and-RA-low-isDiffPeak` | peakTibAnno$`TGFb-and-RA-med-isDiffPeak` | peakTibAnno$`TGFb-and-RA-high-isDiffPeak`, `TGFb-and-RA-high-avgFoldchange` > 1) %>%
  mutate(loc_string = sprintf("%s:%s", chrom, startLocs)) %>%
  pull("loc_string")
  
svg(paste0(output.folder, "allDiffPeaks_RA_TGFb_forSetNumbers.svg"))
venn.plot <- venn.diagram(list(ra.diff.peaks, tgfb.diff.peaks, both.diff.peaks), NULL, fill=c(ra.color, tgfb.color, both.color), alpha=c(alpha.level, alpha.level, alpha.level), cex = 2, cat.fontface=4, category.names=c("ra", "tgfb", "both"))
grid.draw(venn.plot)
dev.off()

svg(paste0(output.folder, "allDiffPeaks_RA_TGFb.svg"))
df.for.venn.plot2 <- data.frame(elements=c(ra.diff.peaks, tgfb.diff.peaks, both.diff.peaks), 
                                sets=c(rep('RA', length(ra.diff.peaks)), rep('TGFb', length(tgfb.diff.peaks)), rep('Both', length(both.diff.peaks))))
venn.plot2 <- venneuler(df.for.venn.plot2)
plot(venn.plot2)
dev.off()


svg(paste0(output.folder, "upregDiffPeaks_RA_TGFb_forSetNumbers.svg"))
venn.plot <- venn.diagram(list(ra.upreg.peaks, tgfb.upreg.peaks, both.upreg.peaks), NULL, fill=c(ra.color, tgfb.color, both.color), alpha=c(alpha.level, alpha.level, alpha.level), cex = 2, cat.fontface=4, category.names=c("ra", "tgfb", "both"))
grid.draw(venn.plot)
dev.off()

svg(paste0(output.folder, "upregDiffPeaks_RA_TGFb.svg"))
df.for.venn.plot2 <- data.frame(elements=c(ra.upreg.peaks, tgfb.upreg.peaks, both.upreg.peaks), 
                                sets=c(rep('RA', length(ra.upreg.peaks)), rep('TGFb', length(tgfb.upreg.peaks)), rep('Both', length(both.upreg.peaks))))
venn.plot2 <- venneuler(df.for.venn.plot2)
plot(venn.plot2)
dev.off()

  
# ### for manually creating them with bioVenn...
# library(data.table)
# fwrite(list(ra.diff.gene.ids), file = here("ra.diff.gene.ids.csv"))
# fwrite(list(ra.upreg.gene.ids), file = here("ra.upreg.gene.ids.csv"))
# fwrite(list(tgfb.diff.gene.ids), file = here("tgfb.diff.gene.ids.csv"))
# fwrite(list(tgfb.upreg.gene.ids), file = here("tgfb.upreg.gene.ids.csv"))
# fwrite(list(both.diff.gene.ids), file = here("both.diff.gene.ids.csv"))
# fwrite(list(both.upreg.gene.ids), file = here("both.upreg.gene.ids.csv"))