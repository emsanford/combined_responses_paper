library(tidyverse)
library(here)
library(grid)
library(gridExtra)

deseq.tib   <- read_tsv(here("extractedData", "DeSeqOutputAllConds.tsv"))
confint.tib <- read_tsv(here("extractedData", "geneExprConfIntsAndAddMultPredictions.tsv"))
dpeaks.tib  <- read_tsv(here("extractedData", "differentialAtacPeaks.tsv"))
##### create TSS tibble and join it with the deseq tibble. 
TSS.loc.tib  <- read_tsv(here('refs', 'EnsgToTssMapping.tsv'))
TSS.loc.tib.extended <- inner_join(mutate(deseq.tib, gene_id = ensg), TSS.loc.tib, by = "gene_id") %>% dplyr::select(-gene_id)
TSS.loc.tib.extended[["tss"]] <- rep(0, nrow(TSS.loc.tib.extended))
TSS.loc.tib.extended[["tss"]][which(TSS.loc.tib.extended$strand == "+")] <- TSS.loc.tib.extended[["tx_start"]][which(TSS.loc.tib.extended$strand == "+")]
TSS.loc.tib.extended[["tss"]][which(TSS.loc.tib.extended$strand == "-")] <- TSS.loc.tib.extended[["tx_end"]][which(TSS.loc.tib.extended$strand == "-")]
##### define set of dose responsive genes in the up and down direction

###### use this block for a dose responsive set
dose.response.factor <- 1.2
gene_set_up <- filter(deseq.tib,
                   `RA-med_isDeGene` == 1,
                   `RA-low_log2fc` > 0,
                   (`RA-low_avgTPM` * dose.response.factor) < (`RA-med_avgTPM`),
                   (`RA-med_avgTPM` * dose.response.factor) < (`RA-high_avgTPM`))$ensg

dose.response.factor <- 1.2
gene_set_down <- filter(deseq.tib,
                      `RA-med_isDeGene` == 1,
                      `RA-low_log2fc` < 0,
                      (`RA-low_avgTPM`) > (`RA-med_avgTPM`  * dose.response.factor),
                      (`RA-med_avgTPM`) > (`RA-high_avgTPM` * dose.response.factor))$ensg

###### use this block for a dose agnostic set
# gene_set_up <- filter(deseq.tib, 
#                       `RA-med_isDeGene` == 1, 
#                       `RA-low_log2fc` > 0)$ensg
# 
# gene_set_down <- filter(deseq.tib, 
#                         `RA-med_isDeGene` == 1, 
#                         `RA-low_log2fc` < 0)$ensg
# 
# TSS.window.radius <- 100000
#####

getNearbyPeakSet <- function(gene_name_vector, TSS.location.tibble, diff.peaks.tibble) {
  this.tsstib <- filter(TSS.location.tibble, ensg %in% gene_name_vector)
  grTssWindowThisGeneSet <-  GRanges(seqnames = c(this.tsstib$chrom),
                                     ranges   = IRanges(start = this.tsstib$tss - TSS.window.radius,
                                                        end   = this.tsstib$tss + TSS.window.radius),
                                     strand   = this.tsstib$strand,
                                     gene     = this.tsstib$ensg)
  
  grDiffPeaks   <-  GRanges(seqnames = diff.peaks.tibble$chrom,
                            ranges = IRanges(start = diff.peaks.tibble$startLocs,
                                             end   = diff.peaks.tibble$endLocs) ,
                            fragCtsControl = diff.peaks.tibble$`EtOH-nlDensity-avgNormFragmentCounts`)
  
  TssWindowPeakOverlaps <- findOverlaps(grTssWindowThisGeneSet, grDiffPeaks)
  nearbyPeakInfo <- diff.peaks.tibble[subjectHits(TssWindowPeakOverlaps), ]
  nearbyPeakInfo[["ensg"]] <- this.tsstib$ensg[queryHits(TssWindowPeakOverlaps)]

  return(nearbyPeakInfo)
}

nearbyUpPeaks   <- getNearbyPeakSet(gene_set_up, TSS.loc.tib.extended, dpeaks.tib) %>%
  filter(`RA-high-isDiffPeak`) %>%
  select(ensg, chrom, startLocs, endLocs, matches("^RA-(low|med|high)-(avg).*")) %>% 
  group_by(ensg) %>%
  mutate(nPeaksNearEnsg = n())

nearbyDownPeaks <- getNearbyPeakSet(gene_set_down, TSS.loc.tib.extended, dpeaks.tib) %>%
  filter(`RA-high-isDiffPeak`) %>%
  select(ensg, chrom, startLocs, endLocs, matches("^RA-(low|med|high)-(avg).*")) %>% 
  group_by(ensg) %>%
  mutate(nPeaksNearEnsg = n())


p1 <- ggplot(data = NULL)
for (gene in sample(unique(nearbyDownPeaks$ensg), 95)) {
  tsubset <- filter(nearbyDownPeaks, ensg == gene)
  # print(tsubset)
  
  for (i in 1:nrow(tsubset)) {
    this.tib <- tibble(x = c(1,2,3), y = c(pull(tsubset[i, ], "RA-low-avgFoldchange")[1], 
                                           pull(tsubset[i, ], "RA-med-avgFoldchange")[1], 
                                           pull(tsubset[i, ], "RA-high-avgFoldchange")[1]))
    tp <- geom_line(data = this.tib, mapping = aes(x = x, y = log2(y)), alpha = 0.12)
    p1 <- p1 + tp
  }
}
p1 <- p1 + ylim(-3, 4) + theme_classic(base_size = 18) + 
  geom_hline(yintercept = 0) + ggtitle("peaks near 95 dose-responsive DOWN genes (RA)") + 
  xlab("low dose <-- medium dose --> high dose") + ylab("log-fold change of nearby peaks") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

p2 <- ggplot(data = NULL)
for (gene in sample(unique(nearbyUpPeaks$ensg), 95)) {
  tsubset <- filter(nearbyUpPeaks, ensg == gene)
  # print(tsubset)
  
  for (i in 1:nrow(tsubset)) {
    this.tib <- tibble(x = c(1,2,3), y = c(pull(tsubset[i, ], "RA-low-avgFoldchange")[1], 
                                           pull(tsubset[i, ], "RA-med-avgFoldchange")[1], 
                                           pull(tsubset[i, ], "RA-high-avgFoldchange")[1]))
    tp <- geom_line(data = this.tib, mapping = aes(x = x, y = log2(y)), alpha = 0.12)
    p2 <- p2 + tp
  }
}
p2 <- p2 + ylim(-3, 4) + theme_classic(base_size = 18) + geom_hline(yintercept = 0) + 
  ggtitle("peaks near 95 dose-responsive UP genes (RA)") + 
  xlab("low dose <-- medium dose --> high dose") + ylab("log-fold change of nearby peaks") +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())

grid.arrange(p1, p2, ncol=2)

nearbyPeakInfo <- diffpeakstib[subjectHits(TssWindowPeakOverlaps), ] %>% 
  mutate(peakmidpoint = (startLocs + endLocs) / 2) %>%
  mutate(relativepeakmidpoint = peakmidpoint - this.tsstib$tss[1])