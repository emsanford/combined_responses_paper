library(tidyverse)
library(here)

joinedTibble <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.tsv'))

subset <- joinedTibble %>%
          filter(n_peaks_nearby_gene == 1, `RA-med_isDeGene` == 1, `RA-med-isDiffPeak`) 
          
  #log2fc is gene, avgfoldchange is peak
qplot(log2(subset$`RA-low-avgFoldchange`), subset$`RA-low_log2fc`) + geom_abline(slope = 1)
qplot(log2(subset$`RA-med-avgFoldchange`), subset$`RA-med_log2fc`) + geom_abline(slope = 1)
qplot(log2(subset$`RA-high-avgFoldchange`), subset$`RA-high_log2fc`) + geom_abline(slope = 1)

subset <- joinedTibble %>%
  filter(n_peaks_nearby_gene == 1, `TGFb-med_isDeGene` == 1, `TGFb-med-isDiffPeak`) 
qplot(log2(subset$`TGFb-low-avgFoldchange`), subset$`TGFb-low_log2fc`) + geom_abline(slope = 1) + geom_abline(slope = -1)
qplot(log2(subset$`TGFb-med-avgFoldchange`), subset$`TGFb-med_log2fc`) + geom_abline(slope = 1) + geom_abline(slope = -1)
qplot(log2(subset$`TGFb-high-avgFoldchange`), subset$`TGFb-high_log2fc`) + geom_abline(slope = 1) + geom_abline(slope = -1)


signal.subset <- 'RA-med'
subsetLogInds <- joinedTibble[[paste0(signal.subset, "-isDiffPeak")]]
joinedTibSubset <- filter(joinedTibble, subsetLogInds)
joinedTibSubset <- joinedTibSubset %>% 
                      group_by(ensg) %>%
                      mutate(n_subset_peaks_nearby_gene = n()) %>%
                      ungroup()

joinedTibSubsetSubset <- joinedTibSubset %>% filter(n_subset_peaks_nearby_gene == 1, joinedTibSubset$`RA-med_isDeGene` == 1)
ggplot(joinedTibSubsetSubset, aes(x=`RA-med_log2fc`, y=log2(`RA-med-avgFoldchange`))) +
  geom_text(mapping = aes(label=gene_name))

ggplot(joinedTibSubsetSubset, aes(y=abs(`RA-med_log2fc`), x=abs(log2(`RA-med-avgFoldchange`)))) +
  geom_point() + geom_abline(slope = 1) + geom_abline(slope = -1) + coord_fixed() + theme_classic(base_size=18) + xlab("absolute value of peak log-fold change") + ylab("absolute value of gene log-fold change") + ggtitle("All RA-responsive genes with one diff peak nearby")

ggplot(joinedTibSubsetSubset, aes(y=`RA-med_log2fc`, x=log2(`RA-med-avgFoldchange`))) +
  geom_point() + geom_abline(slope = 1) + geom_abline(slope = -1) + coord_fixed() + theme_classic(base_size=18) + xlab("peak log-fold change") + ylab("gene log-fold change") + ggtitle("All RA-responsive genes with one diff peak nearby")

p1 <- ggplot(NULL)
#optional do extra rigorous filter
joinedTibSubsetSubset <- filter(joinedTibSubsetSubset, 
                                `RA-med_isDeGene` == 1, 
                                `RA-high_isDeGene` == 1,
                                `RA-high-isDiffPeak`, 
                                `RA-low-isDiffPeak`)
for (g in sample(joinedTibSubsetSubset$ensg, 50)) {
  thisgeneonly <- joinedTibSubsetSubset %>% filter(ensg == g)
  x1 <- log2(thisgeneonly$`RA-low-avgFoldchange`[1])
  x2 <- log2(thisgeneonly$`RA-med-avgFoldchange`[1])
  x3 <- log2(thisgeneonly$`RA-high-avgFoldchange`[1])
  y1 <- thisgeneonly$`RA-low_log2fc`[1]
  y2 <- thisgeneonly$`RA-med_log2fc`[1]
  y3 <- thisgeneonly$`RA-high_log2fc`[1]
  xs <- c(x1, x2, x3)
  ys <- c(y1, y2, y3)
  tibby <- tibble(peaklfc = xs, genelfc = ys)
  tib1 <- tibble(peaklfc = x1, genelfc = y1)
  p1 <- p1 + geom_path(data = tibby, mapping = aes(x = peaklfc, y = genelfc), color = 'light blue')
  p1 <- p1 + geom_point(data = tib1, x = x1, y = y1, color = 'light gray')
  p1 <- p1 + geom_point(data = tib1, x = x2, y = y2, color = 'dark gray')
  p1 <- p1 + geom_point(data = tib1, x = x3, y = y3, color = 'black')
}
p1 + coord_fixed() + theme_classic() + geom_abline(slope=1) + geom_abline(slope= -1)

# same plot copy-pasted but this time absolute value
p1 <- ggplot(NULL)
for (g in sample(joinedTibSubsetSubset$ensg, 50)) {
# for (g in joinedTibSubsetSubset$ensg) {
  thisgeneonly <- joinedTibSubsetSubset %>% filter(ensg == g)
  x1 <- abs(log2(thisgeneonly$`RA-low-avgFoldchange`[1]))
  x2 <- abs(log2(thisgeneonly$`RA-med-avgFoldchange`[1]))
  x3 <- abs(log2(thisgeneonly$`RA-high-avgFoldchange`[1]))
  y1 <- abs(thisgeneonly$`RA-low_log2fc`[1])
  y2 <- abs(thisgeneonly$`RA-med_log2fc`[1])
  y3 <- abs(thisgeneonly$`RA-high_log2fc`[1])
  xs <- c(x1, x2, x3)
  ys <- c(y1, y2, y3)
  tibby <- tibble(peaklfc = xs, genelfc = ys)
  tib1 <- tibble(peaklfc = x1, genelfc = y1)
  p1 <- p1 + geom_path(data = tibby, mapping = aes(x = peaklfc, y = genelfc), color = 'light blue')
  p1 <- p1 + geom_point(data = tib1, x = x1, y = y1, color = 'light gray')
  p1 <- p1 + geom_point(data = tib1, x = x2, y = y2, color = 'dark gray')
  p1 <- p1 + geom_point(data = tib1, x = x3, y = y3, color = 'black')
}
p1 + coord_fixed() + theme_classic(base_size = 18) + geom_abline(slope=1) + xlab("abs value peak log-fold change") + ylab("abs value gene expression log-fold change")

# same plot copy-pasted but this time absolute value and only plot pairings of doses, not all three on same plot
p1 <- ggplot(NULL)
p2 <- ggplot(NULL)
p3 <- ggplot(NULL)
show_labels <- F
show_first_point <- F
show_second_point <- T
# for (g in sample(joinedTibSubsetSubset$ensg, 50)) {
# for (g in filter(joinedTibSubsetSubset, `EtOH-nlDensity_avgTPM` > 5)$ensg) {
# for (g in filter(joinedTibSubsetSubset, `EtOH-nlDensity-avgNormFragmentCounts` > 10)$ensg) {
for (g in filter(joinedTibSubsetSubset, `EtOH-nlDensity-avgNormFragmentCounts` > 10, `EtOH-nlDensity_avgTPM` > 5)$ensg) {
  # for (g in joinedTibSubsetSubset$ensg) {
  thisgeneonly <- joinedTibSubsetSubset %>% filter(ensg == g)
  denomlabel  <- sprintf("%.1f,%.1f", thisgeneonly$`EtOH-nlDensity_avgTPM`[1], thisgeneonly$`EtOH-nlDensity-avgNormFragmentCounts`[1])
  x1 <- abs(log2(thisgeneonly$`RA-low-avgFoldchange`[1]))
  x2 <- abs(log2(thisgeneonly$`RA-med-avgFoldchange`[1]))
  x3 <- abs(log2(thisgeneonly$`RA-high-avgFoldchange`[1]))
  y1 <- abs(thisgeneonly$`RA-low_log2fc`[1])
  y2 <- abs(thisgeneonly$`RA-med_log2fc`[1])
  y3 <- abs(thisgeneonly$`RA-high_log2fc`[1])
  xs1 <- c(x1, x2)
  ys1 <- c(y1, y2)
  xs2 <- c(x2, x3)
  ys2 <- c(y2, y3)
  xs3 <- c(x1, x3)
  ys3 <- c(y1, y3)
  tibby1 <- tibble(peaklfc = xs1, genelfc = ys1)
  tibby2 <- tibble(peaklfc = xs2, genelfc = ys2)
  tibby3 <- tibble(peaklfc = xs3, genelfc = ys3)
  tib1 <- tibble(peaklfc = x1, genelfc = y1)
  p1 <- p1 + geom_path(data = tibby1, mapping = aes(x = peaklfc, y = genelfc), color = 'light blue') #+ xlim(c(0.25, 3.5)) + ylim(c(0.25, 4))
  p2 <- p2 + geom_path(data = tibby2, mapping = aes(x = peaklfc, y = genelfc), color = 'light blue') #+ xlim(c(0.25, 3.5)) + ylim(c(0.25, 4))
  p3 <- p3 + geom_path(data = tibby3, mapping = aes(x = peaklfc, y = genelfc), color = 'light blue') #+ xlim(c(0.25, 3.5)) + ylim(c(0.25, 4))
  
  if (show_labels) {
    p1 <- p1 + geom_text(data = tib1, x=x2, y=y2, label=denomlabel, size=3)
    p2 <- p2 + geom_text(data = tib1, x = x3, y = y3, label = denomlabel, size=3)
    p3 <- p3 + geom_text(data = tib1, x = x3, y = y3, label = denomlabel, size=3)
  }
  if (show_first_point) {
    p1 <- p1 + geom_point(data = tib1, x = x1, y = y1, color = 'light gray')
    p2 <- p2 + geom_point(data = tib1, x = x2, y = y2, color = 'dark gray')
    p3 <- p3 + geom_point(data = tib1, x = x1, y = y1, color = 'dark gray')
    
  }
  if (show_second_point) {
    p1 <- p1 + geom_point(data = tib1, x = x2, y = y2, color = 'black')
    p2 <- p2 + geom_point(data = tib1, x = x3, y = y3, color = 'black')
    p3 <- p3 + geom_point(data = tib1, x = x3, y = y3, color = 'black')
  }
}
p1 <- p1 + coord_fixed() + theme_classic(base_size = 18) + geom_abline(slope=1) + geom_abline(slope= -1) +
  xlab("abs value peak log-fold change") + ylab("abs value gene expression log-fold change") + ggtitle("low -> med dose RA (pt labs: denom TPM, avg atac frag cts)")
p2 <- p2 + coord_fixed() + theme_classic(base_size = 18) + geom_abline(slope=1) + geom_abline(slope= -1) +
  xlab("abs value peak log-fold change") + ylab("abs value gene expression log-fold change") + ggtitle("med -> high dose RA (pt labs: denom TPM, avg atac frag cts)")
p3 <- p3 + coord_fixed() + theme_classic(base_size = 18) + geom_abline(slope=1) + geom_abline(slope= -1) +
  xlab("abs value peak log-fold change") + ylab("abs value gene expression log-fold change") + ggtitle("low -> high dose RA, no small denominators")
p1 <- p1 + xlim(0.25,3) + ylim(0.25,2.5)
p2 <- p2 + xlim(0.25,3) + ylim(0.25,2.5)
p3 <- p3 + xlim(0.25,3) + ylim(0.25,2.5)
library(gridExtra)
grid.arrange(p1,p2,p3, ncol=3)
