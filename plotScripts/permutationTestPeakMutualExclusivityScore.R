library(tidyverse)
library(here)

joinedTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv'))
geneTib   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTib   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))

joinedTib2 <- joinedTib %>% 
  # filter(`TGFb-and-RA-med-avgNormFragmentCounts` > 10) #%>%
  filter(`TGFb-and-RA-low-isDiffPeak` == 1 | `TGFb-and-RA-med-isDiffPeak`  == 1 | `TGFb-and-RA-high-isDiffPeak` == 1 ) #%>%
  #filter(`TGFb-and-RA-med-avgFoldchange` > 1) 

percentile.to.test <- .5

# add desired peak fields to a gene tibble
joinedTib.relevantfields <- joinedTib2 %>% 
  dplyr::select(ensg, PeakMutualExclusivityScoreAdditive) %>%
  group_by(ensg) %>%
  mutate(avgPeakMEscore = quantile(PeakMutualExclusivityScoreAdditive, percentile.to.test)) %>%
  ungroup() %>%
  dplyr::select(ensg, avgPeakMEscore) %>%
  unique()

geneTibAddPeakAggregateAnnotations <- inner_join(geneTib, joinedTib.relevantfields, by = 'ensg')


# 1. Analysis on upregulated genes
genesUp <- geneTibAddPeakAggregateAnnotations %>% 
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0)


subadd_genes    <- filter(genesUp, `integrationCategory-med-dose` == 'sub-additive')
subadd_frac     <- sum(subadd_genes$avgPeakMEscore) / nrow(subadd_genes)
add_genes       <- filter(genesUp, `integrationCategory-med-dose` == 'additive')
add_frac        <- sum(add_genes$avgPeakMEscore) / nrow(add_genes)
mult_genes      <- filter(genesUp, `integrationCategory-med-dose` == 'multiplicative')
mult_frac       <- sum(mult_genes$avgPeakMEscore) / nrow(mult_genes)
supermult_genes <- filter(genesUp, `integrationCategory-med-dose` == 'super-multiplicative')
supermult_frac  <- sum(supermult_genes$avgPeakMEscore) / nrow(supermult_genes)

set.seed(0)
results_vec <- c()
samplesize <- nrow(subadd_genes)
for (ii in 1:12000) {
  this.sample <- sample_n(genesUp, samplesize)
  results_vec <- c(results_vec, sum(this.sample$avgPeakMEscore) / samplesize)
}

qplot(results_vec, bins = 100) + 
  geom_vline(xintercept = subadd_frac, color = "pink") +
  geom_vline(xintercept = add_frac, color = "red") +
  geom_vline(xintercept = mult_frac, color = "blue") +
  geom_vline(xintercept = supermult_frac, color = "dark blue") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of genes, mutual exclusivity scores") +
  xlab(sprintf("%.1fth percentile of mutual exclusivity scores", percentile.to.test * 100))

print(sprintf("sample size: %d", samplesize))
print(c(subadd_frac, add_frac, mult_frac, supermult_frac))
quantile(results_vec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))





