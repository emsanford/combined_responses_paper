library(tidyverse)
library(here)

joinedTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv'))
geneTib   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTib   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))


# add desired peak fields to a gene tibble
joinedTib.relevantfields <- joinedTib %>% 
  dplyr::select(ensg, `geneHasNearbySuperadditivePeak-low`,
                      `geneHasNearbySuperadditivePeak-med`,
                      `geneHasNearbySuperadditivePeak-high`) %>% 
  unique()

geneTibAddPeakAggregateAnnotations <- left_join(geneTib, joinedTib.relevantfields, by = 'ensg')
for (dosage in c("low", "med", "high")) {
  geneTibAddPeakAggregateAnnotations[[paste0("geneHasNearbySuperadditivePeak-", dosage)]] <- replace_na(geneTibAddPeakAggregateAnnotations[[paste0("geneHasNearbySuperadditivePeak-", dosage)]], FALSE)
}

# 1. Analysis on upregulated genes
genesUp <- geneTibAddPeakAggregateAnnotations %>% 
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0)


subadd_genes    <- filter(genesUp, `integrationCategory-low-dose` == 'sub-additive')
subadd_frac     <- sum(subadd_genes$`geneHasNearbySuperadditivePeak-low`) / nrow(subadd_genes)
add_genes       <- filter(genesUp, `integrationCategory-low-dose` == 'additive')
add_frac        <- sum(add_genes$`geneHasNearbySuperadditivePeak-low`) / nrow(add_genes)
mult_genes      <- filter(genesUp, `integrationCategory-low-dose` == 'multiplicative')
mult_frac       <- sum(mult_genes$`geneHasNearbySuperadditivePeak-low`) / nrow(mult_genes)
supermult_genes <- filter(genesUp, `integrationCategory-low-dose` == 'super-multiplicative')
supermult_frac  <- sum(supermult_genes$`geneHasNearbySuperadditivePeak-low`) / nrow(supermult_genes)

set.seed(0)
results_vec <- c()
samplesize <- nrow(supermult_genes)
for (ii in 1:25000) {
  this.sample <- sample_n(genesUp, samplesize)
  results_vec <- c(results_vec, sum(this.sample$`geneHasNearbySuperadditivePeak-low`) / samplesize)
}

qplot(results_vec, bins = 100) + 
  geom_vline(xintercept = subadd_frac) +
  geom_vline(xintercept = add_frac) +
  geom_vline(xintercept = mult_frac) +
  geom_vline(xintercept = supermult_frac, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of genes,\nfraction that are near superadditive peaks") +
  xlab("Fraction of genes with a superadditive peak nearby")

print(sprintf("sample size: %d", samplesize))
print(c(subadd_frac, add_frac, mult_frac, supermult_frac))
quantile(results_vec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))





