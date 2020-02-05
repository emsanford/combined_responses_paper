library(tidyverse)
library(here)

joinedTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv'))
geneTib   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTib   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))

# add desired peak fields to a gene tibble
joinedTib.relevantfields <- joinedTib %>% 
  dplyr::select(ensg, `nPeaksNearby_TGFb-and-RA-low`,
                      `nPeaksNearby_TGFb-and-RA-med`,
                      `nPeaksNearby_TGFb-and-RA-high`) %>% 
  unique()

geneTibAddPeakAggregateAnnotations <- left_join(geneTib, joinedTib.relevantfields, by = 'ensg')
for (dosage in c("low", "med", "high")) {
  geneTibAddPeakAggregateAnnotations[[paste0("nPeaksNearby_TGFb-and-RA-", dosage)]] <- replace_na(geneTibAddPeakAggregateAnnotations[[paste0("nPeaksNearby_TGFb-and-RA-", dosage)]], 0)
}

# 1. Analysis on upregulated genes
genesUp <- geneTibAddPeakAggregateAnnotations %>% 
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0)


subadd_genes    <- filter(genesUp, `integrationCategory-med-dose` == 'sub-additive')
subadd_frac     <- sum(subadd_genes$`nPeaksNearby_TGFb-and-RA-med`) / nrow(subadd_genes)
add_genes       <- filter(genesUp, `integrationCategory-med-dose` == 'additive')
add_frac        <- sum(add_genes$`nPeaksNearby_TGFb-and-RA-med`) / nrow(add_genes)
mult_genes      <- filter(genesUp, `integrationCategory-med-dose` == 'multiplicative')
mult_frac       <- sum(mult_genes$`nPeaksNearby_TGFb-and-RA-med`) / nrow(mult_genes)
supermult_genes <- filter(genesUp, `integrationCategory-med-dose` == 'super-multiplicative')
supermult_frac  <- sum(supermult_genes$`nPeaksNearby_TGFb-and-RA-med`) / nrow(supermult_genes)

set.seed(0)
results_vec <- c()
samplesize <- nrow(subadd_genes)
for (ii in 1:25000) {
  this.sample <- sample_n(genesUp, samplesize)
  results_vec <- c(results_vec, sum(this.sample$`nPeaksNearby_TGFb-and-RA-med`) / samplesize)
}

qplot(results_vec, bins = 100) + 
  geom_vline(xintercept = subadd_frac, color = "pink") +
  geom_vline(xintercept = add_frac, color = "red") +
  geom_vline(xintercept = mult_frac, color = "blue") +
  geom_vline(xintercept = supermult_frac, color = "dark blue") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of genes, number of peaks nearby") +
  xlab("Number of peaks nearby (med dose both signals)")

print(sprintf("sample size: %d", samplesize))
print(c(subadd_frac, add_frac, mult_frac, supermult_frac))
quantile(results_vec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))





