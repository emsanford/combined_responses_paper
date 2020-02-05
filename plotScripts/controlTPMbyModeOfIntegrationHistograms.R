library(tidyverse)
library(here)

geneTibAnno   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))

filter(deseqTibAnno, `RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`)

# filter out diff peaks for which to do signal integration analysis on
siUpregGenes <- geneTibAnno %>% filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`,
                                       `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`,
                                       `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0)  


tpm.int.subset <- siUpregGenes[, c("integrationCategory-low-dose", "integrationCategory-med-dose", "integrationCategory-high-dose", "EtOH-nlDensity_avgTPM")]
ggplot(tpm.int.subset, aes(`EtOH-nlDensity_avgTPM`)) + geom_histogram() + facet_wrap(~`integrationCategory-low-dose`) + ggtitle("low dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 50) + xlab("EtOH control TPM value") 
ggplot(tpm.int.subset, aes(`EtOH-nlDensity_avgTPM`)) + geom_histogram() + facet_wrap(~`integrationCategory-med-dose`) + ggtitle("med dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 50) + xlab("EtOH control TPM value")
ggplot(tpm.int.subset, aes(`EtOH-nlDensity_avgTPM`)) + geom_histogram() + facet_wrap(~`integrationCategory-high-dose`) + ggtitle("high dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 50) + xlab("EtOH control TPM value")


# do the same analysis for peaks

peakTibAnno   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))

# filter out diff peaks for which to do signal integration analysis on
siUpregPeaks <- peakTibAnno %>% filter(`RA-low-isDiffPeak` | `RA-med-isDiffPeak` | `RA-high-isDiffPeak`,
                                       `TGFb-low-isDiffPeak` | `TGFb-med-isDiffPeak` | `TGFb-high-isDiffPeak`,
                                       `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1,
                                       !`peak_integrationCategory-low-dose` %in% c("uncategorized", "super-additive", "sub-multiplicative"),
                                       !`peak_integrationCategory-med-dose` %in% c("uncategorized", "super-additive", "sub-multiplicative"),
                                       !`peak_integrationCategory-high-dose`%in% c("uncategorized", "super-additive", "sub-multiplicative"))
ggplot(siUpregPeaks, aes(`EtOH-nlDensity-avgNormFragmentCounts`)) + geom_histogram() + facet_wrap(~`peak_integrationCategory-low-dose`) + ggtitle("low dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 300) + xlab("EtOH control avgNormFragmentCounts") 
ggplot(siUpregPeaks, aes(`EtOH-nlDensity-avgNormFragmentCounts`)) + geom_histogram() + facet_wrap(~`peak_integrationCategory-med-dose`) + ggtitle("med dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 300) + xlab("EtOH control avgNormFragmentCounts")
ggplot(siUpregPeaks, aes(`EtOH-nlDensity-avgNormFragmentCounts`)) + geom_histogram() + facet_wrap(~`peak_integrationCategory-high-dose`) + ggtitle("high dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 300) + xlab("EtOH control avgNormFragmentCounts")
ggplot(siUpregPeaks, aes(`EtOH-nlDensity-avgNormFragmentCounts`)) + geom_density() + facet_wrap(~`peak_integrationCategory-low-dose`) + ggtitle("low dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 300) + xlab("EtOH control avgNormFragmentCounts") 
ggplot(siUpregPeaks, aes(`EtOH-nlDensity-avgNormFragmentCounts`)) + geom_density() + facet_wrap(~`peak_integrationCategory-med-dose`) + ggtitle("med dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 300) + xlab("EtOH control avgNormFragmentCounts")
ggplot(siUpregPeaks, aes(`EtOH-nlDensity-avgNormFragmentCounts`)) + geom_density() + facet_wrap(~`peak_integrationCategory-high-dose`) + ggtitle("high dose, upregulated genes") + theme_minimal(base_size = 16) + xlim(0, 300) + xlab("EtOH control avgNormFragmentCounts")

# how many peaks will we have left if we toy with thresholds to make them more stringent?
pmax(siUpregPeaks$`TGFb-high-avgNormFragmentCounts`, siUpregPeaks$`RA-high-avgNormFragmentCounts`) %>% qplot()
pmin(siUpregPeaks$`TGFb-high-avgNormFragmentCounts`, siUpregPeaks$`RA-high-avgNormFragmentCounts`) %>% qplot()

minAvgNormFragCounts <- 20
minFoldChange        <- 1.5
sum(((pmax(siUpregPeaks$`TGFb-high-avgNormFragmentCounts`, siUpregPeaks$`RA-high-avgNormFragmentCounts`) ) > minAvgNormFragCounts) & (pmax(siUpregPeaks$`TGFb-high-avgFoldchange`, siUpregPeaks$`RA-high-avgFoldchange`) > minFoldChange))
sum(((pmin(siUpregPeaks$`TGFb-high-avgNormFragmentCounts`, siUpregPeaks$`RA-high-avgNormFragmentCounts`) ) > minAvgNormFragCounts) & (pmin(siUpregPeaks$`TGFb-high-avgFoldchange`, siUpregPeaks$`RA-high-avgFoldchange`) > minFoldChange))

weakpeaks <- filter(siUpregPeaks, `TGFb-high-avgNormFragmentCounts` < minAvgNormFragCounts, siUpregPeaks$`TGFb-high-avgFoldchange` < minFoldChange)


