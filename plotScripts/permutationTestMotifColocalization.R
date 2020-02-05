library(tidyverse)
library(here)

joinedTib <- read_tsv(here('extractedData', 'joinedTablePeaksNearGenes.annotated.tsv'))
geneTib   <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTib   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))

n.permutations.to.make.null.distributions <- 12000

# 1. Analysis on peaks near upregulated genes
genesUp <- joinedTib %>% 
  filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>%
  filter(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0)


# now get the superadditive peaks, unique set (avoid duplicates near genes that are close together)
genesUpPeaksUp <- genesUp %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, peak_startLoc, peak_chrom, `isUpRegulatedPeak-med`, ensg) %>%
  unique()
  # filter(`peakAdditivePredFcResidual-med` > 1.5)


#idea: what kinds of motif matches enrich in the superadditive peaks near these genes? 334 superadditive peaks total...
# null distribution could be the rest of the peaks nearby that aren't superadditive (N = 4,313)
# to avoid double counting, could extract peak features only and unique-ify them
# can test specific motif frequencies, motif pair frequencies...
# first, test frequency of double motif match (all RA, all TGFb)

genesUpPeaksUp.nonSuperAdditive  <- filter(genesUpPeaksUp, ! `isSuperAdditivePeak-med`)
genesUpPeaksUp.superAdditive     <- filter(genesUpPeaksUp, `isSuperAdditivePeak-med`)

dualMotifMatchesNotsuperadditive <- (genesUpPeaksUp.nonSuperAdditive$`group-allTGFB_maxMotifMatchScore` > 1) & (genesUpPeaksUp.nonSuperAdditive$`group-allRA_maxMotifMatchScore` > 1)
dualMotifMatchesSuperadditive    <- (genesUpPeaksUp.superAdditive$`group-allTGFB_maxMotifMatchScore` > 1) & (genesUpPeaksUp.superAdditive$`group-allRA_maxMotifMatchScore` > 1)
dualMotifMatchesAll <- (genesUpPeaksUp$`group-allTGFB_maxMotifMatchScore` > 1) & (genesUpPeaksUp$`group-allRA_maxMotifMatchScore` > 1)

nonsuperadditive.average <- sum(dualMotifMatchesNotsuperadditive) / length(dualMotifMatchesNotsuperadditive)
superadditive.average   <- sum(dualMotifMatchesSuperadditive) / length(dualMotifMatchesSuperadditive)
print(sprintf("fraction dual match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction dual match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(dualMotifMatchesSuperadditive)

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(dualMotifMatchesAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 100) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks near upregulated genes,\nfraction that have dual motif matches") +
  xlab("Fraction of peaks with a TGFb motif match and an RA motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))



# 2. Analysis on all upregulated peaks
#get the superadditive peaks, unique set (avoid duplicates near genes that are close together)
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.nonSuperAdditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

dualMotifMatchesNotsuperadditive <- (peaksUp.nonSuperAdditive$`group-allTGFB_maxMotifMatchScore` > 1) & (peaksUp.nonSuperAdditive$`group-allRA_maxMotifMatchScore` > 1)
dualMotifMatchesSuperadditive    <- (peaksUp.superAdditive$`group-allTGFB_maxMotifMatchScore` > 1) & (peaksUp.superAdditive$`group-allRA_maxMotifMatchScore` > 1)
dualMotifMatchesAll <- (peaksUp$`group-allTGFB_maxMotifMatchScore` > 1) & (peaksUp$`group-allRA_maxMotifMatchScore` > 1)

nonsuperadditive.average <- sum(dualMotifMatchesNotsuperadditive) / length(dualMotifMatchesNotsuperadditive)
superadditive.average   <- sum(dualMotifMatchesSuperadditive) / length(dualMotifMatchesSuperadditive)
print(sprintf("fraction dual match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction dual match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(dualMotifMatchesSuperadditive)
# samplesize <- 350


set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(dualMotifMatchesAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have dual motif matches") +
  xlab("Fraction of peaks with a TGFb motif match and an RA motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))


nGenesRepSuperAdd    <- genesUpPeaksUp.superAdditive$ensg %>% unique() %>% length()
nSuperAddPeaksPerGene.nearUpregulatedWithSuperadditive <- nrow(genesUpPeaksUp.superAdditive) / nGenesRepSuperAdd
nGenesRep.nearUpregulated <- genesUpPeaksUp$ensg %>% unique() %>% length()
nSuperAddPeaksPerGene.nearAllUpregulated <- nrow(filter(genesUpPeaksUp, `isSuperAdditivePeak-med`)) / nGenesRep.nearUpregulated


# 3. Analysis on all upregulated peaks -- FOX motif test at superadditive peaks; hypothesis that FOX will be enriched
#get the superadditive peaks, unique set (avoid duplicates near genes that are close together)
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.notSuperadditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

foxMotifMatchNotsuperadditive <- (peaksUp.notSuperadditive$`group-FOX_maxMotifMatchScore` > 1) 
foxMotifMatchSuperadditive    <- (peaksUp.superAdditive$`group-FOX_maxMotifMatchScore` > 1)
foxMotifMatchAll              <- (peaksUp$`group-FOX_maxMotifMatchScore` > 1) 

nonsuperadditive.average <- sum(foxMotifMatchNotsuperadditive) / length(foxMotifMatchNotsuperadditive)
superadditive.average   <- sum(foxMotifMatchSuperadditive) / length(foxMotifMatchSuperadditive)
print(sprintf("fraction FOX match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction FOX match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(foxMotifMatchSuperadditive)
# samplesize <- 350

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(foxMotifMatchAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have at least one FOX motif match") +
  xlab("Fraction of peaks with a FOX motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))







# 4. Analysis on all upregulated peaks -- AP1  motif test at superadditive peaks
#get the superadditive peaks, unique set (avoid duplicates near genes that are close together)
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.notSuperadditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

ap1MotifMatchNotsuperadditive <- (peaksUp.notSuperadditive$`group-AP1_maxMotifMatchScore` > 1) 
ap1MotifMatchSuperadditive    <- (peaksUp.superAdditive$`group-AP1_maxMotifMatchScore` > 1) 
ap1MotifMatchAll              <- (peaksUp$`group-AP1_maxMotifMatchScore` > 1) 

nonsuperadditive.average <- sum(ap1MotifMatchNotsuperadditive) / length(ap1MotifMatchNotsuperadditive)
superadditive.average   <- sum(ap1MotifMatchSuperadditive) / length(ap1MotifMatchSuperadditive)
print(sprintf("fraction ap1 match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction ap1 match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(ap1MotifMatchSuperadditive)
# samplesize <- 350

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(ap1MotifMatchAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have at least one AP1 motif match") +
  xlab("Fraction of peaks with a AP1 motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))






# 5. Analysis on all upregulated peaks -- SMAD
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.notSuperadditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

smadMotifMatchNotsuperadditive <- (peaksUp.notSuperadditive$`group-SMAD_maxMotifMatchScore` > 1) 
smadMotifMatchSuperadditive    <- (peaksUp.superAdditive$`group-SMAD_maxMotifMatchScore` > 1) 
smadMotifMatchAll              <- (peaksUp$`group-SMAD_maxMotifMatchScore` > 1) 
  
nonsuperadditive.average <- sum(smadMotifMatchNotsuperadditive) / length(smadMotifMatchNotsuperadditive)
superadditive.average   <- sum(smadMotifMatchSuperadditive) / length(smadMotifMatchSuperadditive)
print(sprintf("fraction smad match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction smad match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(smadMotifMatchSuperadditive)
# samplesize <- 350

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(smadMotifMatchAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have at least one smad motif match") +
  xlab("Fraction of peaks with a smad motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))







# 6. Analysis on all upregulated peaks -- RAR
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.notSuperadditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

rarMotifMatchNotsuperadditive <- (peaksUp.notSuperadditive$`group-RAR_maxMotifMatchScore` > 1) 
rarMotifMatchSuperadditive    <- (peaksUp.superAdditive$`group-RAR_maxMotifMatchScore` > 1) 
rarMotifMatchAll              <- (peaksUp$`group-RAR_maxMotifMatchScore` > 1) 

nonsuperadditive.average <- sum(rarMotifMatchNotsuperadditive) / length(rarMotifMatchNotsuperadditive)
superadditive.average   <- sum(rarMotifMatchSuperadditive) / length(rarMotifMatchSuperadditive)
print(sprintf("fraction rar match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction rar match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(rarMotifMatchSuperadditive)
# samplesize <- 350

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(rarMotifMatchAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have at least one RAR motif match") +
  xlab("Fraction of peaks with a RAR motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))






# 7. Analysis on all upregulated peaks -- AP1 + FOX motif test at superadditive peaks; hypothesis that both will be enriched
#get the superadditive peaks, unique set (avoid duplicates near genes that are close together)
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.notSuperadditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

ap1AndFoxMotifMatchNotsuperadditive <- (peaksUp.notSuperadditive$`group-AP1_maxMotifMatchScore` > 1) & (peaksUp.notSuperadditive$`group-FOX_maxMotifMatchScore` > 1) 
ap1AndFoxMotifMatchSuperadditive    <- (peaksUp.superAdditive$`group-AP1_maxMotifMatchScore` > 1) & (peaksUp.superAdditive$`group-FOX_maxMotifMatchScore` > 1) 
ap1AndFoxMotifMatchAll              <- (peaksUp$`group-AP1_maxMotifMatchScore` > 1) & (peaksUp$`group-FOX_maxMotifMatchScore` > 1) 


nonsuperadditive.average <- sum(ap1AndFoxMotifMatchNotsuperadditive) / length(ap1AndFoxMotifMatchNotsuperadditive)
superadditive.average   <- sum(ap1AndFoxMotifMatchSuperadditive) / length(ap1AndFoxMotifMatchSuperadditive)
print(sprintf("fraction ap1AndFox match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction ap1AndFox match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(ap1AndFoxMotifMatchSuperadditive)
# samplesize <- 350

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(ap1AndFoxMotifMatchAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have at least one ap1AndFox motif match") +
  xlab("Fraction of peaks with a ap1AndFox motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))





# 8. Analysis on all upregulated peaks -- SMAD + FOX motif test at superadditive peaks; hypothesis that both will be enriched
#get the superadditive peaks, unique set (avoid duplicates near genes that are close together)
peaksUp <- peakTib %>%
  mutate(`isSuperAdditivePeak-med` = `peakAdditivePredFcResidual-med` > 1.5) %>%
  mutate(`isUpRegulatedPeak-med` = `TGFb-and-RA-med-avgFoldchange` > 1) %>%
  filter(`TGFb-med-avgFoldchange` > 1) %>%
  filter(`RA-med-avgFoldchange` > 1) %>%
  filter(`isUpRegulatedPeak-med`) %>%
  dplyr::select(matches(".*otif.*"), `isSuperAdditivePeak-med`, startLocs, chrom, `isUpRegulatedPeak-med`) %>%
  unique()
# filter(`peakAdditivePredFcResidual-med` > 1.5)

peaksUp.notSuperadditive  <- filter(peaksUp, ! `isSuperAdditivePeak-med`)
peaksUp.superAdditive     <- filter(peaksUp, `isSuperAdditivePeak-med`)

smadAndFoxMotifMatchNotsuperadditive <- (peaksUp.notSuperadditive$`group-SMAD_maxMotifMatchScore` > 1) & (peaksUp.notSuperadditive$`group-FOX_maxMotifMatchScore` > 1) 
smadAndFoxMotifMatchSuperadditive    <- (peaksUp.superAdditive$`group-SMAD_maxMotifMatchScore` > 1) & (peaksUp.superAdditive$`group-FOX_maxMotifMatchScore` > 1) 
smadAndFoxMotifMatchAll              <- (peaksUp$`group-SMAD_maxMotifMatchScore` > 1) & (peaksUp$`group-FOX_maxMotifMatchScore` > 1) 


nonsuperadditive.average <- sum(smadAndFoxMotifMatchNotsuperadditive) / length(smadAndFoxMotifMatchNotsuperadditive)
superadditive.average   <- sum(smadAndFoxMotifMatchSuperadditive) / length(smadAndFoxMotifMatchSuperadditive)
print(sprintf("fraction smadAndFox match, non-superadditive: %0.3f", nonsuperadditive.average))
print(sprintf("fraction smadAndFox match, superadditive: %0.3f", superadditive.average))

# permutation test for significance

samplesize <- length(smadMotifMatchSuperadditive)
# samplesize <- 350

set.seed(0)
resultsvec <- c()
for (ii in 1:n.permutations.to.make.null.distributions) {
  this.sample <- sample(smadAndFoxMotifMatchAll, samplesize)
  fracDual <- sum(this.sample) / samplesize
  resultsvec <- c(resultsvec, fracDual)
}

qplot(resultsvec, bins = 500) +
  geom_vline(xintercept = superadditive.average, color = "red") +
  theme_minimal(base_size = 16) +
  ggtitle("Randomly sampled groups of peaks upregulated by both signals,\nfraction that have at least one smadAndFox motif match") +
  xlab("Fraction of peaks with a smadAndFox motif match")

quantile(resultsvec, c(0, .001, 0.01, 0.025, .25, .50, .75, .975, 0.99, .999, 1))



