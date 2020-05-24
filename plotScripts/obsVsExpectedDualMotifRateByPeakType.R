library(tidyverse)
library(here)
source(here('extractionScripts', 'util.R'))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  upreg.peaks          <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  outputPlotPrefix        <- here("plots", "")
} else {
  upreg.peaks       <- read_tsv(cmdargs[1])
  outputPlotPrefix  <- cmdargs[2]
}

ra.dominant.factors   <- c("RARA", 'FOXA1', 'FOXA2', 'FOXA3', 'FOXC2', 'FOXD3', 'SPI', 'SPIB', "SPIC", "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5")
tgfb.dominant.factors <- c("SMAD3", "SMAD4", "SMAD9", "JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", "BACH1", "BACH2", "BATF",
                           "SMARCC1", "NFE2", "NFE2L2", "MAFF", "MAFK")

n.random.permutations  <- 1000
n.bootstrap.samples    <- 1000
peaktypes              <- sapply(upreg.peaks$`peak_integrationCategory-med-dose`, convertUpregCvalCatToDvalCat)
num.motif.match.matrix <- upreg.peaks %>% dplyr::select(matches("numMotifMatches"))

# build list
peak.motif.match.tib <- NULL

for (ii in 1:nrow(num.motif.match.matrix)) {
  peak.id   <- paste0(upreg.peaks[ii, "chrom"], ":", upreg.peaks[ii, "startLocs"],"-", upreg.peaks[ii, "endLocs"]) 
  peak.type <- peaktypes[ii]
  peak.motif.matches <- c()
  is.ra.motif        <- c()
  is.tgfb.motif      <- c()
  for (jj in 1:ncol(num.motif.match.matrix)) {
    motif.name <- strsplit(colnames(num.motif.match.matrix)[jj], "_", 1)[[1]][1]
    is.RA.dominant.motif   <- motif.name %in% ra.dominant.factors
    is.TGFb.dominant.motif <- motif.name %in% tgfb.dominant.factors
    peak.motif.matches <- c(peak.motif.matches, rep(motif.name,             num.motif.match.matrix[ii, jj]))
    is.ra.motif        <- c(is.ra.motif,        rep(is.RA.dominant.motif,   num.motif.match.matrix[ii, jj]))
    is.tgfb.motif      <- c(is.tgfb.motif,      rep(is.TGFb.dominant.motif, num.motif.match.matrix[ii, jj]))
  }
  peak.motif.match.tib <- rbind(peak.motif.match.tib, tibble(peak.id = peak.id, peak.type = peak.type, motif.name = peak.motif.matches, is.ra.motif = is.ra.motif, is.tgfb.motif = is.tgfb.motif))
}

# next, reduce related motifs into "groups"--e.g. since the related AP-1 motifs have similar motifs, we don't want to treat them as individual occurences
# but rather one observation. so instead of 1 jun, 1 fos, 1 jdp, and 1 fosl motif, we would just have 1 ap-1 motif.

# new tibble has one row per motif group match rather than one row per motif match
# retinoic acid factor groups
rar.factors <- c("RARA")
fox.factors <- c('FOXA1', 'FOXA2', 'FOXA3', 'FOXC2', 'FOXD3')
ets.factors <- c('SPI', 'SPIB', "SPIC", "EHF", "ELF1", "ELF2", "ELF3", "ELF4", "ELF5")
ra.factor.groups <- list(rar.factors, fox.factors, ets.factors)
ra.factor.group.names <- c("RARA", "FOX", "ETS")
#tgf-b factor groups
smad.factors <- c("SMAD3", "SMAD4", "SMAD9")
ap1.and.related.factors <- c("JUN", "JUNB", "JUND", "JDP2", "FOS", "FOSB", "FOSL1", "FOSL2", "BACH1", "BACH2", "BATF")
smarcc1.factor <- c("SMARCC1")
nfe.factors <- c("NFE2", "NFE2L2")
maf.factors <- c("MAFF", "MAFK")
tgfb.factor.groups <- list(smad.factors, ap1.and.related.factors, smarcc1.factor, nfe.factors, maf.factors)
tgfb.factor.group.names <- c("SMAD", "AP1", "SMARCC1", "NFE2", "MAF")
#neither factor groups
hox.factors  <- c("HOXA13", 'HOXB13', 'HOXC10', 'HOXC12', 'HOXC13', 'HOXD13')
nfkb.factors <- c("NFKB1", 'REL', 'RELA')
cdx.factors  <- c("CDX1", "CDX2")
ctcf.factors <- c("CTCF")
bcl.factors  <- c("BCL11A", "BCL11B")
grhl.factors <- c("GRHL1")
other.factor.groups      <- list(hox.factors, nfkb.factors, cdx.factors, ctcf.factors, bcl.factors, grhl.factors)
other.factor.group.names <- list("HOX", "NFKB", "CDX", "CTCF", "BCL", "GRHL1")

countMaxMotifCountInGroup <- function(motif.vector, motifs.in.group) {
  group.related.motifs <- motif.vector[motif.vector %in% motifs.in.group]
  if (length(group.related.motifs) > 0) {
    return(max(table(group.related.motifs)))
  } else {
    return(0)
  }
}

grouped.motif.tib <- NULL
peaks.with.motifs <- peak.motif.match.tib$peak.id %>% unique()
for (peakid in peaks.with.motifs) {
  this.peak.tib  <- peak.motif.match.tib %>% filter(peak.id == peakid)
  this.peak.type <- this.peak.tib$peak.type[1]
  this.peak.motifs <- this.peak.tib$motif.name
  
  for (ii in 1:length(ra.factor.groups)) {
    n.this.group <- countMaxMotifCountInGroup(this.peak.motifs, ra.factor.groups[[ii]])
    grouped.motif.tib <- rbind(grouped.motif.tib, tibble(peak.id = peakid, peak.type = this.peak.type, motif.group = rep(ra.factor.group.names[[ii]], n.this.group), is.ra.motif = TRUE, is.tgfb.motif = FALSE))
  }
  
  for (ii in 1:length(tgfb.factor.groups)) {
    n.this.group <- countMaxMotifCountInGroup(this.peak.motifs, tgfb.factor.groups[[ii]])
    grouped.motif.tib <- rbind(grouped.motif.tib, tibble(peak.id = peakid, peak.type = this.peak.type, motif.group = rep(tgfb.factor.group.names[[ii]], n.this.group), is.ra.motif = FALSE, is.tgfb.motif = TRUE))
  }
  
  for (ii in 1:length(other.factor.groups)) {
    n.this.group <- countMaxMotifCountInGroup(this.peak.motifs, other.factor.groups[[ii]])
    grouped.motif.tib <- rbind(grouped.motif.tib, tibble(peak.id = peakid, peak.type = this.peak.type, motif.group = rep(other.factor.group.names[[ii]], n.this.group), is.ra.motif = FALSE, is.tgfb.motif = FALSE))
  }
}

calculateNumDualMotifMatchRate <- function(motif.tib) {
  filt.tib.one.row.per.peak <- motif.tib %>%
    group_by(peak.id) %>%
    mutate(peak.has.dual.motif = any(is.ra.motif) & any(is.tgfb.motif)) %>%
    dplyr::select(peak.id, peak.type, peak.has.dual.motif) %>%
    unique()
  n.dual.motifs <- sum(filt.tib.one.row.per.peak$peak.has.dual.motif)
  return(n.dual.motifs)
}

subadd.peak.grouped.motif.tib   <- grouped.motif.tib %>% filter(peak.type == "sub-additive")
add.peak.grouped.motif.tib      <- grouped.motif.tib %>% filter(peak.type == "additive")
superadd.peak.grouped.motif.tib <- grouped.motif.tib %>% filter(peak.type == "super-additive")

n.subadd.peaks <- sum(peaktypes == 'sub-additive')
n.add.peaks <- sum(peaktypes == 'additive')
n.superadd.peaks <- sum(peaktypes == 'super-additive')

subadd.permutation.dual.motif.match.rates   <- c()
add.permutation.dual.motif.match.rates      <- c()
superadd.permutation.dual.motif.match.rates <- c()

set.seed(0)
for (ii in 1:n.random.permutations) {
  permuted.subadd.peak.grouped.motif.tib        <- subadd.peak.grouped.motif.tib
  permuted.subadd.peak.grouped.motif.tib[, 3:5] <- subadd.peak.grouped.motif.tib[sample(nrow(subadd.peak.grouped.motif.tib)), 3:5]
  subadd.permutation.dual.motif.match.rates     <- c(subadd.permutation.dual.motif.match.rates, calculateNumDualMotifMatchRate(permuted.subadd.peak.grouped.motif.tib) / n.subadd.peaks)
  
  permuted.add.peak.grouped.motif.tib        <- add.peak.grouped.motif.tib
  permuted.add.peak.grouped.motif.tib[, 3:5] <- add.peak.grouped.motif.tib[sample(nrow(add.peak.grouped.motif.tib)), 3:5]
  add.permutation.dual.motif.match.rates     <- c(add.permutation.dual.motif.match.rates, calculateNumDualMotifMatchRate(permuted.add.peak.grouped.motif.tib) / n.add.peaks)
  
  permuted.superadd.peak.grouped.motif.tib        <- superadd.peak.grouped.motif.tib
  permuted.superadd.peak.grouped.motif.tib[, 3:5] <- superadd.peak.grouped.motif.tib[sample(nrow(superadd.peak.grouped.motif.tib)), 3:5]
  superadd.permutation.dual.motif.match.rates     <- c(superadd.permutation.dual.motif.match.rates, calculateNumDualMotifMatchRate(permuted.superadd.peak.grouped.motif.tib) / n.superadd.peaks)
}
obs.rate.subadd   <- calculateNumDualMotifMatchRate(subadd.peak.grouped.motif.tib)   / n.subadd.peaks
obs.rate.add      <- calculateNumDualMotifMatchRate(add.peak.grouped.motif.tib)      / n.add.peaks
obs.rate.superadd <- calculateNumDualMotifMatchRate(superadd.peak.grouped.motif.tib) / n.superadd.peaks

mean(subadd.permutation.dual.motif.match.rates) - obs.rate.subadd
mean(add.permutation.dual.motif.match.rates) - obs.rate.add
mean(superadd.permutation.dual.motif.match.rates) - obs.rate.superadd

obs.rate.subadd / mean(subadd.permutation.dual.motif.match.rates)
obs.rate.add / mean(add.permutation.dual.motif.match.rates)
obs.rate.superadd / mean(superadd.permutation.dual.motif.match.rates)

# library(patchwork)
# p1 <- qplot(subadd.permutation.dual.motif.match.rates, binwidth = 0.001) + xlim(.15, .55) + geom_vline(xintercept = obs.rate.subadd) + ggtitle("sub-additive peaks: expected distribution vs observed (vertical line) rate of dual-motif matches")
# p2 <- qplot(add.permutation.dual.motif.match.rates, binwidth = 0.001) + xlim(.15, .55) + geom_vline(xintercept = obs.rate.add) + ggtitle("additive peaks: expected distribution vs observed (vertical line) rate of dual-motif matches")
# p3 <- qplot(superadd.permutation.dual.motif.match.rates, binwidth = 0.001) + xlim(.15, .55) + geom_vline(xintercept = obs.rate.superadd) + ggtitle("super-additive peaks: expected distribution vs observed (vertical line) rate of dual-motif matches")
# p1 / p2 / p3

### are super-additive peaks more likely to contain dual motifs?
getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks <- function(peak.tib.anno, peak.integration.category, peak.dose) {
  int.categories <- sapply(pull(peak.tib.anno, paste0("peak_integrationCategory-", peak.dose, "-dose")), convertUpregCvalCatToDvalCat)
  int.cat.inds   <- int.categories == peak.integration.category
  filt.peak.tib.anno <- peak.tib.anno[int.cat.inds, ]
  
  dual.motif.tib <- filt.peak.tib.anno %>%
    mutate(hasTGFbMatch = `group-TGFbdominant_maxMotifMatchScore` > 0,
           hasRAMatch   = `group-RAdominant_maxMotifMatchScore` > 0) %>%
    mutate(hasDualMotifMatch = hasTGFbMatch & hasRAMatch)
  
  n.peaks.this.cat <- nrow(dual.motif.tib)
  measured_frac_dual_motif = sum(dual.motif.tib$hasDualMotifMatch) / n.peaks.this.cat

  return(measured_frac_dual_motif)
}

frac.dual.motif.matches.measured.subadditive   <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(upreg.peaks, "sub-additive", "med")
frac.dual.motif.matches.measured.additive      <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(upreg.peaks, "additive", "med")
frac.dual.motif.matches.measured.superadditive <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(upreg.peaks, "super-additive", "med")

# do bootstrap right here, add CI's
set.seed(0)
frac.dual.motif.matches.measured.subadditive.bootstrap.values <- c()
frac.dual.motif.matches.measured.additive.bootstrap.values      <- c()
frac.dual.motif.matches.measured.superadditive.bootstrap.values <- c()

for (ii in 1:n.bootstrap.samples) {
  upregPeakBootstrapSample <- sample_n(upreg.peaks, nrow(upreg.peaks), replace = TRUE)
  
  bootstrap.frac.dual.motif.matches.measured.subadditive   <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(upregPeakBootstrapSample, "sub-additive", "med") - frac.dual.motif.matches.measured.subadditive
  bootstrap.frac.dual.motif.matches.measured.additive      <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(upregPeakBootstrapSample, "additive", "med") - frac.dual.motif.matches.measured.additive
  bootstrap.frac.dual.motif.matches.measured.superadditive <- getExpectedAndMeasuredDualMotifMatchesAtAnnotatedPeaks(upregPeakBootstrapSample, "super-additive", "med") - frac.dual.motif.matches.measured.superadditive

  frac.dual.motif.matches.measured.subadditive.bootstrap.values   <- c(frac.dual.motif.matches.measured.subadditive.bootstrap.values, bootstrap.frac.dual.motif.matches.measured.subadditive)
  frac.dual.motif.matches.measured.additive.bootstrap.values      <- c(frac.dual.motif.matches.measured.additive.bootstrap.values, bootstrap.frac.dual.motif.matches.measured.additive)
  frac.dual.motif.matches.measured.superadditive.bootstrap.values <- c(frac.dual.motif.matches.measured.superadditive.bootstrap.values, bootstrap.frac.dual.motif.matches.measured.superadditive)
}
# now build tibble for plot, row-by-row
reduced.dual.motif.analysis.tib <- tibble(expected.vs.measured = c("measured","expected", "measured","expected", "measured", "expected"),
                                          peak.category        = factor(c('sub-additive', 'sub-additive', 'additive', 'additive', 'super-additive', 'super-additive'), levels = c("sub-additive", "additive", "super-additive")),
                                          frac.dual.motifs = c(frac.dual.motif.matches.measured.subadditive, mean(subadd.permutation.dual.motif.match.rates),
                                                               frac.dual.motif.matches.measured.additive, mean(add.permutation.dual.motif.match.rates),
                                                               frac.dual.motif.matches.measured.superadditive, mean(superadd.permutation.dual.motif.match.rates)),
                                          upper.ci = c(frac.dual.motif.matches.measured.subadditive - quantile(frac.dual.motif.matches.measured.subadditive.bootstrap.values, .05), 
                                                       quantile(subadd.permutation.dual.motif.match.rates, .95),
                                                       frac.dual.motif.matches.measured.additive - quantile(frac.dual.motif.matches.measured.additive.bootstrap.values, .05), 
                                                       quantile(add.permutation.dual.motif.match.rates, .95),
                                                       frac.dual.motif.matches.measured.superadditive - quantile(frac.dual.motif.matches.measured.superadditive.bootstrap.values, .05), 
                                                       quantile(superadd.permutation.dual.motif.match.rates, .95)),
                                          lower.ci = c(frac.dual.motif.matches.measured.subadditive - quantile(frac.dual.motif.matches.measured.subadditive.bootstrap.values, .95), 
                                                       quantile(subadd.permutation.dual.motif.match.rates, .05),
                                                       frac.dual.motif.matches.measured.additive - quantile(frac.dual.motif.matches.measured.additive.bootstrap.values, .95), 
                                                       quantile(add.permutation.dual.motif.match.rates, .05),
                                                       frac.dual.motif.matches.measured.superadditive - quantile(frac.dual.motif.matches.measured.superadditive.bootstrap.values, .95), 
                                                       quantile(superadd.permutation.dual.motif.match.rates, .05))
)

dual.motif.analysis.plot <- reduced.dual.motif.analysis.tib %>%
  ggplot(aes(x = peak.category, y = frac.dual.motifs, fill = expected.vs.measured, ymin = lower.ci, ymax = upper.ci)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(position = position_dodge(width=0.9), width = 0) +
  theme_classic()


ggsave(paste0(outputPlotPrefix, "motif_analysis_freq_dual_motif_matches_by_peakIntCategory.svg"), plot = dual.motif.analysis.plot, width = 12, height = 12)

