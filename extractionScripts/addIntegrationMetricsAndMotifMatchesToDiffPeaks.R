library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)
library(chromVARmotifs)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(2019)
register(MulticoreParam(4, progressbar = TRUE))

data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper

## user-defined parameters
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  diffPeakTib    <- read_tsv(here('extractedData', 'differentialAtacPeaks.tsv'))
  dpOutputTibLoc <- here('extractedData', 'differentialAtacPeaks.annotated.tsv')
  selected.PWM.objects <- read_rds(here('extractedData', 'most_variable_motifs_Robject.rds'))
} else {
  diffPeakTib    <- read_tsv(cmdargs[1])
  selected.PWM.objects <- read_rds(cmdargs[2])
  dpOutputTibLoc <- cmdargs[3]
}

## function definitions
calcAdditivePeakMutualExclusivityScore <- function(etohCts, raCts, tgfbCts) {
  raEffectMagnitude   <- abs(raCts - etohCts)
  tgfbEffectMagnitude <- abs(tgfbCts - etohCts)
  lesserEffect <- min(raEffectMagnitude, tgfbEffectMagnitude)
  totalEffect  <- raEffectMagnitude + tgfbEffectMagnitude
  mutualExclusivityScore <- lesserEffect / totalEffect
  return(mutualExclusivityScore)
}

# the asymmetric score is close to 0 if a peak is TGFb-dominant and close to 1 if a peak is RA-dominant
calcAsymmmetricPeakMutualExclusivityScore <- function(etohCts, raCts, tgfbCts) {
  raEffectMagnitude   <- abs(raCts - etohCts)
  tgfbEffectMagnitude <- abs(tgfbCts - etohCts)
  totalEffect  <- raEffectMagnitude + tgfbEffectMagnitude
  mutualExclusivityScore <- raEffectMagnitude / totalEffect
  return(mutualExclusivityScore)
}

calcMultiplicativePeakMutualExclusivityScore <- function(etohCts, raCts, tgfbCts) {
  # add pseudocount to all if any of the conditions have zero counts. Oct-31-2019: leads to 132/83,091 peaks getting pseudocounts
  if (any(c(etohCts, raCts, tgfbCts) == 0)) {
    etohCts <- etohCts + 0.5
    raCts   <- raCts + 0.5
    tgfbCts <- tgfbCts + 0.5
    print("pseudocounts added to each peak in the mutual exclusivity score function")
  } 
  raFoldChange   <- raCts / etohCts
  tgfbFoldChange <- tgfbCts / etohCts
  raEffectMagnitude <- raFoldChange
  if (raEffectMagnitude < 1) {
    raEffectMagnitude <- 1 / raEffectMagnitude
  }
  tgfbEffectMagnitude <- tgfbFoldChange
  if (tgfbEffectMagnitude < 1) {
    tgfbEffectMagnitude <- 1 / tgfbEffectMagnitude
  }
  lesserEffect <- min(raEffectMagnitude, tgfbEffectMagnitude) - 1
  totalEffect  <- raEffectMagnitude + tgfbEffectMagnitude - 2         # effects in the fold-change space cannot be less than 1 and we want to capture the ratio of the relative effects
  mutualExclusivityScore <- lesserEffect / totalEffect
  return(mutualExclusivityScore)
}

addMutualExclusivityScoresToPeakTib <- function(peaktib) {
  addAvgMutualExclusivityScores      <- c()
  multAvgMutualExclusivityScores     <- c()
  asymmAddAvgMutualExclusivityScores <- c()
  for (ii in 1:length(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`)) {
    addlow  <- calcAdditivePeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                      peaktib$`RA-low-avgNormFragmentCounts`[ii], 
                                                      peaktib$`TGFb-low-avgNormFragmentCounts`[ii])
    addmed  <- calcAdditivePeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                      peaktib$`RA-med-avgNormFragmentCounts`[ii], 
                                                      peaktib$`TGFb-med-avgNormFragmentCounts`[ii])
    addhigh <- calcAdditivePeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                      peaktib$`RA-high-avgNormFragmentCounts`[ii], 
                                                      peaktib$`TGFb-high-avgNormFragmentCounts`[ii])
    addAvgMutualExclusivityScore <- (addlow + addmed + addhigh) / 3
    
    aaddlow  <- calcAsymmmetricPeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                          peaktib$`RA-low-avgNormFragmentCounts`[ii], 
                                                          peaktib$`TGFb-low-avgNormFragmentCounts`[ii])
    aaddmed  <- calcAsymmmetricPeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                          peaktib$`RA-med-avgNormFragmentCounts`[ii], 
                                                          peaktib$`TGFb-med-avgNormFragmentCounts`[ii])
    aaddhigh <- calcAsymmmetricPeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                          peaktib$`RA-high-avgNormFragmentCounts`[ii], 
                                                          peaktib$`TGFb-high-avgNormFragmentCounts`[ii])
    asymmAddAvgMutualExclusivityScore <- (aaddlow + aaddmed + aaddhigh) / 3
      
    multlow  <- calcMultiplicativePeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                             peaktib$`RA-low-avgNormFragmentCounts`[ii], 
                                                             peaktib$`TGFb-low-avgNormFragmentCounts`[ii])
    multmed  <- calcMultiplicativePeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                             peaktib$`RA-med-avgNormFragmentCounts`[ii], 
                                                             peaktib$`TGFb-med-avgNormFragmentCounts`[ii])
    multhigh <- calcMultiplicativePeakMutualExclusivityScore(peaktib$`EtOH-nlDensity-avgNormFragmentCounts`[ii], 
                                                             peaktib$`RA-high-avgNormFragmentCounts`[ii], 
                                                             peaktib$`TGFb-high-avgNormFragmentCounts`[ii])
    multAvgMutualExclusivityScore <- (multlow + multmed + multhigh) / 3
    
    addAvgMutualExclusivityScores  <- c(addAvgMutualExclusivityScores, addAvgMutualExclusivityScore)
    asymmAddAvgMutualExclusivityScores <- c(asymmAddAvgMutualExclusivityScores, asymmAddAvgMutualExclusivityScore)
    multAvgMutualExclusivityScores <- c(multAvgMutualExclusivityScores, multAvgMutualExclusivityScore)
  }
  
  peaktib[["PeakMutualExclusivityScoreAdditive"]]           <- addAvgMutualExclusivityScores
  peaktib[["PeakMutualExclusivityScoreAsymmetricAdditive"]] <- asymmAddAvgMutualExclusivityScores
  peaktib[["PeakMutualExclusivityScoreMultiplicative"]]     <- multAvgMutualExclusivityScores
  return(peaktib)
}

calcAdditivePredDiff <- function(peaktib, dosage) {
  normcts.ra    <- pull(peaktib, var = paste0("TGFb-", dosage, "-avgNormFragmentCounts"))
  normcts.tgfb  <- pull(peaktib, var = paste0("RA-", dosage, "-avgNormFragmentCounts"))
  normcts.etoh  <- pull(peaktib, var = `EtOH-nlDensity-avgNormFragmentCounts`)
  print(sprintf("Pseudocts added to %d peaks in the calc additive pred diff function", sum(normcts.etoh == 0)))
  normcts.etoh[normcts.etoh == 0] <- 0.5
  normcts.ra[normcts.ra == 0]     <- 0.5
  normcts.tgfb[normcts.tgfb == 0] <- 0.5
  print(sprintf("Pseudocts added to %d etoh peaks", sum(normcts.etoh == 0)))
  print(sprintf("Pseudocts added to %d ra peaks", sum(normcts.ra == 0)))
  print(sprintf("Pseudocts added to %d tgfb peaks", sum(normcts.tgfb == 0)))
  
  addpred.cts   <- normcts.ra + normcts.tgfb - normcts.etoh
  addpred.fc    <- addpred.cts / normcts.etoh 
  both.fc <- pull(peaktib, var = paste0("TGFb-and-RA-", dosage, "-avgNormFragmentCounts")) / normcts.etoh
  fc.pred.diff <- addpred.fc - both.fc 
  return(fc.pred.diff)
}

calcMultiplicativePredDiff <- function(peaktib, dosage) {
  normcts.ra    <- pull(peaktib, var = paste0("TGFb-", dosage, "-avgNormFragmentCounts"))
  normcts.tgfb  <- pull(peaktib, var = paste0("RA-", dosage, "-avgNormFragmentCounts"))
  normcts.etoh  <- pull(peaktib, var = `EtOH-nlDensity-avgNormFragmentCounts`)
  normcts.etoh[normcts.etoh == 0] <- 0.5
  normcts.ra[normcts.ra == 0]     <- 0.5
  normcts.tgfb[normcts.tgfb == 0] <- 0.5
  print(sprintf("Pseudocts added to %d etoh peaks", sum(normcts.etoh == 0)))
  print(sprintf("Pseudocts added to %d ra peaks", sum(normcts.ra == 0)))
  print(sprintf("Pseudocts added to %d tgfb peaks", sum(normcts.tgfb == 0)))
  
  multpred.cts   <- normcts.etoh * (normcts.ra / normcts.etoh) * (normcts.tgfb / normcts.etoh)
  multpred.fc    <- multpred.cts / normcts.etoh 
  both.fc <- pull(peaktib, var = paste0("TGFb-and-RA-", dosage, "-avgNormFragmentCounts")) / normcts.etoh
  fc.pred.diff <- multpred.fc - both.fc 
  return(fc.pred.diff)
}


# motif stuff...

## body of script

dpOutputTib <- addMutualExclusivityScoresToPeakTib(diffPeakTib)

for (dose in c("low", "med", "high")) {
  peak.add.pred.diff <- calcAdditivePredDiff(dpOutputTib, dose)
  dpOutputTib[[paste0('peakAdditivePredFcResidual-', dose)]] <- peak.add.pred.diff
  
  peak.mult.pred.diff <- calcMultiplicativePredDiff(dpOutputTib, dose)
  dpOutputTib[[paste0('peakMultiplicativePredFcResidual-', dose)]] <- peak.mult.pred.diff
}

diffpeak_gr <- GRanges(seqnames = dpOutputTib$chrom, 
                       ranges = IRanges(start = dpOutputTib$startLocs,
                                        end   = dpOutputTib$endLocs))

# note: details on match scores are here https://github.com/jhkorhonen/MOODS/wiki/Brief-theoretical-introduction
# "Intuitively, this score compares the probability that the model specified by the original PWM generated the sequence 
# versus the probability that the background model generated the sequence"
# note #2: the score sums the log probability ratio of PWM value over background value for each position in the motif. because of this,
# i think it's important to remember that longer motifs can obtain higher scores than shorter motifs.
selectedMotifsMatches <- matchMotifs(selected.PWM.objects, 
                                     diffpeak_gr, 
                                     genome = BSgenome.Hsapiens.UCSC.hg38, 
                                     bg = "genome",
                                     out = "scores") 
motif.names <- sapply(strsplit(colnames(selectedMotifsMatches), "_"), function (x) x[3])
counter <- 1
for (motif.name in motif.names) {
  dpOutputTib[[paste0(motif.name, '_motifMatchScore')]] <- motifScores(selectedMotifsMatches)[, counter]
  counter <- counter + 1
}

# now categorize motifs into sets and add motifsets
# these were the top 72 motifs selected by this script, which i then manually categorized into groups

# factors that go down in both RA and TGFb (i.e. EtOH factors)
klf.factors        <- c("KLF1", "KLF2", "KLF3", "KLF4", "KLF5", "KLF7")
ap2.factors        <- c("TFAP2B", "TFAP2C", "TFAP2E", "TFAP2A")
other.etoh.factors <- c("GRHL1", "SP3")
all.etoh.factors   <- c(klf.factors, ap2.factors, other.etoh.factors)

# NFKB factors, which go up a little in both RA and TGFb treatment
nfkb.factors <- c("REL", "RELA", "NFKB1")

# factors that go up in RA treatment
fox.factors      <- c("FOXJ2", "FOXA2", "FOXA1", "FOXA3", "FOXC2", "FOXD2", "FOXD3")
hox.factors      <- c("HOXC13", "HOXD13", "HOXB13", "HOXC10", "HOXA13")
rar.factors      <- c("RARA")
irf.factors.ra   <- c("IRF2", "IRF3", "IRF8")
elf.factors      <- c("ELF1", "ELF2", "ELF3", "ELF4", "ELF5")
other.ra.factors <- c("BCL11A", "BCL11B", "CDX1", "CDX2", "EHF", "GABPA", "NR2C2", "ETV7", "SPI1", "SPIB", "SPIC", "ELK4", "ETS2", "STAT2")
all.ra.factors   <- c(fox.factors, hox.factors, rar.factors, irf.factors.ra, elf.factors, other.ra.factors)

# factors that go up in TGFb treatment
ap1.factors        <- c("FOS", "FOSL1", "FOSL2", "FOSB", "JUN", "JUNB", "JUND", "JDP2", "BATF")
smad.factors       <- c("SMAD3", "SMAD4", "SMAD9")
irf.factors.tgfb   <- c("IRF4")
other.tgfb.factors <- c("BACH1", "BACH2", "SMARCC1", "NFE2", "NFE2L2", "MAFF", "MAFK")
all.tgfb.factors   <- c(ap1.factors, smad.factors, irf.factors.tgfb, other.tgfb.factors)

# other factor groups 
ctcf.factors <- c("CTCF", "CTCFL")  # slightly associated with RA but weakly associated with TGFb as well, previously found to be variable in controls

factor.group.list  <- list(klf.factors, ap2.factors, other.etoh.factors, all.etoh.factors, nfkb.factors, fox.factors, hox.factors, rar.factors, irf.factors.ra, elf.factors, other.ra.factors, all.ra.factors, ap1.factors, smad.factors, irf.factors.tgfb, other.tgfb.factors, all.tgfb.factors, ctcf.factors)
factor.group.names <-   c('group-KLF', 'group-AP2', 'group-otherEtOH',  'group-allEtOH',  'group-NFKB', 'group-FOX', 'group-HOX', 'group-RAR', 'group-raIRF',  'group-ELF', 'group-otherRA',  'group-allRA',  'group-AP1', 'group-SMAD', 'group-tgfbIRF',  'group-otherTGFB',  'group-allTGFB',  'groupCTCF')
stopifnot(length(factor.group.list) == length(factor.group.names))

for (ii in 1:length(factor.group.names)) {
  factor.group.name <- factor.group.names[ii]
  factor.group <- factor.group.list[[ii]]
  motif_selection_regex <- paste0("(", paste(factor.group, collapse="|"), ")", "_motifMatchScore")
  motifGroupMatchTib <- dplyr::select(dpOutputTib, matches(motif_selection_regex))
  motifGroupMaxScoresAsVector <- sapply(1:nrow(motifGroupMatchTib), function(x) max(as_vector(motifGroupMatchTib[x, ])))
  dpOutputTib[[paste0(factor.group.name, '_maxMotifMatchScore')]] <- motifGroupMaxScoresAsVector
}

#todo remove this
test2 <- dpOutputTib
dpOutputTib <- test2
# at this point, we do an almost identical categorial "peak integration type" calculation. most of this code is copy-pasted and lightly edited
# from the "addIntegrationMetricsToDeGenes.R" script.

# We assign a categorical variable to signal integration type
# parameters used for defining the range of estimated normal distribution for observed data.
#    additive/multiplicative predictions' locations with respect to the confidence interval 
#    of this distribution determine whether genes get classified as "additive", "multiplicative",
#    "subadditive", "supermultiplicative", "in-between", or "unable to tell"
tail.size <- 0.10
tail.test.lower <- 0 + tail.size
tail.test.upper <- 1 - tail.size
##########################################################################################

### helper functions for categorical signal integration assignment ###
rowVars <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

rowSse <- function(x) {
  rowSums((x - rowMeans(x))^2)
}

rowCVs_unbiased <- function(x) {
  N <- dim(x)[2]
  row.vars <- rowSums((x - rowMeans(x))^2)/(N - 1)
  row.stds <- sqrt(row.vars)
  row.cvs.biased  <- row.stds / rowMeans(x)
  row.cvs.unbiased <- (1 + 1/(4 * N)) * row.cvs.biased
  return(row.cvs.unbiased)
}

oneConditionEstMeanAndVarANFK <- function(cond, peak.tib, relevant.metadata) {
  filt.metadata <- filter(relevant.metadata, condition == cond)
  sample.names  <- filt.metadata$`condition-replicate`
  anfk.columns   <- peak.tib[, paste0(sample.names, "-normFragmentCounts")]
  mean_estimates <- rowMeans(anfk.columns)
  variance_estimates <- rowVars(anfk.columns)
  sum_squared_error <- rowSse(anfk.columns)
  cv_estimates <- rowCVs_unbiased(anfk.columns)
  tib.estimates <- tibble(peakID = paste0(peak.tib$chrom, peak.tib$startLocs, "-", peak.tib$endLocs))
  tib.estimates[[paste0(cond, "_peak_meanEstANFK")]] <- mean_estimates
  tib.estimates[[paste0(cond, "_peak_varEstANFK")]]  <- variance_estimates
  tib.estimates[[paste0(cond, "_peak_sseANFK")]]     <- sum_squared_error
  tib.estimates[[paste0(cond, "_peak_cvEstANFK")]]   <- cv_estimates
  return(tib.estimates)
}

determineCombEffectDirection <- function(a_effect, b_effect) {
  if (all(a_effect > 0, b_effect > 0)) {
    return("both up")
  } else if (xor(a_effect > 0, b_effect > 0)) {
    return("opposing")
  } else if (all(a_effect < 0, b_effect < 0)) {
    return("both down")
  } else {
    print(sprintf("no effect direction categorization. A effect is %.3f, B effect is %.3f", a_effect, b_effect))
    return("other")
  }
}

boundfxBothUp <- function(lowerbound, upperbound, addpred, multpred) {
  # print(paste0(lowerbound,"  ", upperbound))
  stopifnot(lowerbound <= upperbound)
  addPredBelow    <- addpred < lowerbound 
  addPredBetween  <- between(addpred, lowerbound, upperbound)
  addPredAbove    <- addpred  > upperbound 
  multPredBelow   <- multpred < lowerbound 
  multPredBetween <- between(multpred, lowerbound, upperbound)
  multPredAbove   <- multpred > upperbound 
  if (addPredBelow & multPredBelow) {
    return("super-multiplicative")
  } else if (addPredBelow & multPredBetween) {
    return("multiplicative")
  } else if (addPredBetween & multPredBetween) {
    return("ambiguous") 
  } else if (addPredBetween & multPredAbove) {
    return("additive")
  } else if (addPredAbove & multPredAbove) {
    return("sub-additive")
  } else if (addPredBelow & multPredAbove) {
    return("between-add-and-mult")
  } else {
    return("uncategorized")
  }
}

boundfxBothDown <- function(lowerbound, upperbound, addpred, multpred) {
  # print(paste0(lowerbound,"  ", upperbound))
  stopifnot(lowerbound <= upperbound)
  addPredBelow    <- addpred < lowerbound 
  addPredBetween  <- between(addpred, lowerbound, upperbound)
  addPredAbove    <- addpred  > upperbound 
  multPredBelow   <- multpred < lowerbound 
  multPredBetween <- between(multpred, lowerbound, upperbound)
  multPredAbove   <- multpred > upperbound 
  if (addPredBelow & multPredBelow) {
    return("sub-multiplicative")
  } else if (addPredBelow & multPredBetween) {
    return("multiplicative")
  } else if (addPredBetween & multPredBetween) {
    return("ambiguous") 
  } else if (addPredBetween & multPredAbove) {
    return("additive")
  } else if (addPredAbove & multPredAbove) {
    return("super-additive")
  } else if (addPredBelow & multPredAbove) {
    return("between-add-and-mult")
  } else {
    return("uncategorized")
  }
}

boundfxOpposing <- function(lowerbound, upperbound, addpred, multpred) {
  # print(paste0(lowerbound, "  ", upperbound))
  stopifnot(lowerbound <= upperbound)
  addPredBelow    <- addpred < lowerbound 
  addPredBetween  <- between(addpred, lowerbound, upperbound)
  addPredAbove    <- addpred  > upperbound 
  multPredBelow   <- multpred < lowerbound 
  multPredBetween <- between(multpred, lowerbound, upperbound)
  multPredAbove   <- multpred > upperbound 
  if (addPredBelow & multPredBelow) {
    return("super-additive")
  } else if (multPredBelow & addPredBetween) {
    return("additive")
  } else if (addPredBetween & multPredBetween) {
    return("ambiguous") 
  } else if (multPredBetween & addPredAbove) {
    return("multiplicative")
  } else if (addPredAbove & multPredAbove) {
    return("sub-multiplicative")
  } else if (multPredBelow & addPredAbove) {
    return("between-add-and-mult")
  } else {
    return("uncategorized")
  }
}
#################################################################################

# add additive model and multiplicative model predictions to output tibble. also add integrationConstant,
#   which is the value of the coefficient of the add/mult difference term that fits the observed data.
dosages <- c("low", "med", "high")
for (dosage in dosages) {
  etohAvgNormFragmentCounts <- dpOutputTib[["EtOH-nlDensity-avgNormFragmentCounts"]]
  raAvgNormFragmentCounts   <- dpOutputTib[[paste0("RA-", dosage, "-avgNormFragmentCounts")]]
  tgfbAvgNormFragmentCounts <- dpOutputTib[[paste0("TGFb-", dosage, "-avgNormFragmentCounts")]]
  bothAvgNormFragmentCounts <- dpOutputTib[[paste0("TGFb-and-RA-", dosage, "-avgNormFragmentCounts")]]
  
  print(sprintf("pseudocounts added to %d etoh peaks", sum(etohAvgNormFragmentCounts == 0)))
  etohAvgNormFragmentCounts[etohAvgNormFragmentCounts==0] <- 0.5
  print(sprintf("pseudocounts added to %d ra %s peaks", sum(raAvgNormFragmentCounts == 0), dosage))
  raAvgNormFragmentCounts[raAvgNormFragmentCounts==0] <- 0.5
  print(sprintf("pseudocounts added to %d tgfb %s peaks", sum(tgfbAvgNormFragmentCounts == 0), dosage))
  tgfbAvgNormFragmentCounts[tgfbAvgNormFragmentCounts==0] <- 0.5
  print(sprintf("pseudocounts added to %d both %s peaks", sum(bothAvgNormFragmentCounts == 0), dosage))
  bothAvgNormFragmentCounts[bothAvgNormFragmentCounts==0] <- 0.5
  
  #model: effect of RA = A, effect of TGFb = B
  A <- raAvgNormFragmentCounts   - etohAvgNormFragmentCounts
  B <- tgfbAvgNormFragmentCounts - etohAvgNormFragmentCounts
  addPred  <- etohAvgNormFragmentCounts + A + B
  multPred <- (raAvgNormFragmentCounts / etohAvgNormFragmentCounts) * (tgfbAvgNormFragmentCounts / etohAvgNormFragmentCounts) * etohAvgNormFragmentCounts  # this can be simplified but is more clear in this form
  integrationConstant <- (etohAvgNormFragmentCounts*(bothAvgNormFragmentCounts - addPred))/(A * B)
  dpOutputTib[[paste0("peak_RA-avgNormFragCountEffect-", dosage)]] <- A
  dpOutputTib[[paste0("peak_TGFb-avgNormFragCountEffect-", dosage)]] <- B
  dpOutputTib[[paste0("peak_addPred-", dosage)]]  <- addPred
  dpOutputTib[[paste0("peak_multPred-", dosage)]] <- multPred
  dpOutputTib[[paste0("peak_integrationConstant-", dosage)]] <- integrationConstant
  dpOutputTib[[paste0("peak_addMultPredFcDiff-", dosage)]] <- (multPred - addPred) / (etohAvgNormFragmentCounts)
}

##############################################################################################################################
# Make a pooled coefficient of variation (CV) estimate across low, med, and high doses of each condition
#   (this pooled estimate allows us to take advantage of similar variance profiles across different dosages)
#   ANFK = average normalized fragment count
addPseudoConfIntervals <- function(gextib, metadata, condname1, condname2, condname3, groupCondName) {
  # usually, condname1 will be low dose and condname 3 will be high dose, but they can also be low density or high density for the EtOH case
  nDosesTested <- length(dosages)
  anfkStatsLow  <- oneConditionEstMeanAndVarANFK(condname1, gextib, metadata) 
  anfkStatsMed  <- oneConditionEstMeanAndVarANFK(condname2, gextib, metadata) %>% dplyr::select(-peakID)
  anfkStatsHigh <- oneConditionEstMeanAndVarANFK(condname3, gextib, metadata) %>% dplyr::select(-peakID)
  anfkStatsAll  <- cbind(anfkStatsLow, anfkStatsMed, anfkStatsHigh)
  anfkStatsAllSharedCvEst <- anfkStatsAll %>% 
    dplyr::select(matches("_peak_cvEstANFK")) %>% 
    mutate('peak_cvEstPooled' = rowSums(dplyr::select(anfkStatsAll, matches("_peak_cvEstANFK"))) / nDosesTested)
  
  # sometimes a variance value is zero which ruins the CV estimate. in this case we take the zero-variance
  # estimates out of the CV estimate.
  invalidCVEst     <- !is.finite(anfkStatsAllSharedCvEst$peak_cvEstPooled)
  invalidCVEstInds <- which(invalidCVEst)
  for (ii in invalidCVEstInds) {
    stopifnot(is.nan(anfkStatsAllSharedCvEst[["peak_cvEstPooled"]][ii]))
    lowCvVal  <- anfkStatsAllSharedCvEst[[paste0(condname1, "_peak_cvEstANFK")]][ii]
    medCvVal  <- anfkStatsAllSharedCvEst[[paste0(condname2, "_peak_cvEstANFK")]][ii]
    highCvVal <- anfkStatsAllSharedCvEst[[paste0(condname3, "_peak_cvEstANFK")]][ii]
    cvValVec  <- c(lowCvVal, medCvVal, highCvVal)
    lowIsNaN  <- is.nan(anfkStatsAllSharedCvEst[[paste0(condname1, "_peak_cvEstANFK")]][ii])
    medIsNaN  <- is.nan(anfkStatsAllSharedCvEst[[paste0(condname2, "_peak_cvEstANFK")]][ii])
    highIsNaN <- is.nan(anfkStatsAllSharedCvEst[[paste0(condname3, "_peak_cvEstANFK")]][ii])
    isNaNvec  <- c(lowIsNaN, medIsNaN, highIsNaN)
    denominatorForAverage <- sum(!isNaNvec)
    numeratorForAverage   <- sum(cvValVec[!isNaNvec])
    newAvgVal <- numeratorForAverage / denominatorForAverage
    anfkStatsAllSharedCvEst[["peak_cvEstPooled"]][ii] <- newAvgVal
  }
  # if all CV estimates are indeed zero, then just leave it as zero
  if (!all(!is.nan(anfkStatsAllSharedCvEst[["peak_cvEstPooled"]]))) {
    invalidCVEst <- !is.finite(anfkStatsAllSharedCvEst$peak_cvEstPooled)
    print(sprintf("warning: %d genes have CV estimates of zero for %s", sum(invalidCVEst), groupCondName))
    invalidCVEstInds <- which(invalidCVEst)
    anfkStatsAllSharedCvEst[["peak_cvEstPooled"]][invalidCVEstInds] <- 0
  }
  stopifnot(all(!is.nan(anfkStatsAllSharedCvEst[["peak_cvEstPooled"]])))
  gextib[[paste0(groupCondName,"_peak_cvEstPooled")]] <- anfkStatsAllSharedCvEst[["peak_cvEstPooled"]]
  
  ## now use a normal distribution and the CV estimates to create pseudo-confidence intervals for the measurements of ANFK in both signals
  for(condname in c(condname1, condname2, condname3)) {
    mean_vec <- anfkStatsAll[[paste0(condname, "_peak_meanEstANFK")]]
    sd_vec   <- mean_vec * gextib[[paste0(groupCondName,"_peak_cvEstPooled")]]
    
    gextib[[paste0(condname, "_peak_pseudoConfintLowerBoundaryANFK")]] <- qnorm(tail.test.lower, 
                                                                          mean = mean_vec, 
                                                                          sd = sd_vec)
    gextib[[paste0(condname, "_peak_pseudoConfintUpperBoundaryANFK")]] <- qnorm(tail.test.upper, 
                                                                          mean = mean_vec, 
                                                                          sd = sd_vec)
  }
  return(gextib)
}
##############################################################################################################################

samplemetadata   <- read_tsv(here('sampleMetadata_SI2-SI4.txt'))
relevantmetadata <- filter(samplemetadata, as.integer(substr(sampleID,1,2)) %in% c(1:36))  # discard RNA-seq samples that failed library prep (17, 18) or extra samples from second run (46-50) 

dpOutputTib <- addPseudoConfIntervals(dpOutputTib, relevantmetadata, "TGFb-and-RA-low", "TGFb-and-RA-med", "TGFb-and-RA-high", "TGFb-and-RA")
dpOutputTib <- addPseudoConfIntervals(dpOutputTib, relevantmetadata, "TGFb-low", "TGFb-med", "TGFb-high", "TGFb")
dpOutputTib <- addPseudoConfIntervals(dpOutputTib, relevantmetadata, "RA-low", "RA-med", "RA-high", "RA")
dpOutputTib <- addPseudoConfIntervals(dpOutputTib, relevantmetadata, "EtOH-halfDensity", "EtOH-nlDensity", "EtOH-highDensity", "EtOH")


## add categorical signal integration value for each gene. requires a model of combined treatment variability.
tgfbcols <- paste0('TGFb-', dosages, "-avgNormFragmentCounts")
racols   <- paste0('RA-', dosages, "-avgNormFragmentCounts")
tgfbAvgEffectAllDoses <- rowMeans(dpOutputTib[, tgfbcols]) - dpOutputTib[["EtOH-nlDensity-avgNormFragmentCounts"]]
raAvgEffectAllDoses   <- rowMeans(dpOutputTib[, racols]) - dpOutputTib[["EtOH-nlDensity-avgNormFragmentCounts"]]
# initialize values of integration category field
for (dosage in dosages) {
  dpOutputTib[[paste0("peak_integrationCategory-", dosage,"-dose")]] <- NA
}
# this for loop iterates over the index row and requires that the tibbles be arranged in the same gene order
for(ii in 1:nrow(dpOutputTib)) {
  effectDir <- determineCombEffectDirection(tgfbAvgEffectAllDoses[ii], raAvgEffectAllDoses[ii]) #  <- function(a_effect, b_effect) {
  dpOutputTib[[paste0("peak_integrationEffectDir")]][ii] <- effectDir
  for(dosage in dosages) {
    thisDoseCategoricalIntegrationMode <- NA
    if (effectDir == "both up") {
      thisDoseCategoricalIntegrationMode <- boundfxBothUp(dpOutputTib[[paste0("TGFb-and-RA-", dosage,"_peak_pseudoConfintLowerBoundaryANFK")]][ii],
                                                          dpOutputTib[[paste0("TGFb-and-RA-", dosage,"_peak_pseudoConfintUpperBoundaryANFK")]][ii],
                                                          dpOutputTib[[paste0("peak_addPred-", dosage)]][ii],
                                                          dpOutputTib[[paste0("peak_multPred-", dosage)]][ii])
    } else if (effectDir == "opposing") {
      thisDoseCategoricalIntegrationMode <- boundfxOpposing(dpOutputTib[[paste0("TGFb-and-RA-", dosage,"_peak_pseudoConfintLowerBoundaryANFK")]][ii],
                                                            dpOutputTib[[paste0("TGFb-and-RA-", dosage,"_peak_pseudoConfintUpperBoundaryANFK")]][ii],
                                                            dpOutputTib[[paste0("peak_addPred-", dosage)]][ii],
                                                            dpOutputTib[[paste0("peak_multPred-", dosage)]][ii])
    } else if (effectDir == "both down") {
      thisDoseCategoricalIntegrationMode <- boundfxBothDown(dpOutputTib[[paste0("TGFb-and-RA-", dosage,"_peak_pseudoConfintLowerBoundaryANFK")]][ii],
                                                            dpOutputTib[[paste0("TGFb-and-RA-", dosage,"_peak_pseudoConfintUpperBoundaryANFK")]][ii],
                                                            dpOutputTib[[paste0("peak_addPred-", dosage)]][ii],
                                                            dpOutputTib[[paste0("peak_multPred-", dosage)]][ii])
    } else {
      thisDoseCategoricalIntegrationMode <- "uncategorized"
      effectDir <- "uncategorized"
    }
    dpOutputTib[[paste0("peak_integrationCategory-", dosage,"-dose")]][ii] <- thisDoseCategoricalIntegrationMode
  }
}


###

write_tsv(dpOutputTib, dpOutputTibLoc)
