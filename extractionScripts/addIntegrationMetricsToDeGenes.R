library("tidyverse")
library(here)

##########################################################################################


cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  deseqTib  <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.tsv'))
  samplemetadata   <- read_tsv(here('sampleMetadata_SI2-SI4.txt'))
  outputLoc <- here('extractedData', 'DeSeqOutputAllConds.annotated.tsv')
} else {
  deseqTib <- read_tsv(cmdargs[1])
  samplemetadata <- read_tsv(cmdargs[2])
  outputLoc <- cmdargs[3]
}

outputTib <- deseqTib  # we will progressively add to this tibble
relevantmetadata <- filter(samplemetadata, as.integer(substr(sampleID,1,2)) %in% c(1:16, 19:36, 51, 52))  # discard RNA-seq samples that failed library prep (17, 18) or extra samples from second run (46-50) 

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

oneConditionEstMeanAndVarTPM <- function(cond, deseq.tib, relevant.metadata) {
  filt.metadata <- filter(relevant.metadata, condition == cond)
  sample.names  <- filt.metadata$sampleID
  tpm.columns   <- deseq.tib[, paste0(sample.names, "_tpm")]
  mean_estimates <- rowMeans(tpm.columns)
  variance_estimates <- rowVars(tpm.columns)
  sum_squared_error <- rowSse(tpm.columns)
  cv_estimates <- rowCVs_unbiased(tpm.columns)
  tib.estimates <- tibble(ensg = deseq.tib$ensg)
  tib.estimates[[paste0(cond, '_meanEstTPM')]] <- mean_estimates
  tib.estimates[[paste0(cond, '_varEstTPM')]]  <- variance_estimates
  tib.estimates[[paste0(cond, '_sseTPM')]]     <- sum_squared_error
  tib.estimates[[paste0(cond, '_cvEstTPM')]]   <- cv_estimates
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
  etohTPM <- deseqTib[["EtOH-nlDensity_avgTPM"]]
  
  print(sprintf("pseudocounts added to %d genes", sum(etohTPM == 0)))
  etohTPM[etohTPM==0] <- 0.0001
  
  raTPM   <- deseqTib[[paste0("RA-", dosage, "_avgTPM")]]
  tgfbTPM <- deseqTib[[paste0("TGFb-", dosage, "_avgTPM")]]
  bothTPM <- deseqTib[[paste0("TGFb-and-RA-", dosage, "_avgTPM")]]
  #model: effect of RA = A, effect of TGFb = B
  A <- raTPM   - etohTPM
  B <- tgfbTPM - etohTPM
  addPred  <- etohTPM + A + B
  multPred <- (raTPM / etohTPM) * (tgfbTPM / etohTPM) * etohTPM  # this can be simplified but is more clear in this form
  integrationConstant <- (etohTPM*(bothTPM - addPred))/(A * B)
  outputTib[[paste0("RA-tpmEffect-", dosage)]] <- A
  outputTib[[paste0("TGFb-tpmEffect-", dosage)]] <- B
  outputTib[[paste0("addPred-", dosage)]]  <- addPred
  outputTib[[paste0("multPred-", dosage)]] <- multPred
  outputTib[[paste0("integrationConstant-", dosage)]] <- integrationConstant
  outputTib[[paste0("addMultPredFcDiff-", dosage)]] <- (multPred - addPred) / (etohTPM)
}

##############################################################################################################################
# Make a pooled coefficient of variation (CV) estimate across low, med, and high doses of each condition
#   (this pooled estimate allows us to take advantage of similar variance profiles across different dosages)
addPseudoConfIntervals <- function(gextib, metadata, condname1, condname2, condname3, groupCondName) {
  # usually, condname1 will be low dose and condname 3 will be high dose, but they can also be low density or high density for the EtOH case
  nDosesTested <- length(dosages)
  tpmStatsLow  <- oneConditionEstMeanAndVarTPM(condname1, gextib, metadata) 
  tpmStatsMed  <- oneConditionEstMeanAndVarTPM(condname2, gextib, metadata) %>% dplyr::select(-ensg)
  tpmStatsHigh <- oneConditionEstMeanAndVarTPM(condname3, gextib, metadata) %>% dplyr::select(-ensg)
  tpmStatsAll  <- cbind(tpmStatsLow, tpmStatsMed, tpmStatsHigh)
  tpmStatsAllSharedCvEst <- tpmStatsAll %>% 
    dplyr::select(matches('_cvEstTPM')) %>% 
    mutate('cvEstPooled' = rowSums(dplyr::select(tpmStatsAll, matches('_cvEstTPM'))) / nDosesTested)
  
  # sometimes a variance value is zero which ruins the CV estimate. in this case we take the zero-variance
  # estimates out of the CV estimate.
  invalidCVEst     <- !is.finite(tpmStatsAllSharedCvEst$cvEstPooled)
  invalidCVEstInds <- which(invalidCVEst)
  for (ii in invalidCVEstInds) {
    stopifnot(is.nan(tpmStatsAllSharedCvEst[["cvEstPooled"]][ii]))
    lowCvVal  <- tpmStatsAllSharedCvEst[[paste0(condname1, "_cvEstTPM")]][ii]
    medCvVal  <- tpmStatsAllSharedCvEst[[paste0(condname2, "_cvEstTPM")]][ii]
    highCvVal <- tpmStatsAllSharedCvEst[[paste0(condname3, "_cvEstTPM")]][ii]
    cvValVec  <- c(lowCvVal, medCvVal, highCvVal)
    lowIsNaN  <- is.nan(tpmStatsAllSharedCvEst[[paste0(condname1, "_cvEstTPM")]][ii])
    medIsNaN  <- is.nan(tpmStatsAllSharedCvEst[[paste0(condname2, "_cvEstTPM")]][ii])
    highIsNaN <- is.nan(tpmStatsAllSharedCvEst[[paste0(condname3, "_cvEstTPM")]][ii])
    isNaNvec  <- c(lowIsNaN, medIsNaN, highIsNaN)
    denominatorForAverage <- sum(!isNaNvec)
    numeratorForAverage   <- sum(cvValVec[!isNaNvec])
    newAvgVal <- numeratorForAverage / denominatorForAverage
    tpmStatsAllSharedCvEst[["cvEstPooled"]][ii] <- newAvgVal
  }
  # if all CV estimates are indeed zero, then just leave it as zero
  if (!all(!is.nan(tpmStatsAllSharedCvEst[["cvEstPooled"]]))) {
    invalidCVEst <- !is.finite(tpmStatsAllSharedCvEst$cvEstPooled)
    print(sprintf("warning: %d genes have CV estimates of zero for %s", sum(invalidCVEst), groupCondName))
    invalidCVEstInds <- which(invalidCVEst)
    tpmStatsAllSharedCvEst[["cvEstPooled"]][invalidCVEstInds] <- 0
  }
  stopifnot(all(!is.nan(tpmStatsAllSharedCvEst[["cvEstPooled"]])))
  gextib[[paste0(groupCondName,"_cvEstPooled")]] <- tpmStatsAllSharedCvEst[["cvEstPooled"]]
  
  ## now use a normal distribution and the CV estimates to create pseudo-confidence intervals for the measurements of TPM in both signals
  for(condname in c(condname1, condname2, condname3)) {
    mean_vec <- tpmStatsAll[[paste0(condname, "_meanEstTPM")]]
    sd_vec   <- mean_vec * gextib[[paste0(groupCondName,"_cvEstPooled")]]
    
    gextib[[paste0(condname, "_pseudoConfintLowerBoundaryTPM")]] <- qnorm(tail.test.lower, 
                                                                          mean = mean_vec, 
                                                                          sd = sd_vec)
    gextib[[paste0(condname, "_pseudoConfintUpperBoundaryTPM")]] <- qnorm(tail.test.upper, 
                                                                          mean = mean_vec, 
                                                                          sd = sd_vec)
  }
  return(gextib)
}
##############################################################################################################################
outputTib <- addPseudoConfIntervals(outputTib, relevantmetadata, "TGFb-and-RA-low", "TGFb-and-RA-med", "TGFb-and-RA-high", "TGFb-and-RA")
outputTib <- addPseudoConfIntervals(outputTib, relevantmetadata, "TGFb-low", "TGFb-med", "TGFb-high", "TGFb")
outputTib <- addPseudoConfIntervals(outputTib, relevantmetadata, "RA-low", "RA-med", "RA-high", "RA")
outputTib <- addPseudoConfIntervals(outputTib, relevantmetadata, "EtOH-halfDensity", "EtOH-nlDensity", "EtOH-highDensity", "EtOH")


## add categorical signal integration value for each gene. requires a model of combined treatment variability.
tgfbcols <- paste0('TGFb-', dosages, '_avgTPM')
racols   <- paste0('RA-', dosages, '_avgTPM')
tgfbAvgEffectAllDoses <- rowMeans(deseqTib[, tgfbcols]) - deseqTib[["EtOH-nlDensity_avgTPM"]]
raAvgEffectAllDoses   <- rowMeans(deseqTib[, racols]) - deseqTib[["EtOH-nlDensity_avgTPM"]]
# initialize values of integration category field
for (dosage in dosages) {
  outputTib[[paste0("integrationCategory-", dosage,"-dose")]] <- NA
}
# this for loop iterates over the index row and requires that the tibbles be arranged in the same gene order
for(ii in 1:nrow(outputTib)) {
  effectDir <- determineCombEffectDirection(tgfbAvgEffectAllDoses[ii], raAvgEffectAllDoses[ii]) #  <- function(a_effect, b_effect) {
  outputTib[[paste0("integrationEffectDir")]][ii] <- effectDir
  for(dosage in dosages) {
    thisDoseCategoricalIntegrationMode <- NA
    if (effectDir == "both up") {
      thisDoseCategoricalIntegrationMode <- boundfxBothUp(outputTib[[paste0("TGFb-and-RA-", dosage,"_pseudoConfintLowerBoundaryTPM")]][ii],
                                                          outputTib[[paste0("TGFb-and-RA-", dosage,"_pseudoConfintUpperBoundaryTPM")]][ii],
                                                          outputTib[[paste0("addPred-", dosage)]][ii],
                                                          outputTib[[paste0("multPred-", dosage)]][ii])
    } else if (effectDir == "opposing") {
      thisDoseCategoricalIntegrationMode <- boundfxOpposing(outputTib[[paste0("TGFb-and-RA-", dosage,"_pseudoConfintLowerBoundaryTPM")]][ii],
                                                            outputTib[[paste0("TGFb-and-RA-", dosage,"_pseudoConfintUpperBoundaryTPM")]][ii],
                                                            outputTib[[paste0("addPred-", dosage)]][ii],
                                                            outputTib[[paste0("multPred-", dosage)]][ii])
    } else if (effectDir == "both down") {
      thisDoseCategoricalIntegrationMode <- boundfxBothDown(outputTib[[paste0("TGFb-and-RA-", dosage,"_pseudoConfintLowerBoundaryTPM")]][ii],
                                                            outputTib[[paste0("TGFb-and-RA-", dosage,"_pseudoConfintUpperBoundaryTPM")]][ii],
                                                            outputTib[[paste0("addPred-", dosage)]][ii],
                                                            outputTib[[paste0("multPred-", dosage)]][ii])
    } else {
      thisDoseCategoricalIntegrationMode <- "uncategorized"
      effectDir <- "uncategorized"
    }
    outputTib[[paste0("integrationCategory-", dosage,"-dose")]][ii] <- thisDoseCategoricalIntegrationMode
  }
}

# add consensus signal integration categorical value; unambiguous values override ambiguous ones, 
#   and conflicting unambiguous values lead to a determination of "inconsistent"
outputTib[["consensusIntegrationCategory"]] <- NA
for(ii in 1:nrow(outputTib)) {
  consensusIntegrationCat <- NA
  integrationCategoriesThisGene <- as_vector(outputTib[ii, paste0("integrationCategory-", dosages,"-dose")])
  if (length(unique(integrationCategoriesThisGene)) == 1) {
    consensusIntegrationCat <- integrationCategoriesThisGene[1]
    # optional change to this logic: do not let one unambiguous cat override two ambiguous cats, i.e. if two ambiguous in the cat call it inconsistent or ambiguous
  } else if ((length(unique(integrationCategoriesThisGene)) == 2) & ("ambiguous" %in% integrationCategoriesThisGene) & !("uncategorized" %in% integrationCategoriesThisGene)) {
    consensusIntegrationCat <- integrationCategoriesThisGene[which(integrationCategoriesThisGene != "ambiguous")[1]]
  } else {
    consensusIntegrationCat <- "inconsistent"
  }
  outputTib[["consensusIntegrationCategory"]][ii] <- consensusIntegrationCat
}

# write output table to file
write_tsv(outputTib, outputLoc)
