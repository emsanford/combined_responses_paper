library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)

##########################################################################################
deseq.tib <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.tsv'))
samplemetadata   <- read_tsv(here('sampleMetadata_SI2-SI4.txt'))
relevantmetadata <- filter(samplemetadata, as.integer(substr(sampleID,1,2)) %in% c(1:16, 19:36, 51, 52))  # discard RNA-seq samples that failed library prep (17, 18) or extra samples from second run (46-50) 
minimum.foldchange.highcond <- 0.5
minimum.TPM.control <- 2
nDosesTested <- 3

## signal effect direction should be in {"both up", "both down", "opposing"}
signal.effect.directions <- c("both down", "opposing", "both up")
make.beeswarm.plots <- F

# parameters used for defining the range of estimated normal distribution for observed data.
#    additive/multiplicative predictions' locations with respect to the confidence interval 
#    of this distribution determine whether genes get classified as "additive", "multiplicative",
#    "subadditive", "supermultiplicative", "in-between", or "unable to tell"
tail.size <- 0.10
tail.test.lower <- 0 + tail.size
tail.test.upper <- 1 - tail.size
##########################################################################################


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


conf.int.tib <- tibble(ensg = deseq.tib$ensg)

for (condition in c('EtOH', 'TGFb', 'RA', 'TGFb-and-RA')) {
  
  if (condition == 'EtOH') {
    condLow  <- oneConditionEstMeanAndVarTPM(paste0(condition, '-halfDensity'), deseq.tib, relevantmetadata)
    condMed  <- oneConditionEstMeanAndVarTPM(paste0(condition, '-nlDensity'), deseq.tib, relevantmetadata) %>% dplyr::select(-ensg)
    condHigh <- oneConditionEstMeanAndVarTPM(paste0(condition, '-highDensity'), deseq.tib, relevantmetadata) %>% dplyr::select(-ensg)
  } else {
    condLow  <- oneConditionEstMeanAndVarTPM(paste0(condition, '-low'), deseq.tib, relevantmetadata)
    condMed  <- oneConditionEstMeanAndVarTPM(paste0(condition, '-med'), deseq.tib, relevantmetadata) %>% dplyr::select(-ensg)
    condHigh <- oneConditionEstMeanAndVarTPM(paste0(condition, '-high'), deseq.tib, relevantmetadata) %>% dplyr::select(-ensg)
  }

  condAll  <- cbind(condLow, condMed, condHigh)
  
  # with na.rm = TRUE, throw away zero-valued variance estimates instead of the whole CV estimate
  bothAllSharedCvEst <- condAll %>% 
    dplyr::select(matches('_cvEstTPM')) %>% 
    mutate('cvEstPooled' = rowMeans(dplyr::select(condAll, matches('_cvEstTPM')), na.rm = TRUE))
  
  # calculate upper and lower confidence intervals for all conditions, assuming shared CV estimates within conditions
  if (condition == 'EtOH') {
    for (density in c("halfDensity", "nlDensity", "highDensity")) {
      mean_vec <- condAll[[paste0(condition, "-", density, "_meanEstTPM")]]
      sd_vec   <- mean_vec * bothAllSharedCvEst$cvEstPooled
      conf.int.tib[[paste0(condition, "-", density, "_lowerBoundary")]] <- qnorm(tail.test.lower, 
                                                                                mean = mean_vec, 
                                                                                sd = sd_vec)
      conf.int.tib[[paste0(condition, "-", density, "_upperBoundary")]] <- qnorm(tail.test.upper, 
                                                                                mean = mean_vec, 
                                                                                sd = sd_vec)
    }
    #do something similar
  } else {
    for (dosage in c("low", "med", "high")) {
      mean_vec <- condAll[[paste0(condition, "-", dosage, "_meanEstTPM")]]
      sd_vec   <- mean_vec * bothAllSharedCvEst$cvEstPooled
      conf.int.tib[[paste0(condition, "-", dosage, "_lowerBoundary")]] <- qnorm(tail.test.lower, 
                                                                                mean = mean_vec, 
                                                                                sd = sd_vec)
      conf.int.tib[[paste0(condition, "-", dosage, "_upperBoundary")]] <- qnorm(tail.test.upper, 
                                                                                mean = mean_vec, 
                                                                                sd = sd_vec)
    }
  }
}

tpm.only.subset <- deseq.tib %>% dplyr::select(matches("^(TGFb|RA|TGFb-and-RA)-(low|med|high)_avgTPM$"))
tpm.only.subset[["EtOH-nlDensity_avgTPM"]] <- deseq.tib[["EtOH-nlDensity_avgTPM"]]

# add additive and multiplicative predictions to table
for (dosage in c("low", "med", "high")) {
  etohTPM <- tpm.only.subset[["EtOH-nlDensity_avgTPM"]]
  
  print(sprintf("pseudocounts added to %d genes to avoid zero-valued denominators in fold-change calculations", sum(etohTPM == 0)))
  etohTPM[etohTPM==0] <- 0.1
  
  raTPM   <- tpm.only.subset[[paste0("RA-", dosage, "_avgTPM")]]
  tgfbTPM <- tpm.only.subset[[paste0("TGFb-", dosage, "_avgTPM")]]
  # model: effect of RA = A, effect of TGFb = B
  A <- raTPM - etohTPM
  B <- tgfbTPM - etohTPM
  addPred  <- etohTPM + A + B
  multPred <- (raTPM / etohTPM) * (tgfbTPM / etohTPM) * etohTPM  # this expression can be simplified but is more clear in this form
  
  fcDiffMultVsAddPred <- (A*B)/(etohTPM * etohTPM)
  
  conf.int.tib[[paste0("addPred-", dosage, "-dose")]]  <- addPred
  conf.int.tib[[paste0("multPred-", dosage,  "-dose")]] <- multPred
  conf.int.tib[[paste0("fcDiffMultVsAddPred-", dosage,  "-dose")]] <- fcDiffMultVsAddPred
  
  # now do the same prediction backwards; what if we treated this as if we were removing drugs?
  # we then try to predict the etoh value from the drugged values
  bothTPM <- tpm.only.subset[[paste0("TGFb-and-RA-", dosage, "_avgTPM")]]
  A2 <- raTPM - bothTPM
  B2 <- tgfbTPM - bothTPM
  addPredBackwards <- bothTPM + A2 + B2
  multPredBackwards <- (raTPM / bothTPM) * (tgfbTPM / bothTPM) * bothTPM
  fcDiffMultVsAddPredBackwards <- (A2*B2)/(bothTPM * bothTPM)
  conf.int.tib[[paste0("addPredBackwards-", dosage, "-dose")]]  <- addPredBackwards
  conf.int.tib[[paste0("multPredBackwards-", dosage,  "-dose")]] <- multPredBackwards
  conf.int.tib[[paste0("fcDiffMultVsAddPredBackwards-", dosage,  "-dose")]] <- fcDiffMultVsAddPredBackwards
}

write_tsv(conf.int.tib, here("extractedData", "geneExprConfIntsAndAddMultPredictions.tsv"), col_names = T)


