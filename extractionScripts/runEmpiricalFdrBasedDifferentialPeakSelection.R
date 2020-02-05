library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)

############################## user-defined parameters go in this block of code ##################################
cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  absoluteMinPeakCoverage  <- 30
  lowerBoundFoldChange <- 1.5 
  fragmentCounts <- readRDS(here('extractedData', 'atacFragmentCountsAllCondsMergedDifferentialPeaks.rds'))  # ran
  outputFileDifferentialPeaks <- here('extractedData', 'differentialAtacPeaks.tsv')
  outputFileFDR_plot_prefix <- here('plots', 'differentialAtacPeaks_FDR_plot')
} else {
  absoluteMinPeakCoverage <- as.numeric(cmdargs[1])
  lowerBoundFoldChange    <- as.numeric(cmdargs[2])
  fragmentCounts <- readRDS(cmdargs[3])  
  outputFileDifferentialPeaks <- cmdargs[4]
  outputFileFDR_plot_prefix <- cmdargs[5]

}

sampleMetadata        <- read_tsv(here('sampleMetadata_SI2-SI4.txt'), col_names=TRUE)
tn5insertionPointsBed <- paste(sampleMetadata$ATAC_analysisDir, sampleMetadata$sampleID, '.Tn5_insertion_points.tagAlign.gz', sep="")
sampleMetadata        <- mutate(sampleMetadata, tn5insertionPointsBed=tn5insertionPointsBed)
filtMetadataOrdered   <- filter(sampleMetadata, between(as.integer(substr(sampleID,1,2)),1,36)) # discard metadata for RNA-seq technical replicates
conditionsToCount     <- unique(filtMetadataOrdered$condition)
nConditions           <- length(conditionsToCount)
nReplicates           <- 3

countMatrix           <- as.matrix(counts(fragmentCounts))
colnames(countMatrix) <- filtMetadataOrdered$`condition-replicate`
nSamples              <- length(colnames(countMatrix))
nPeaks                <- nrow(countMatrix)

sampleSizeFactors     <- nSamples * colSums(countMatrix) / sum(countMatrix)
normalizedCountMatrix <- t(t(countMatrix) / sampleSizeFactors)  # normalization is based on the number of reads in peaks

peakMeanCoverages        <- rowMeans(countMatrix)       
fdrGridLowerPeakCoverage <- max(absoluteMinPeakCoverage, quantile(peakMeanCoverages, .01))
fdrGridUpperPeakCoverage <- quantile(peakMeanCoverages, .99)
minCovgsToTest           <- 2^(seq(log2(fdrGridLowerPeakCoverage), log2(fdrGridUpperPeakCoverage), length.out=50)) # logarithmic steps for coverage
upperBoundFoldChange <- 10
foldChangeRangesToTest   <- 2^(seq(log2(lowerBoundFoldChange), log2(upperBoundFoldChange), length.out=50))
# minCovgsToTest         <- c(20, 100)
# foldChangeRangesToTest <- c(1.8, 3.6)

referenceControlName   <- c("EtOH-nlDensity")
additionalControlNames <- c("EtOH-halfDensity", "EtOH-highDensity")
experimentNames        <- sort(setdiff(conditionsToCount, c(referenceControlName, additionalControlNames)))
# the maximum fdr of a parameter set to be included in the "group" of acceptable fold-change/min normalized coverage parameter pairings
singleton.fdr.threshold <- 0.0025
pool.the.replicates     <- T
if (pool.the.replicates) {
  nReplicates <- 1
}
############################## note there is a nother manual code chunk at the bottom ###########################


# define function to get an FDR and list of differential peak indices from controls and experimental conditions
empiricalFdrByCountAndFoldChange <- function(normalizedCountMatrix, filtMetadataOrdered, primaryControlName, additionalControlNames, experimentCondNames, minFoldChange, minNormalizedCovg, pool_replicates = F) {

  
  allConditionNames <- c(primaryControlName, additionalControlNames, experimentCondNames)
  allRepsNumDiffPeaksControls      <- c()
  allRepsControlDiffPeakInds       <- c()
  allRepsNumDiffPeaksExperiments   <- c()
  allRepsDiffPeakIndsExperiments   <- c()
  allRepsExperimentDiffPeakIndList <- c(list(), list(), list())
  allRepsAdtlCntrlDiffPeakIndList  <- c(list(), list(), list())
  
  replicates <- unique(filtMetadataOrdered$replicate)
  
  if (pool_replicates) {
    nExperimentSamples <- length(experimentCondNames)
    nConds             <- length(allConditionNames)
    condAvgNormalizedCountMatrix <- matrix(NA, nrow(normalizedCountMatrix), nConds)
    pooledConditionNames <- paste0(allConditionNames, '-repsPooled')
    colnames(condAvgNormalizedCountMatrix) <- pooledConditionNames
    for (condName in allConditionNames) {
      pooledCondName <- paste0(condName, '-repsPooled')
      conditionAllReps <- filter(filtMetadataOrdered, condition == condName)$`condition-replicate`
      condAvgNormalizedCountMatrix[, pooledCondName] <- rowMeans(normalizedCountMatrix[, conditionAllReps])
    }
    pooledMetadata <- tibble(`condition-replicate` = pooledConditionNames, 
                             condition = allConditionNames, 
                             replicate = rep('repsPooled', nConds))
    diffpeak.result.list <- empiricalFdrByCountAndFoldChangeOneReplicate(condAvgNormalizedCountMatrix, 
                                                                         pooledMetadata, 
                                                                         primaryControlName, 
                                                                         additionalControlNames, 
                                                                         experimentCondNames, 
                                                                         minFoldChange, 
                                                                         minNormalizedCovg)
    allRepsNumDiffPeaksControls           <- diffpeak.result.list[[1]]
    allRepsControlDiffPeakInds            <- diffpeak.result.list[[2]]
    allRepsAdtlCntrlDiffPeakIndList[[1]]  <- diffpeak.result.list[[3]]
    allRepsNumDiffPeaksExperiments        <- diffpeak.result.list[[4]]
    allRepsDiffPeakIndsExperiments        <- diffpeak.result.list[[5]]
    allRepsExperimentDiffPeakIndList[[1]] <- diffpeak.result.list[[6]]
    
  } else {
    nExperimentSamples <- nrow(filter(filtMetadataOrdered, condition %in% experimentCondNames))
    replCounter <- 1
    for (repl in replicates) {
      # summary of one iteration: 
      #  1. filter metadata to single replicate
      #  2. call the one replicate function
      #  3. append results to the "allreps" variables
      filtReplMetadataOrdered <- filter(filtMetadataOrdered, replicate == repl)
      repl.diffpeak.result.list <- empiricalFdrByCountAndFoldChangeOneReplicate(normalizedCountMatrix, 
                                                                                filtReplMetadataOrdered, 
                                                                                primaryControlName, 
                                                                                additionalControlNames, 
                                                                                experimentCondNames, 
                                                                                minFoldChange, 
                                                                                minNormalizedCovg)
      ndiffpeaksControls             <- repl.diffpeak.result.list[[1]]
      controlDiffPeakIndsUnion       <- repl.diffpeak.result.list[[2]]
      indAddtControlDiffpeakIndsList <- repl.diffpeak.result.list[[3]]
      ndiffpeaksExperiments          <- repl.diffpeak.result.list[[4]]
      experimentDiffPeakIndsUnion    <- repl.diffpeak.result.list[[5]]
      indExperimentDiffpeakIndsList  <- repl.diffpeak.result.list[[6]]

      allRepsNumDiffPeaksControls      <- c(allRepsNumDiffPeaksControls, ndiffpeaksControls)
      allRepsControlDiffPeakInds       <- union(allRepsControlDiffPeakInds, controlDiffPeakIndsUnion)
      allRepsNumDiffPeaksExperiments   <- c(allRepsNumDiffPeaksExperiments, ndiffpeaksExperiments)
      allRepsDiffPeakIndsExperiments   <- union(allRepsDiffPeakIndsExperiments, experimentDiffPeakIndsUnion)
      allRepsExperimentDiffPeakIndList[[replCounter]] <- indExperimentDiffpeakIndsList
      allRepsAdtlCntrlDiffPeakIndList[[replCounter]]  <- indAddtControlDiffpeakIndsList
      
      replCounter <- replCounter + 1
    }
  }
  
  est.num.false.discoveries.in.experiments <- mean(allRepsNumDiffPeaksControls) * nExperimentSamples
  est.fdr <- (est.num.false.discoveries.in.experiments) / (sum(allRepsNumDiffPeaksExperiments))
  
  print(sprintf('Fold-change cutoff is %f, min nl covg is %f', minFoldChange, minNormalizedCovg))
  print(sprintf('Num differential peaks in controls (with duplicates): %d', sum(allRepsNumDiffPeaksControls)))
  print(sprintf('Num differential peaks in experiments (with duplicates): %d', sum(allRepsNumDiffPeaksExperiments)))
  print(sprintf('est FDR: %f', est.fdr))
  
  numDiffPeaksWithoutDuplicates <- length(allRepsDiffPeakIndsExperiments)
  return(list(est.fdr, numDiffPeaksWithoutDuplicates, allRepsControlDiffPeakInds, allRepsDiffPeakIndsExperiments, allRepsExperimentDiffPeakIndList, allRepsAdtlCntrlDiffPeakIndList))
}

empiricalFdrByCountAndFoldChangeOneReplicate <- function(normalizedCountMatrix, 
                                                         filtReplMetadataOrdered, 
                                                         primaryControlName, 
                                                         additionalControlNames, 
                                                         experimentCondNames, 
                                                         minFoldChange, 
                                                         minNormalizedCovg) {
  # Computes an empirical false discovery rate for one replicate (or one pooled replicate if pooling is done upstream)
  stopifnot(length(unique(filtReplMetadataOrdered$replicate)) == 1) # throw error if more than one replicate present; the sample metadata file must be filtered or modified upstream of this function
  ndiffpeaksControls <- c()
  controlDiffPeakIndsUnion <- c()
  indAddtControlDiffpeakIndsList <- list()
  counter <- 1
  for (addtlControlName in additionalControlNames) {
    diffPeakInds <- calcPairwiseDiffPeaksDataWranglingWrapper(primaryControlName, addtlControlName, minNormalizedCovg, minFoldChange, normalizedCountMatrix, filtReplMetadataOrdered)
    ndiffpeaksControls <- c(ndiffpeaksControls, length(diffPeakInds))
    controlDiffPeakIndsUnion <- union(controlDiffPeakIndsUnion, diffPeakInds)
    indAddtControlDiffpeakIndsList[[counter]] <- diffPeakInds
    counter <- counter + 1
  }
  
  ndiffpeaksExperiments <- c()
  experimentDiffPeakIndsUnion <- c()
  indExperimentDiffpeakIndsList <- list()
  counter <- 1
  for (experimentCondName in experimentCondNames) {
    diffPeakInds <- calcPairwiseDiffPeaksDataWranglingWrapper(primaryControlName, experimentCondName, minNormalizedCovg, minFoldChange, normalizedCountMatrix, filtReplMetadataOrdered)
    ndiffpeaksExperiments <- c(ndiffpeaksExperiments, length(diffPeakInds))
    experimentDiffPeakIndsUnion <- union(experimentDiffPeakIndsUnion, diffPeakInds)
    indExperimentDiffpeakIndsList[[counter]] <- diffPeakInds
    counter <- counter + 1
  }
  
  return(list(ndiffpeaksControls, controlDiffPeakIndsUnion, indAddtControlDiffpeakIndsList, ndiffpeaksExperiments, experimentDiffPeakIndsUnion, indExperimentDiffpeakIndsList))
}

calcPairwiseDiffPeaksDataWranglingWrapper <- function(controlName, sampleName, minNormalizedCovg, minFoldChange, normalizedCountMatrix, filtReplMetadataOrdered) {
  # use metadata to get the condition-replicate name (which matches the matrix name). then extract the normalized counts using column-name based indexing, then food to calcPairwiseDiffPeaks and return its result TK todo
  controlSampleName <- filter(filtReplMetadataOrdered, condition == controlName)$`condition-replicate`[1]
  otherSampleName   <- filter(filtReplMetadataOrdered, condition == sampleName)$`condition-replicate`[1]
  controlNormalizedCounts <- normalizedCountMatrix[,controlSampleName]
  otherNormalizedCounts   <- normalizedCountMatrix[,otherSampleName]
  return(calcPairwiseDiffPeaks(controlNormalizedCounts, otherNormalizedCounts, minNormalizedCovg, minFoldChange))
}

calcPairwiseDiffPeaks <- function(sample.1.counts.normalized, sample.2.counts.normalized, minNormalizedCovg, minFoldChange) {
  sample1.foldChanges <- sample.1.counts.normalized / sample.2.counts.normalized
  sample2.foldChanges <- sample.2.counts.normalized / sample.1.counts.normalized
  
  sample1.peaks.with.sufficient.depth <- sample.1.counts.normalized >= minNormalizedCovg
  sample2.peaks.with.sufficient.depth <- sample.2.counts.normalized >= minNormalizedCovg
  peaks.with.sufficient.depth <- sample1.peaks.with.sufficient.depth | sample2.peaks.with.sufficient.depth
  
  s1.enriched.candPeakIndices <- which((sample1.foldChanges > minFoldChange) & peaks.with.sufficient.depth)  ## todo: decide what to do with infinite values downstream...
  s2.enriched.candPeakIndices <- which((sample2.foldChanges > minFoldChange) & peaks.with.sufficient.depth)
  diffPeakIndices <- union(s1.enriched.candPeakIndices, s2.enriched.candPeakIndices)
  return(diffPeakIndices)
}

### End of function definitions, beginning of script 

# first, perform a grid search of FDRs over a range of foldChange and minimum normalized coverage thresholds
fdrGridSearchMatx <- matrix(NA, length(foldChangeRangesToTest), length(minCovgsToTest))
for (ii in 1:length(foldChangeRangesToTest)) {
  for (jj in 1:length(minCovgsToTest)) {
    
    diffPeakResultVec <- empiricalFdrByCountAndFoldChange(normalizedCountMatrix, 
                                                          filtMetadataOrdered, 
                                                          referenceControlName, 
                                                          additionalControlNames, 
                                                          experimentNames, 
                                                          foldChangeRangesToTest[ii], 
                                                          minCovgsToTest[jj], 
                                                          pool_replicates = pool.the.replicates)
    this.fdr <- diffPeakResultVec[[1]]
    this.nDiffPeak <- diffPeakResultVec[[2]]
    fdrGridSearchMatx[ii,jj] <- this.fdr
    
    if (ii == 1 & jj == 1) {
      tibbleForSettingFdrThreshold <- tibble(foldChanges=foldChangeRangesToTest[ii], minCovgs=minCovgsToTest[jj], fdr=this.fdr, nDiffPeak=this.nDiffPeak)
    } else {
      tibbleForSettingFdrThreshold <- add_row(tibbleForSettingFdrThreshold, foldChanges=foldChangeRangesToTest[ii], minCovgs=minCovgsToTest[jj], fdr=this.fdr, nDiffPeak=this.nDiffPeak)
    }
  }
}

## plot the FDR matrix as a color/heatmap grid
fdrGridPlot.linear <- ggplot(tibbleForSettingFdrThreshold, aes(foldChanges,minCovgs)) + geom_tile(aes(fill= log10(fdr) ))  # plot linear scale axes
fdrGridPlot.log <- ggplot(tibbleForSettingFdrThreshold, aes(log2(foldChanges),log2(minCovgs))) + geom_tile(aes(fill= log10(fdr) ), interpolate=TRUE)  # plot log scale axes
ggsave(paste0(outputFileFDR_plot_prefix, '_linscale.svg'), plot = fdrGridPlot.linear)
ggsave(paste0(outputFileFDR_plot_prefix, '_logscale.svg'), plot = fdrGridPlot.log)

# Now, based off a single FDR threshold, select multiple parameter sets at different foldChange / minCovg amounts that give this FDR, and combine into one final set of peaks
saved.fold.changes.ordered <- c()
saved.min.normalized.covg.ordered <- c()
corresponding.fdrs <- c()
for (ii in 1:length(foldChangeRangesToTest)) {
  fdr.diff.vec <- fdrGridSearchMatx[ii,] - singleton.fdr.threshold
  if (any(fdr.diff.vec < 0)){
    first.index.under.thresh <- which(fdr.diff.vec < 0)[1]
    fold.change.to.save <- foldChangeRangesToTest[ii]
    min.covg.to.save <- minCovgsToTest[first.index.under.thresh]
    saved.fold.changes.ordered <- c(saved.fold.changes.ordered, fold.change.to.save)
    saved.min.normalized.covg.ordered <- c(saved.min.normalized.covg.ordered, min.covg.to.save)
    corresponding.fdrs <- c(corresponding.fdrs, fdrGridSearchMatx[ii, first.index.under.thresh])
  }
}

# now collect all indices corresponding to selected parameter pairs under fdr threshold, calculate union and see what the combined fdr is 
# (the assumption here is that each condition contributes an equal number of false peaks, and that they can be estimated by the number of false peaks in the controls)
saved.diff.control.inds.union    <- c()
saved.diff.experiment.inds.union <- c()
ind.experiment.diffpeaks.list    <- rep(list(list()), nReplicates)
addtl.control.diffpeaks.list     <- rep(list(list()), nReplicates)
for (ii in 1:length(saved.fold.changes.ordered)) {
  diffPeakResultVec <- empiricalFdrByCountAndFoldChange(normalizedCountMatrix, 
                                                        filtMetadataOrdered, 
                                                        referenceControlName, 
                                                        additionalControlNames, 
                                                        experimentNames, 
                                                        saved.fold.changes.ordered[ii], 
                                                        saved.min.normalized.covg.ordered[ii], 
                                                        pool_replicates = pool.the.replicates)
  this.param.set.diff.control.inds                <- diffPeakResultVec[[3]]
  this.param.set.diff.experiment.inds             <- diffPeakResultVec[[4]]
  this.param.set.ind.experiment.diffpeak.ind.list <- diffPeakResultVec[[5]]
  this.param.set.ind.addtlcntrl.diffpeak.ind.list <- diffPeakResultVec[[6]]
  
  print(sprintf("Fold change = %.2f, ndiffpeak = %f", saved.fold.changes.ordered[ii], length(this.param.set.diff.experiment.inds)))
  saved.diff.control.inds.union <- union(saved.diff.control.inds.union, this.param.set.diff.control.inds)
  saved.diff.experiment.inds.union <- union(saved.diff.experiment.inds.union, this.param.set.diff.experiment.inds)
  condCounter <- 1
  for (sn in experimentNames) {
    for (replicateIndex in 1:nReplicates) {
      this.experiment.diffpeak.inds <- this.param.set.ind.experiment.diffpeak.ind.list[[replicateIndex]][[condCounter]]
      if (length(ind.experiment.diffpeaks.list[[replicateIndex]]) < condCounter) {
        ind.experiment.diffpeaks.list[[replicateIndex]][[condCounter]] <- this.experiment.diffpeak.inds
      } else {
        ind.experiment.diffpeaks.list[[replicateIndex]][[condCounter]] <- union(ind.experiment.diffpeaks.list[[replicateIndex]][[condCounter]], this.experiment.diffpeak.inds)
      }
    }
    condCounter <- condCounter + 1
  }
  condCounter <- 1
  for (sn in additionalControlNames) {
    for (replicateIndex in 1:nReplicates) {
      this.adtlControl.diffpeak.inds <- this.param.set.ind.addtlcntrl.diffpeak.ind.list[[replicateIndex]][[condCounter]]
      if (length(addtl.control.diffpeaks.list[[replicateIndex]]) < condCounter) {
        addtl.control.diffpeaks.list[[replicateIndex]][[condCounter]] <- this.adtlControl.diffpeak.inds
      } else {
        addtl.control.diffpeaks.list[[replicateIndex]][[condCounter]] <- union(addtl.control.diffpeaks.list[[replicateIndex]][[condCounter]], this.adtlControl.diffpeak.inds)
      }
    }
    condCounter <- condCounter + 1
  }
  
}


# 1. count the total number of differential peaks across the additional controls and the experimental conditions, 
#        (use the "union" of replicate calls at selected fold-change/coverage combinations)
# 2. use the peaks found in  additional controls to estimate the number of false positive peaks per sample,
# 3. multiply the number of false positive peaks per sample by the number of experimental samples to get
#        an estimate of the number of false positive peaks in the experimental conditions. 
# 4. Divide the number in *3) by the total number of called differential peaks in experiments to estimate the total
#        FDR for all selected peaks. 


if (pool.the.replicates) {
  sum.each.exp.replicate.union         <- Reduce("+", lapply(ind.experiment.diffpeaks.list, function(x) Reduce("+", lapply(x, function(y) length(y)))))
  sum.each.adtlcontrol.replicate.union <- Reduce("+", lapply(addtl.control.diffpeaks.list,  function(x) Reduce("+", lapply(x, function(y) length(y)))))
  nAddtlControls <- length(additionalControlNames)
  nExpermntConds <- length(experimentNames)
} else {
  sum.each.exp.replicate.union         <- Reduce("+", lapply(ind.experiment.diffpeaks.list, function(x) Reduce("+", lapply(x, function(y) length(y)))))
  sum.each.adtlcontrol.replicate.union <- Reduce("+", lapply(addtl.control.diffpeaks.list,  function(x) Reduce("+", lapply(x, function(y) length(y)))))
  nAddtlControls <- nReplicates * length(additionalControlNames)
  nExpermntConds <- nReplicates * length(experimentNames)
}

est.num.false.peaks.per.sample <- sum.each.adtlcontrol.replicate.union / nAddtlControls
final.combined.est.fdr <- (est.num.false.peaks.per.sample * nExpermntConds) / sum.each.exp.replicate.union
total.number.experiment.diff.peaks.union <- length(saved.diff.experiment.inds.union)
print(sprintf("Total number of differential peaks (union across all replicates): %d", total.number.experiment.diff.peaks.union))
print(sprintf("Final estimated combined FDR: %0.4f", final.combined.est.fdr))

#####################################################################################################################################
#####################################################################################################################################
# create tibble to write differential peaks to output file. indices are in BED format (start is 0-based, end is 1-based)
sortedDiffExperimentIndices <- sort(saved.diff.experiment.inds.union)
fragCountsGRanges <- granges(fragmentCounts)
chrom     <- as.vector(seqnames(fragCountsGRanges[sortedDiffExperimentIndices]))
startLocs <- start(fragCountsGRanges[sortedDiffExperimentIndices]) - 1
endLocs   <- end(fragCountsGRanges[sortedDiffExperimentIndices])
diffPeaksTibble <- tibble(chrom, startLocs, endLocs)

condCounter <- 1
for (condName in experimentNames) {
  replicateNames   <- c("-rep1", "-rep2", "-rep3")
  condRepNames     <- paste0(condName, replicateNames)
  replicateCounter <- 1
  # add replicate level fields
  for (replicate in replicateNames) {
    replicateName    <- paste0(condName, replicate)
    normFragCtName   <- paste0(replicateName, "-normFragmentCounts")
    foldChgName      <- paste0(replicateName, "-foldChange")
    diffPeakName     <- paste0(replicateName, "-isDiffPeak")
    replicateNormCvg <- normalizedCountMatrix[sortedDiffExperimentIndices, replicateName]
    refContrlNormCvg <- normalizedCountMatrix[sortedDiffExperimentIndices, paste0('EtOH-nlDensity', replicate)]
    foldchangeVector <- replicateNormCvg / refContrlNormCvg
    if (!pool.the.replicates) {
      replicateDiffPeakIndices <- ind.experiment.diffpeaks.list[[replicateCounter]][[condCounter]]
      stopifnot(all(replicateDiffPeakIndices %in% sortedDiffExperimentIndices))
      isDiffPeakThisRep <- sortedDiffExperimentIndices %in% replicateDiffPeakIndices
      diffPeaksTibble[[diffPeakName]]   <- isDiffPeakThisRep
    }
    
    diffPeaksTibble[[normFragCtName]] <- replicateNormCvg
    diffPeaksTibble[[foldChgName]]    <- foldchangeVector
    
    
    replicateCounter <- replicateCounter + 1
  }
  # now add condition-level fields; avg frag counts, avg fold change, isDiffPeak(Any/Two/All)Rep(s)
  
  condAverageNormFragCounts <- rowMeans(normalizedCountMatrix[sortedDiffExperimentIndices, c(condRepNames)])
  diffPeaksTibble[[paste0(condName, "-avgNormFragmentCounts")]] <- condAverageNormFragCounts
  
  etohAverageNormFragCounts <- rowMeans(normalizedCountMatrix[sortedDiffExperimentIndices, c(paste0('EtOH-nlDensity', replicateNames))])
  foldchangeVectorAvg       <- condAverageNormFragCounts / etohAverageNormFragCounts
  diffPeaksTibble[[paste0(condName, "-avgFoldchange")]] <- foldchangeVectorAvg
  
  if (pool.the.replicates) {
    diffPeaksTibble[[paste0(condName, "-isDiffPeak")]] <- sortedDiffExperimentIndices %in% ind.experiment.diffpeaks.list[[1]][[condCounter]]
  } else {
    # count how many times each different peak shows up across replicates using plyr count function (outputs data frame)
    counted <- plyr::count(c(ind.experiment.diffpeaks.list[[1]][[condCounter]], 
                             ind.experiment.diffpeaks.list[[2]][[condCounter]], 
                             ind.experiment.diffpeaks.list[[3]][[condCounter]]))
    diffPeakIndsAnyRep         <- filter(counted, freq >= 1)$x
    diffPeakIndsAtLeastTwoReps <- filter(counted, freq >= 2)$x
    diffPeakIndsAllReps        <- filter(counted, freq >= 3)$x
    diffPeaksTibble[[paste0(condName, "-isDiffPeakAnyRep")]]         <- sortedDiffExperimentIndices %in% diffPeakIndsAnyRep
    diffPeaksTibble[[paste0(condName, "-isDiffPeakAtLeastTwoReps")]] <- sortedDiffExperimentIndices %in% diffPeakIndsAtLeastTwoReps
    diffPeaksTibble[[paste0(condName, "-isDiffPeakAllReps")]]        <- sortedDiffExperimentIndices %in% diffPeakIndsAllReps
  }
  condCounter <- condCounter + 1
}


## copy-pasted form of above-block of code but for the additonal controls
condCounter <- 1
for (condName in additionalControlNames) {
  replicateNames   <- c("-rep1", "-rep2", "-rep3")
  condRepNames     <- paste0(condName, replicateNames)
  replicateCounter <- 1
  # add replicate level fields
  for (replicate in replicateNames) {
    replicateName    <- paste0(condName, replicate)
    normFragCtName   <- paste0(replicateName, "-normFragmentCounts")
    foldChgName      <- paste0(replicateName, "-foldChange")
    diffPeakName     <- paste0(replicateName, "-isDiffPeak")
    replicateNormCvg <- normalizedCountMatrix[sortedDiffExperimentIndices, replicateName]
    refContrlNormCvg <- normalizedCountMatrix[sortedDiffExperimentIndices, paste0('EtOH-nlDensity', replicate)]
    foldchangeVector <- replicateNormCvg / refContrlNormCvg
    if (! pool.the.replicates) {
      replicateDiffPeakIndices <- addtl.control.diffpeaks.list[[replicateCounter]][[condCounter]]
      isDiffPeakThisRep <- sortedDiffExperimentIndices %in% replicateDiffPeakIndices
      diffPeaksTibble[[diffPeakName]]   <- isDiffPeakThisRep
    }
    diffPeaksTibble[[normFragCtName]] <- replicateNormCvg
    diffPeaksTibble[[foldChgName]]    <- foldchangeVector
    
    replicateCounter <- replicateCounter + 1
  }
  # now add condition-level fields; avg frag counts, avg fold change, isDiffPeak(Any/Two/All)Rep(s)
  
  condAverageNormFragCounts <- rowMeans(normalizedCountMatrix[sortedDiffExperimentIndices, c(condRepNames)])
  diffPeaksTibble[[paste0(condName, "-avgNormFragmentCounts")]] <- condAverageNormFragCounts
  
  etohAverageNormFragCounts <- rowMeans(normalizedCountMatrix[sortedDiffExperimentIndices, c(paste0('EtOH-nlDensity', replicateNames))])
  foldchangeVectorAvg       <- condAverageNormFragCounts / etohAverageNormFragCounts
  diffPeaksTibble[[paste0(condName, "-avgFoldchange")]] <- foldchangeVectorAvg
  
  
  if (pool.the.replicates) {
    diffPeaksTibble[[paste0(condName, "-isDiffPeak")]] <- sortedDiffExperimentIndices %in% addtl.control.diffpeaks.list[[1]][[condCounter]]
  } else {
    # count how many times each different peak shows up across replicates using plyr count function (outputs data frame)
    counted <- plyr::count(c(addtl.control.diffpeaks.list[[1]][[condCounter]], 
                             addtl.control.diffpeaks.list[[2]][[condCounter]], 
                             addtl.control.diffpeaks.list[[3]][[condCounter]]))
    diffPeakIndsAnyRep         <- filter(counted, freq >= 1)$x
    diffPeakIndsAtLeastTwoReps <- filter(counted, freq >= 2)$x
    diffPeakIndsAllReps        <- filter(counted, freq >= 3)$x
    diffPeaksTibble[[paste0(condName, "-isDiffPeakAnyRep")]]         <- sortedDiffExperimentIndices %in% diffPeakIndsAnyRep
    diffPeaksTibble[[paste0(condName, "-isDiffPeakAtLeastTwoReps")]] <- sortedDiffExperimentIndices %in% diffPeakIndsAtLeastTwoReps
    diffPeaksTibble[[paste0(condName, "-isDiffPeakAllReps")]]        <- sortedDiffExperimentIndices %in% diffPeakIndsAllReps
  }
  condCounter <- condCounter + 1
}

# add the fragment counts for etoh control normal density
for (replicate in replicateNames) {
  replicateName    <- paste0('EtOH-nlDensity', replicate)
  normFragCtName   <- paste0(replicateName, "-normFragmentCounts")
  replicateNormCvg <- normalizedCountMatrix[sortedDiffExperimentIndices, replicateName]
  diffPeaksTibble[[normFragCtName]] <- replicateNormCvg
}
etohAverageNormFragCounts <- rowMeans(normalizedCountMatrix[sortedDiffExperimentIndices, c(paste0('EtOH-nlDensity', replicateNames))])
diffPeaksTibble[['EtOH-nlDensity-avgNormFragmentCounts']] <- etohAverageNormFragCounts

write_tsv(diffPeaksTibble, outputFileDifferentialPeaks, col_names = TRUE)
#####################################################################################################################################
#####################################################################################################################################
