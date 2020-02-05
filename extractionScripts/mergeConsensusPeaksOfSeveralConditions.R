# requires bedtools to be installed on the command line

mergeOrder <- c(here('extractedData', 'consensusPeakFiles', 'TGFb-and-RA-low.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'TGFb-and-RA-med.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'TGFb-and-RA-high.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'TGFb-low.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'TGFb-med.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'TGFb-high.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'RA-low.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'RA-med.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'RA-high.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'EtOH-nlDensity.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'EtOH-halfDensity.consensusSummits.bed'),  
                here('extractedData', 'consensusPeakFiles', 'EtOH-highDensity.consensusSummits.bed'))

mergeOrder <- map(mergeOrder, function (x) sprintf('"%s"', x)) #sprintf functions needed for command line parsing

outputMergedPeaksFile <- sprintf('"%s"', here('extractedData', 'mergedConsensusPeakFiles', 'allCondsMergedSummitWindows.bed'))

tempBedFile.forMerging  <- sprintf('"%s"', here('extractedData', 'mergedConsensusPeakFiles', 'tempMergedFile.bed'))
tempBedFile.forMerging2 <- sprintf('"%s"', here('extractedData', 'mergedConsensusPeakFiles', 'tempMergedFile2.bed'))
tempBedFile.forDiffs    <- sprintf('"%s"', here('extractedData', 'mergedConsensusPeakFiles', 'tempDiffPeaks.bed'))
tempBedFile.forSorting  <- sprintf('"%s"', here('extractedData', 'mergedConsensusPeakFiles', 'tempMergedSortedFile.bed'))
#step 1: merge overlapping segments in the first BED file
system(paste('bedtools', 'merge', '-c', '4,5', '-o', 'collapse,max', '-i', mergeOrder[1], '>', tempBedFile.forSorting))  
#next steps: in order, add new conditions to existing peaks, only if they don't already overlap a current peak
for (fileToMergeIn in mergeOrder[2:length(mergeOrder)]) {
  system(paste('bedtools', 'merge', '-c', '4,5', '-o', 'collapse,max', '-i', fileToMergeIn, '>', tempBedFile.forMerging2))
  system(paste('bedtools', 'intersect', '-v', '-sorted', '-a', tempBedFile.forMerging2, '-b', tempBedFile.forSorting, '>', tempBedFile.forDiffs))
  system(paste('cat', tempBedFile.forDiffs, tempBedFile.forSorting, ">", tempBedFile.forMerging))
  system(paste('sort', '-k1,1', '-k2,2n', tempBedFile.forMerging, '>', tempBedFile.forSorting))
}
system(paste('cp', tempBedFile.forSorting, outputMergedPeaksFile))
system(paste('rm', tempBedFile.forMerging, tempBedFile.forMerging2, tempBedFile.forDiffs, tempBedFile.forSorting))