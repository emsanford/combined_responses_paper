library(tidyverse)
library(here)


# this is janky, but you can switch between haigh and luow dosages using ctrl-F plus replacing
dosage_selection <- "low" 

subaddPeaksTib    <- read_tsv(here('extractedData', 'savedPeakSetsNearGeneSets', paste0('sub-additive_', dosage_selection ,'_dosage.tsv')))
addPeaksTib       <- read_tsv(here('extractedData', 'savedPeakSetsNearGeneSets', paste0('additive_', dosage_selection ,'_dosage.tsv')))
multPeaksTib      <- read_tsv(here('extractedData', 'savedPeakSetsNearGeneSets', paste0('multiplicative_', dosage_selection ,'_dosage.tsv')))
supermultPeaksTib <- read_tsv(here('extractedData', 'savedPeakSetsNearGeneSets', paste0('super-multiplicative_', dosage_selection ,'_dosage.tsv')))


motifColNames <- colnames(addPeaksTib)[grepl("_motifMatchScore", colnames(addPeaksTib))]
motifNames <- sapply(strsplit(motifColNames, "_"), function(x) x[1])
AP1.motifs <- c('FOSL2', 'FOSB', 'JUND', 'FOS', 'JUNB', 'FOSL1', 'JUN')
ELF.motifs <- c('ELF5')
FOXA.motifs <- c('FOXA2', 'FOXA1', 'FOXA3')
RAR.motifs  <- c('RARA', 'RARB', 'RARG', 'RXRA', 'RXRB', 'RXRG')
HOXA.motifs <- c('HOXA11', 'HOXA9', 'HOXA1', 'HOXA2', 'HOXA3', 'HOXA5', 'HOXA6', 'HOXA13', 'HOXA7', 'HOXA4')
SMAD.motifs <- c('SMAD4', 'SMAD3', 'SMAD1', 'SMAD2', 'SMAD5', 'SMAD9')
TGFB.motifs <- union(AP1.motifs, SMAD.motifs)
RA.motifs <- union(FOXA.motifs, union(HOXA.motifs, union(RAR.motifs, ELF.motifs)))

calcAdditivePeakMutualExclusivityScore <- function(etohCts, raCts, tgfbCts) {
  raEffectMagnitude   <- abs(raCts - etohCts)
  tgfbEffectMagnitude <- abs(tgfbCts - etohCts)
  lesserEffect <- min(raEffectMagnitude, tgfbEffectMagnitude)
  totalEffect  <- raEffectMagnitude + tgfbEffectMagnitude
  mutualExclusivityScore <- -1 * log2(lesserEffect / totalEffect)
  return(mutualExclusivityScore)
}

calcMultiplicativePeakMutualExclusivityScore <- function(etohCts, raCts, tgfbCts) {
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
  lesserEffect <- min(raEffectMagnitude, tgfbEffectMagnitude)
  totalEffect  <- raEffectMagnitude + tgfbEffectMagnitude
  mutualExclusivityScore <- -1 * log2(lesserEffect / totalEffect)
  return(mutualExclusivityScore)
}

addMutualExclusivityScoresToPeakTib <- function(peaktib) {
  addAvgMutualExclusivityScores <- c()
  multAvgMutualExclusivityScores <- c()
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
    multAvgMutualExclusivityScores <- c(multAvgMutualExclusivityScores, multAvgMutualExclusivityScore)
  }
  
  peaktib[["PeakMutualExclusivityScoreAdditive"]]       <- addAvgMutualExclusivityScores
  peaktib[["PeakMutualExclusivityScoreMultiplicative"]] <- multAvgMutualExclusivityScores
  return(peaktib)
}

calcAdditivePredDiff <- function(peaktib) {
  normcts.ra    <- peaktib$`TGFb-low-avgNormFragmentCounts`
  normcts.tgfb  <- peaktib$`RA-low-avgNormFragmentCounts`
  normcts.etoh  <- peaktib$`EtOH-nlDensity-avgNormFragmentCounts`
  addpred.cts   <- normcts.ra + normcts.tgfb - normcts.etoh
  addpred.fc    <- addpred.cts / normcts.etoh 
  
  both.fc <- peaktib$`TGFb-and-RA-low-avgFoldchange`
  fc.pred.diff <- both.fc - addpred.fc 
  
  return(fc.pred.diff)
}

clipEdges <- function(vec, minval, maxval) {
  vec[vec > maxval] <- maxval
  vec[vec < minval] <- minval
  return(vec)
}

motifSets <- list(AP1.motifs, ELF.motifs, FOXA.motifs, HOXA.motifs, SMAD.motifs, TGFB.motifs, RA.motifs, RAR.motifs)
motifSetNames <- c("AP1", "ELF", "FOXA", "HOXA", "SMAD", "RAR", "TGFB", "RA")
counter <- 1
for (motifSet in motifSets) {
  motifSetColNames <- paste0(motifSet, '_motifMatchScore')
  motifSetName <- motifSetNames[counter]

  subaddPeaksTib[[paste0(motifSetName, '_motifMatchInSet')]] <- rowSums(subaddPeaksTib[, motifSetColNames]) > 0   
  addPeaksTib[[paste0(motifSetName, '_motifMatchInSet')]] <- rowSums(addPeaksTib[, motifSetColNames]) > 0  
  multPeaksTib[[paste0(motifSetName, '_motifMatchInSet')]] <- rowSums(multPeaksTib[, motifSetColNames]) > 0
  supermultPeaksTib[[paste0(motifSetName, '_motifMatchInSet')]] <- rowSums(supermultPeaksTib[, motifSetColNames]) > 0 
  
  counter <- counter + 1
}

## TODO 
getGeneCentricMotifTibble <- function(somePeaksTib) {
  genes <- unique(somePeaksTib$ensg)
  genecentricOutputTibbleRows <- list()
  counter <- 1
  for (gene_id in genes) {
    motifSetMatches <- c()
    for (motifSetName in motifSetNames) {
      peaks.thisgene  <- filter(somePeaksTib, ensg==gene_id)
      motifSetMatches <- c(motifSetMatches, any(peaks.thisgene[[paste0(motifSetName, "_motifMatchInSet")]]))
    }
    genecentricOutputTibbleRows[[counter]] <- c(motifSetMatches)
    counter <- counter + 1
  }
  tempmatxthing <- NULL
  for (i in 1:length(genecentricOutputTibbleRows)) {
    tempmatxthing <- rbind(tempmatxthing, genecentricOutputTibbleRows[[i]])
  }
  tempmatxthing <- as.matrix(tempmatxthing)
  genecentricOutputTibble <- as_tibble(tempmatxthing)
  colnames(genecentricOutputTibble) <- motifSetNames
  genecentricOutputTibble[["ensg"]] <- genes
  return(genecentricOutputTibble)
}

subaddPeaksTibGeneCentric <- getGeneCentricMotifTibble(subaddPeaksTib)
addPeaksTibGeneCentric <- getGeneCentricMotifTibble(addPeaksTib)
multPeaksTibGeneCentric <- getGeneCentricMotifTibble(multPeaksTib)
supermultPeaksTibGeneCentric <- getGeneCentricMotifTibble(supermultPeaksTib)

colSums(subaddPeaksTibGeneCentric[,1:8]) / nrow(subaddPeaksTibGeneCentric)
colSums(addPeaksTibGeneCentric[,1:8]) / nrow(addPeaksTibGeneCentric)
colSums(multPeaksTibGeneCentric[,1:8]) / nrow(multPeaksTibGeneCentric)
colSums(supermultPeaksTibGeneCentric[,1:8]) / nrow(supermultPeaksTibGeneCentric)

c1 <- cor(subaddPeaksTibGeneCentric[,1:6])
c2 <- cor(addPeaksTibGeneCentric[,1:6])
c3 <- cor(multPeaksTibGeneCentric[,1:6])
c4 <- cor(supermultPeaksTibGeneCentric[,1:6])
library("matlab")
# c1["TGFB", ] <- 0
# c1[,"TGFB"] <- 0
# c1["TGFB", "TGFB"] <- 0
c1[6,6] <- -1
c2[6,6] <- -1
c3[6,6] <- -1
c4[6,6] <- -1
imagesc(c1, col = gray(seq(0,1, by = 0.01)), main='subadd, n=11')
imagesc(c2, col = gray(seq(0,1, by = 0.01)), main='add, n=34')
imagesc(c3, col = gray(seq(0,1, by = 0.01)), main='mult, n=45')
imagesc(c4, col = gray(seq(0,1, by = 0.01)), main='supermult, n=45')

# print frac TGFb matches across all peaks
# print frac RA matches
# print frac RA and TGFb matches (okay to be in different peaks)



library(gridExtra)
# avgFoldchange refers to peak, _log2fc refers to gene
p1 <- ggplot(subaddPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point() +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
  theme_minimal(base_size = 10) + theme(legend.position = "none") + ggtitle("peaks near sub-additive genes") +
  xlim(-2, 4) + ylim(-2, 4) +
  facet_wrap(TGFB_motifMatchInSet ~ RA_motifMatchInSet, nrow = 1)

p2 <- ggplot(addPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point() +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
  theme_minimal(base_size = 10) + theme(legend.position = "none") + ggtitle("peaks near additive genes") +
  xlim(-2, 4) + ylim(-2, 4) +
  facet_wrap(TGFB_motifMatchInSet ~ RA_motifMatchInSet, nrow = 1)

p3 <- ggplot(multPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point() +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
  theme_minimal(base_size = 10) + theme(legend.position = "none") + ggtitle("peaks near multiplicative genes") +
  xlim(-2, 4) + ylim(-2, 4) +
  facet_wrap(TGFB_motifMatchInSet ~ RA_motifMatchInSet, nrow = 1)

p4 <- ggplot(supermultPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point() +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
  theme_minimal(base_size = 10) + theme(legend.position = "none") + ggtitle("peaks near supermultiplicative genes") +
  xlim(-2, 4) + ylim(-2, 4) +
  facet_wrap(TGFB_motifMatchInSet ~ RA_motifMatchInSet, nrow = 1)

grid.arrange(p1, p2, p3, p4, nrow = 4, ncol = 1)
# 
# p1 <- ggplot(subaddPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
#   geom_point(mapping = aes(color = RA_motifMatchInSet + TGFB_motifMatchInSet)) +
#   xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
#   theme_minimal(base_size = 19) + theme(legend.position = "none") + ggtitle("peaks near sub-additive genes") + scale_fill_brewer(palette = "Reds") +
#   xlim(-2, 4) + ylim(-2, 4) + scale_colour_gradient2(low = "gray", mid = "pink", high = "red", midpoint = 1)
# 
# p2 <- ggplot(addPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
#   geom_point(mapping = aes(color = RA_motifMatchInSet + TGFB_motifMatchInSet)) +
#   xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
#   theme_minimal(base_size = 19) + theme(legend.position = "none") + ggtitle("peaks near additive genes") + scale_fill_brewer(palette = "Reds") +
#   xlim(-2, 4) + ylim(-2, 4) + scale_colour_gradient2(low = "gray", mid = "pink", high = "red", midpoint = 1)
# 
# p3 <- ggplot(multPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
#   geom_point(mapping = aes(color = RA_motifMatchInSet + TGFB_motifMatchInSet)) +
#   xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
#   theme_minimal(base_size = 19) + theme(legend.position = "none") + ggtitle("peaks near multiplicative genes") + scale_fill_brewer(palette = "Reds") +
#   xlim(-2, 4) + ylim(-2, 4) + scale_colour_gradient2(low = "gray", mid = "pink", high = "red", midpoint = 1)
# 
# p4 <- ggplot(supermultPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
#   geom_point(mapping = aes(color = RA_motifMatchInSet + TGFB_motifMatchInSet)) +
#   xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in retinoic acid") + 
#   theme_minimal(base_size = 19) + theme(legend.position = "none") + ggtitle("peaks near super-mult genes") + scale_fill_brewer(palette = "Reds") +
#   xlim(-2, 4) + ylim(-2, 4) + scale_colour_gradient2(low = "gray", mid = "pink", high = "red", midpoint = 1)


clipmin <- -5.0
clipmax <- 5.0
basefontsize <- 19
p5 <- ggplot(subaddPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = clipEdges(calcAdditivePredDiff(subaddPeaksTib), clipmin, clipmax))) +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in RA") + 
  theme_minimal(base_size = basefontsize) + theme(legend.position = "none") +
  labs(color = "peak diff pred", title = "sub-additive") + xlim(-2, 4) + ylim(-2, 4)  + scale_color_gradientn(colours = rainbow(5), limits = c(clipmin, clipmax))

p6 <- ggplot(addPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = clipEdges(calcAdditivePredDiff(addPeaksTib), clipmin, clipmax))) +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in RA") + 
  theme_minimal(base_size = basefontsize) + theme(legend.position = "none") +
  labs(color = "peak diff pred", title = "additive") + xlim(-2, 4) + ylim(-2, 4)  + scale_color_gradientn(colours = rainbow(5), limits = c(clipmin, clipmax))

p7 <- ggplot(multPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = clipEdges(calcAdditivePredDiff(multPeaksTib), clipmin, clipmax))) +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in RA") + 
  theme_minimal(base_size = basefontsize) + theme(legend.position = "none") +
  labs(color = "peak diff pred", title = "multiplicative") + xlim(-2, 4) + ylim(-2, 4)  + scale_color_gradientn(colours = rainbow(5), limits = c(clipmin, clipmax))

p8 <- ggplot(supermultPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = clipEdges(calcAdditivePredDiff(supermultPeaksTib), clipmin, clipmax))) +
  xlab("log(fold-change) in TGF-beta") + ylab("log(fold-change) in RA") + 
  theme_minimal(base_size = basefontsize) + theme(legend.position = "right") +
  labs(color = "FC-diff from\nadditive pred", title = "super-multiplicative") + xlim(-2, 4) + ylim(-2, 4)  + scale_color_gradientn(colours = rainbow(5), limits = c(clipmin, clipmax))


grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
grid.arrange(p5, p6, p7, p8, nrow = 2, ncol = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow=2, ncol=4)

###### which of these peaks are "pioneer opening peaks"?
normcvgthresh <- 4.0
p9  <- ggplot(subaddPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = subaddPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts` <= normcvgthresh)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "sub-additive") + xlim(-2, 4) + ylim(-2, 4)  

p10 <- ggplot(addPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = addPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts` <= normcvgthresh)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "additive") + xlim(-2, 4) + ylim(-2, 4)  

p11 <- ggplot(multPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = multPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts` <= normcvgthresh)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "multiplicative") + xlim(-2, 4) + ylim(-2, 4)  

p12 <-ggplot(supermultPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(color = supermultPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts` <= normcvgthresh)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "super-multiplicative") + xlim(-2, 4) + ylim(-2, 4)  

grid.arrange(p9, p10, p11, p12, nrow = 2, ncol = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow=3, ncol=4)

####### let's visualize fold change in the both condition with the size of the dots
p13  <- ggplot(subaddPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(size = subaddPeaksTib$`TGFb-and-RA-low-avgFoldchange`)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "sub-additive") + xlim(-2, 4) + ylim(-2, 4)  

p14 <- ggplot(addPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(size = addPeaksTib$`TGFb-and-RA-low-avgFoldchange`)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "additive") + xlim(-2, 4) + ylim(-2, 4)  

p15 <- ggplot(multPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(size = multPeaksTib$`TGFb-and-RA-low-avgFoldchange`)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "multiplicative") + xlim(-2, 4) + ylim(-2, 4)  

p16 <-ggplot(supermultPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(size = supermultPeaksTib$`TGFb-and-RA-low-avgFoldchange`)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "super-multiplicative") + xlim(-2, 4) + ylim(-2, 4)  

grid.arrange(p13, p14, p15, p16)


####### how dose-responsive are these peaks?
p16  <- ggplot(subaddPeaksTib, aes(x = log2(`TGFb-low-avgFoldchange`), y = log2(`RA-low-avgFoldchange`))) +
  geom_point(mapping = aes(size = subaddPeaksTib$`TGFb-and-RA-low-avgFoldchange`)) +
  labs(color = sprintf("EtOH nlCounts <= %0.1f", normcvgthresh), title = "sub-additive") + xlim(-2, 4) + ylim(-2, 4)  


#######
# print out statistics for frequency of "pioneer peaks" for each category
tibnames <- c("subaddPeaks", "addPeaks", "multPeaks", "supermultPeaks")
counter <- 1
for (tib in list(subaddPeaksTib, addPeaksTib, multPeaksTib, supermultPeaksTib)) {
  npeaks <- length(tib$ensg)
  pioneerPeakIndices <- tib$`EtOH-nlDensity-avgNormFragmentCounts` <= normcvgthresh
  npeaks_pioneer <- sum(pioneerPeakIndices)
  ngenes_uniq          <- length(unique(tib$ensg))
  ngenes_wPioneerPeaks <- length(unique(tib[pioneerPeakIndices,]$ensg))
  print(sprintf("%s: %d pioneer peaks, %d total peaks, %0.3f fraction of total peaks found near %d genes of %d total genes", 
                tibnames[counter], npeaks_pioneer, npeaks, npeaks_pioneer/npeaks, ngenes_wPioneerPeaks, ngenes_uniq))
  counter <- counter + 1
}
#######

#######
# print out statistics for frequency of motif matching for each category. Are we seeing colocalization of diff. motifs in the same peak?
tibnames <- c("subadd genes", "add genes", "mult genes", "supermult genes")
tibForSummary <- NULL
counter <- 1
for (tib in list(subaddPeaksTib, addPeaksTib, multPeaksTib, supermultPeaksTib)) {
  npeaks <- length(tib$ensg)
  npeaks_0_matches <- sum((!tib$RA_motifMatchInSet) & (!tib$TGFB_motifMatchInSet))
  npeaks_1_match   <- sum(xor(tib$RA_motifMatchInSet, tib$TGFB_motifMatchInSet))
  npeaks_2_matches <- sum((tib$RA_motifMatchInSet) & (tib$TGFB_motifMatchInSet))
  print(sprintf("%s: %0.3f zero matches, %0.3f one match, %0.3f two matches", 
                tibnames[counter], npeaks_0_matches/npeaks, npeaks_1_match/npeaks, npeaks_2_matches/npeaks))
  tibForSummary <- bind_rows(tibForSummary, tibble(integrMode = tibnames[counter], 
                                                   fracZeroMotifMatches = npeaks_0_matches/npeaks,
                                                   fracOneMotifMatch = npeaks_1_match/npeaks,
                                                   fracTwoMotifMatches = npeaks_2_matches/npeaks))
  counter <- counter + 1
}

# Kinetic model predicts multiplicative behavior doesn't require colocalization.
# are we seeing multiple classes of TFs activated by diff. signals for super-mult genes?


# make summary plot for number of motif matches in each peak near gene sets of each category
tibForSummary <- gather(tibForSummary, key = "PeakFracCat", value = "FracPeaks", c("fracZeroMotifMatches", "fracOneMotifMatch", "fracTwoMotifMatches"))
tibForSummary$PeakFracCat <- factor(tibForSummary$PeakFracCat, levels = c("fracZeroMotifMatches", "fracOneMotifMatch", "fracTwoMotifMatches"))
tibForSummary$integrMode <- factor(tibForSummary$integrMode, levels = c("subadd genes", "add genes", "mult genes", "supermult genes"))
pMotifFracs <- ggplot(tibForSummary, aes(x = integrMode, y = FracPeaks, fill = PeakFracCat)) + 
  geom_bar(stat="identity", position = "dodge") +
  xlab("mode of integration") + ylab("fraction of peaks within\n100kb of genes") + 
  theme_minimal(base_size = 24) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  xlab("") +
  scale_fill_brewer(palette = "Reds") + theme(legend.position = "bottom", legend.direction="vertical")
pMotifFracs

# test for correlation between integration mode and whether peaks are super-additive
counter <- 1
thresholdPredDiff <- c(1.0, 1.5, 2.0)
genesWithSuperAddPeakTib <- tibble(geneSet = character(), 
                                   geneFracWithSuperAddPeak = double(), 
                                   numGenesWithSuperAddPeak = integer(), 
                                   num_genes = integer(), 
                                   threshold = double())

origGeneSets <- here('extractedData', 'savedGeneSets', c('sub-additive_low_dosage', 
                                                         'additive_low_dosage', 
                                                         'multiplicative_low_dosage', 
                                                         'super-multiplicative_low_dosage'))
for (tib in list(subaddPeaksTib, addPeaksTib, multPeaksTib, supermultPeaksTib)) {
  ngenestot <- length(read_tsv(origGeneSets[counter], col_names = T)$ensg);
  
  ngeneswithsuperaddpeaks1 <- tib[calcAdditivePredDiff(tib) > thresholdPredDiff[1],] %>% 
    filter(`TGFb-and-RA-low-avgFoldchange` > 0) %>%
    select(ensg) %>% 
    unique() %>% 
    pull(ensg) %>% 
    length()
  ngeneswithsuperaddpeaks2 <- tib[calcAdditivePredDiff(tib) > thresholdPredDiff[2],] %>% 
    filter(`TGFb-and-RA-low-avgFoldchange` > 0) %>%
    select(ensg) %>% 
    unique() %>% 
    pull(ensg) %>% 
    length()
  ngeneswithsuperaddpeaks3 <- tib[calcAdditivePredDiff(tib) > thresholdPredDiff[3],] %>% 
    filter(`TGFb-and-RA-low-avgFoldchange` > 0) %>%
    select(ensg) %>% 
    unique() %>% 
    pull(ensg) %>% 
    length()
  
  print(sprintf("%s: %d genes with superaddpeaks, %d genes total, %0.2f percent", 
                tibnames[counter], ngeneswithsuperaddpeaks, ngenestot, ngeneswithsuperaddpeaks/ngenestot))
  
  genesWithSuperAddPeakTib[counter * 3 - 2,] <- list(tibnames[counter], ngeneswithsuperaddpeaks1/ngenestot, ngeneswithsuperaddpeaks1, ngenestot, thresholdPredDiff[1])
  genesWithSuperAddPeakTib[counter * 3 - 1,] <- list(tibnames[counter], ngeneswithsuperaddpeaks2/ngenestot, ngeneswithsuperaddpeaks2, ngenestot, thresholdPredDiff[2])
  genesWithSuperAddPeakTib[counter * 3,] <- list(tibnames[counter], ngeneswithsuperaddpeaks3/ngenestot, ngeneswithsuperaddpeaks3, ngenestot, thresholdPredDiff[3])
  
  counter <- counter + 1
}
genesWithSuperAddPeakTib$threshold <- factor(genesWithSuperAddPeakTib$threshold)
# plot bar chart showing "dose-dependent" increase of frequency of genes having superadditive peaks nearby
ggplot(genesWithSuperAddPeakTib, aes(x=factor(geneSet, levels = c("subadd genes", "add genes", "mult genes", "supermult genes")), 
                                     y=geneFracWithSuperAddPeak,
                                     fill=threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("fraction of genes with at least one\nsuperadditive peak nearby") + 
  xlab("") + 
  theme_minimal(base_size = 20) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_brewer(palette = "Blues") + 
  labs(fill = '"superadditive"\nthreshold')

# could do a similar analysis for just raw fold change

# test if number of differential elements nearby also predicts integration mode
avgNearbyPeaksUpTib <- tibble(geneSet = character(), 
                              avgNearbyPeaksTGFb = double(), 
                              avgNearbyPeaksRA = double(), 
                              avgNearbyPeaksBoth = double(), 
                              num_genes = integer())
counter <- 1
peak.up.FC.threshold <- 1
for (tib in list(subaddPeaksTib, addPeaksTib, multPeaksTib, supermultPeaksTib)) {
  genes.only.tib <- read_tsv(origGeneSets[counter], col_names = T)
  ngenestot <- length(genes.only.tib$ensg);
  
  nearbyPeaksUpTGFb <- tib %>%
    select(ensg, numNearbyTGFbDiffAccessibleElements, `TGFb-low-avgFoldchange`, `TGFb-low_isDiffPeakAtLeastTwoReps`) %>%
    filter(`TGFb-low-avgFoldchange` > peak.up.FC.threshold) %>%
    group_by(ensg) %>%
    mutate(nPeaksUpNear = n()) %>%
    select(ensg, nPeaksUpNear) %>%
    right_join(select(genes.only.tib, ensg), by = "ensg") %>%
    unique() 
  nearbyPeaksUpTGFb[is.na(nearbyPeaksUpTGFb)]  <- 0
  avgNumNearbyPeaksTGFb <- mean(pull(nearbyPeaksUpTGFb, nPeaksUpNear))
  
  nearbyPeaksUpRA <- tib %>%
    select(ensg, numNearbyRADiffAccessibleElements, `RA-low-avgFoldchange`, `RA-low_isDiffPeakAtLeastTwoReps`) %>%
    filter(`RA-low-avgFoldchange` > peak.up.FC.threshold) %>%
    group_by(ensg) %>%
    mutate(nPeaksUpNear = n()) %>%
    select(ensg, nPeaksUpNear) %>%
    right_join(select(genes.only.tib, ensg), by = "ensg") %>%
    unique() 
  nearbyPeaksUpRA[is.na(nearbyPeaksUpRA)]  <- 0
  avgNumNearbyPeaksRA <- mean(pull(nearbyPeaksUpRA, nPeaksUpNear))
  
  nearbyPeaksUpBoth <- tib %>%
    select(ensg, numNearbyBothDiffAccessibleElements, `TGFb-and-RA-low-avgFoldchange`, `TGFb-and-RA-low_isDiffPeakAtLeastTwoReps`) %>%
    filter(`TGFb-and-RA-low-avgFoldchange` > peak.up.FC.threshold) %>%
    group_by(ensg) %>%
    mutate(nPeaksUpNear = n()) %>%
    select(ensg, nPeaksUpNear) %>%
    right_join(select(genes.only.tib, ensg), by = "ensg") %>%
    unique() 
  nearbyPeaksUpBoth[is.na(nearbyPeaksUpBoth)]  <- 0
  avgNumNearbyPeaksBoth <- mean(pull(nearbyPeaksUpBoth, nPeaksUpNear))

  print(sprintf("%s: %0.2f avg peaks TGFb, %0.2f avg peaks RA, %0.2f avg peaks Both, %d genes total", 
                tibnames[counter], avgNumNearbyPeaksTGFb, avgNumNearbyPeaksRA, avgNumNearbyPeaksBoth, ngenestot))
  # print(nearbyPeaksUpBoth$nPeaksUpNear)
  avgNearbyPeaksUpTib[counter,] <- list(tibnames[counter], avgNumNearbyPeaksTGFb, avgNumNearbyPeaksRA, avgNumNearbyPeaksBoth, ngenestot)
  
  # avgNearbyPeaksUpTib <- rbind(avgNearbyPeaksUpTib, c(tibnames[counter], avgNumNearbyPeaksTGFb, avgNumNearbyPeaksRA, avgNumNearbyPeaksBoth, ngenestot))
  
  counter <- counter + 1
}
avgNearbyPeaksUpTib
longAvgPeakUpTib <- gather(avgNearbyPeaksUpTib, key = "GeneCategory", value = "AvgNumPeaks", c("avgNearbyPeaksTGFb",  "avgNearbyPeaksRA",  "avgNearbyPeaksBoth"))
ggplot(longAvgPeakUpTib, aes(x = factor(geneSet, c("subadd genes", "add genes", "mult genes", "supermult genes")), 
                             y = AvgNumPeaks, 
                             fill = factor(GeneCategory, levels = c("avgNearbyPeaksTGFb", "avgNearbyPeaksRA", "avgNearbyPeaksBoth")))) +
  geom_bar(stat="identity", position = "dodge") + 
  theme_minimal(base_size = 20) +
  labs(fill = "", x = "", y = "average number of upregulated\npeaks within 100kb of gene") +
  scale_fill_brewer(palette = "Blues")
  


### what are the new peaks getting recruiting by the "both" treatment looking like?
raSinglePeaks   <- multPeaksTib$`RA-med_isDiffPeakAtLeastTwoReps`
tgfbSinglePeaks <- multPeaksTib$`TGFb-med_isDiffPeakAtLeastTwoReps`
bothSinglePeaks <- multPeaksTib$`TGFb-and-RA-med_isDiffPeakAtLeastTwoReps`
newPeaksInBoth  <- bothSinglePeaks & !(raSinglePeaks | tgfbSinglePeaks)
newPeaksSubset <- multPeaksTib[newPeaksInBoth,]
print(multPeaksTib[newPeaksInBoth,])

# do we see more "pioneer activity" (evidenced by peaks with low etoh counts) in any category?
p50 <- qplot(subaddPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts`, xlim = c(0, 250), binwidth =5)
p51 <- qplot(addPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts`, xlim = c(0, 250), binwidth =5)
p52 <- qplot(multPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts`, xlim = c(0, 250), binwidth =5)
p53 <- qplot(supermultPeaksTib$`EtOH-nlDensity-avgNormFragmentCounts`, xlim = c(0, 250), binwidth =5)
grid.arrange(p50, p51, p52, p53, nrow=2, ncol=2)

## next: 
# # saved code snippets:

# 
# break()

subaddPeaksTib    <- addMutualExclusivityScoresToPeakTib(subaddPeaksTib)
addPeaksTib       <- addMutualExclusivityScoresToPeakTib(addPeaksTib)
multPeaksTib      <- addMutualExclusivityScoresToPeakTib(multPeaksTib)
supermultPeaksTib <- addMutualExclusivityScoresToPeakTib(supermultPeaksTib)

# histogram plots of mutual exclusivity scores
grid.arrange(qplot(subaddPeaksTib$PeakMutualExclusivityScoreAdditive) + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near sub-additive genes') + theme_classic(base_size = 16),
             qplot(addPeaksTib$PeakMutualExclusivityScoreAdditive) + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near additive genes') + theme_classic(base_size = 16),
             qplot(multPeaksTib$PeakMutualExclusivityScoreAdditive) + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near multiplicative genes') + theme_classic(base_size = 16),
             qplot(supermultPeaksTib$PeakMutualExclusivityScoreAdditive)  + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near super-multiplicative genes') + theme_classic(base_size = 16), 
             ncol = 4)

grid.arrange(qplot(subaddPeaksTib$PeakMutualExclusivityScoreAdditive, subaddPeaksTib$`TGFb-and-RA-med-avgFoldchange`) + xlim(0.5, 6.5) + ylim(0, 20),
             qplot(addPeaksTib$PeakMutualExclusivityScoreAdditive, addPeaksTib$`TGFb-and-RA-med-avgFoldchange`) + xlim(0.5, 6.5) + ylim(0, 20),
             qplot(multPeaksTib$PeakMutualExclusivityScoreAdditive, multPeaksTib$`TGFb-and-RA-med-avgFoldchange`) + xlim(0.5, 6.5) + ylim(0, 20),
             qplot(supermultPeaksTib$PeakMutualExclusivityScoreAdditive, supermultPeaksTib$`TGFb-and-RA-med-avgFoldchange`)  + xlim(0.5, 6.5) + ylim(0, 20), 
             ncol = 4)


subaddPeaksTib1     <- subaddPeaksTib %>% filter(`RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1)
addPeaksTib1        <- addPeaksTib %>% filter(`RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1)
multPeaksTib1       <- multPeaksTib %>% filter(`RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1)
supermultPeaksTib1  <- supermultPeaksTib %>% filter(`RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1)

# histogram plots of mutual exclusivity scores but only looking at peaks where both signals upregulate the gene
grid.arrange(qplot(subaddPeaksTib1$PeakMutualExclusivityScoreAdditive) + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near sub-additive genes') + theme_classic(base_size = 16),
             qplot(addPeaksTib1$PeakMutualExclusivityScoreAdditive) + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near additive genes') + theme_classic(base_size = 16),
             qplot(multPeaksTib1$PeakMutualExclusivityScoreAdditive) + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near multiplicative genes') + theme_classic(base_size = 16),
             qplot(supermultPeaksTib1$PeakMutualExclusivityScoreAdditive)  + xlim(0.5, 6.5) + xlab('RA/TGFb mutual exclusivity score,\npeaks near super-multiplicative genes') + theme_classic(base_size = 16), 
             ncol = 4)

grid.arrange(qplot(subaddPeaksTib1$PeakMutualExclusivityScoreAdditive, subaddPeaksTib1$`TGFb-and-RA-med-avgFoldchange`) + xlim(0.5, 6.5) + ylim(0, 20),
             qplot(addPeaksTib1$PeakMutualExclusivityScoreAdditive, addPeaksTib1$`TGFb-and-RA-med-avgFoldchange`) + xlim(0.5, 6.5) + ylim(0, 20),
             qplot(multPeaksTib1$PeakMutualExclusivityScoreAdditive, multPeaksTib1$`TGFb-and-RA-med-avgFoldchange`) + xlim(0.5, 6.5) + ylim(0, 20),
             qplot(supermultPeaksTib1$PeakMutualExclusivityScoreAdditive, supermultPeaksTib1$`TGFb-and-RA-med-avgFoldchange`)  + xlim(0.5, 6.5) + ylim(0, 20), 
             ncol = 4)


npeaksAdd <- length(addPeaksTib$TGFB_motifMatchInSet)
peakMatchFracTGFbAdd <- sum(addPeaksTib$TGFB_motifMatchInSet) / npeaksAdd
peakMatchFracRAAdd   <- sum(addPeaksTib$RA_motifMatchInSet) / npeaksAdd
peakMatchFracBothAdd <- sum(addPeaksTib$RA_motifMatchInSet & addPeaksTib$TGFB_motifMatchInSet) / npeaksAdd
peakMatchFracRarAndAp1Add <- sum(addPeaksTib$RAR_motifMatchInSet & addPeaksTib$AP1_motifMatchInSet) / npeaksAdd

npeaksMult <- length(multPeaksTib$ensg)
peakMatchFracTGFbMult <- sum(multPeaksTib$TGFB_motifMatchInSet) / npeaksMult
peakMatchFracRAMult   <- sum(multPeaksTib$RA_motifMatchInSet) / npeaksMult
peakMatchFracBothMult <- sum(multPeaksTib$TGFB_motifMatchInSet & multPeaksTib$RA_motifMatchInSet) / npeaksMult
peakMatchFracRarAndAp1Mult <- sum(multPeaksTib$RAR_motifMatchInSet & multPeaksTib$AP1_motifMatchInSet) / npeaksMult