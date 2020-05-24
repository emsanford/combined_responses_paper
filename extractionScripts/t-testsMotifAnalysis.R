library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)
library(chromVARmotifs)
library(BiocParallel)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
source(here('extractionScripts', 'util.R'))

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  upreg.peaks          <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  mostVariableMotifSet <- read_rds("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/mostVariableMotifs_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
  outputPlotPrefix     <- here("plots", "")
} else {
  upreg.peaks           <- read_tsv(cmdargs[1])
  mostVariableMotifSet <- read_rds(cmdargs[2])
  outputPlotPrefix     <- cmdargs[3]
}

n.bootstrap.samples <- 1000

#0 define motif set (top 50 or all in cisBP data set)
data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
cisbp_motifs <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper

calcAvgNumMotifsPerPeakType <- function(upreg.peaks.granges, motif_ix, motifsetname) {
  restib <- NULL
  
  # 3: count number of motif matches 
  motif.count.matx <- as.matrix(motifCounts(motif_ix))
  motif.sum.each.peak <- rowSums(motif.count.matx)
  peak.types <- names(motif.sum.each.peak)
  
  inds.subadditive   <- peak.types == "sub-additive"
  inds.additive      <- peak.types == "additive"
  inds.superadditive <- peak.types == "super-additive"
  
  peak_sizes <- width(ranges(upreg.peaks.granges))
  subadd.total.peak.size     <- sum(peak_sizes[inds.subadditive])  
  add.total.peak.size        <- sum(peak_sizes[inds.additive])  
  superadd.total.peak.size   <- sum(peak_sizes[inds.superadditive]) 
  
  subadd.avg.peak.size     <- sum(peak_sizes[inds.subadditive]) / sum(inds.subadditive)
  add.avg.peak.size        <- sum(peak_sizes[inds.additive]) / sum(inds.additive)  
  superadd.avg.peak.size   <- sum(peak_sizes[inds.superadditive]) / sum(inds.superadditive)  
  
  subadd.num.matches.per.150.peak.bp   <- 150 * sum(motif.sum.each.peak[inds.subadditive])   / subadd.total.peak.size
  add.num.matches.per.150.peak.bp      <- 150 * sum(motif.sum.each.peak[inds.additive])      / add.total.peak.size
  superadd.num.matches.per.150.peak.bp <- 150 * sum(motif.sum.each.peak[inds.superadditive]) / superadd.total.peak.size
  
  # 4: store results in result tibble
  restib <- tibble(avg_motifs_per_peak = c(subadd.num.matches.per.150.peak.bp, add.num.matches.per.150.peak.bp, superadd.num.matches.per.150.peak.bp), 
                   avg_peak_size       = c(subadd.avg.peak.size, add.avg.peak.size, superadd.avg.peak.size),
                   peak_type = factor(c("sub-additive", "additive", "super-additive"), levels = c("sub-additive", "additive", "super-additive")), 
                   motif_set = rep(motifsetname, 3))
  
  return(list(restib, peak_sizes, motif.sum.each.peak, peak.types))
}

# 1: make genomic ranges for upreg peaks (medium dose)

peaktypes        <- sapply(upreg.peaks$`peak_integrationCategory-med-dose`, convertUpregCvalCatToDvalCat)
names(peaktypes) <- peaktypes

upreg.peaks.granges <- GRanges(seqnames = upreg.peaks$chrom,
                               ranges = IRanges(start = upreg.peaks$startLocs,
                                                end   = upreg.peaks$endLocs),
                               peaktype = peaktypes)

# 2: match motifs to these genomic ranges (may need hg38 reference seq). pick full motif set or top50 motifs
motif_ix_mostVariable <- matchMotifs(mostVariableMotifSet, upreg.peaks.granges, genome = BSgenome.Hsapiens.UCSC.hg38, out = "scores")
motif_ix_allCisBP     <- matchMotifs(cisbp_motifs, upreg.peaks.granges, genome = BSgenome.Hsapiens.UCSC.hg38, out = "scores")

result.list.topmotifs <- calcAvgNumMotifsPerPeakType(upreg.peaks.granges, motif_ix_mostVariable, "most variable motifs")
result.list.allmotifs <- calcAvgNumMotifsPerPeakType(upreg.peaks.granges, motif_ix_allCisBP, "all cisBP motifs")
peak.sizes            <- result.list.allmotifs[[2]]

inds.subadditive   <- result.list.topmotifs[[4]] == "sub-additive"
inds.additive      <- result.list.topmotifs[[4]] == "additive"
inds.superadditive <- result.list.topmotifs[[4]] == "super-additive"

print("all cisBP motifs: t-test sub-additive vs. additive peak motif density")
subadd.motif.densities   <- result.list.allmotifs[[3]][inds.subadditive]   / peak.sizes[inds.subadditive]
add.motif.densities      <- result.list.allmotifs[[3]][inds.additive]      / peak.sizes[inds.additive]
superadd.motif.densities <- result.list.allmotifs[[3]][inds.superadditive] / peak.sizes[inds.superadditive]
print(t.test(subadd.motif.densities, add.motif.densities))
print("all cisBP motifs: t-test super-additive vs. additive peak motif density")
print(t.test(superadd.motif.densities, add.motif.densities))

print("enriched motifs: t-test sub-additive vs. additive peak motif density")
subadd.motif.densities   <- result.list.topmotifs[[3]][inds.subadditive]   / peak.sizes[inds.subadditive]
add.motif.densities      <- result.list.topmotifs[[3]][inds.additive]      / peak.sizes[inds.additive]
superadd.motif.densities <- result.list.topmotifs[[3]][inds.superadditive] / peak.sizes[inds.superadditive]
print(t.test(subadd.motif.densities, add.motif.densities))
print("enriched motifs: t-test super-additive vs. additive peak motif density")
print(t.test(superadd.motif.densities, add.motif.densities))



# Now do t-tests for individual motif matches
ind.motif.counts <- motifCounts(motif_ix_mostVariable)
smarcc1.motif.name <- "ENSG00000173473_LINE2867_SMARCC1_D_N2"
smarcc1.densities.subadditive <- ind.motif.counts[inds.subadditive, smarcc1.motif.name] / peak.sizes[inds.subadditive]
smarcc1.densities.additive    <- ind.motif.counts[inds.additive,    smarcc1.motif.name] / peak.sizes[inds.additive]
smarcc1.densities.superadditive    <- ind.motif.counts[inds.superadditive, smarcc1.motif.name] / peak.sizes[inds.superadditive]
print("smarcc1 motif density t test, subadditive vs. additive")
print(t.test(smarcc1.densities.subadditive, smarcc1.densities.additive))
# print(t.test(smarcc1.densities.superadditive, smarcc1.densities.additive))

jun.motif.name     <- "ENSG00000177606_LINE517_JUN_D_N5"
jun.densities.subadditive <- ind.motif.counts[inds.subadditive, jun.motif.name] / peak.sizes[inds.subadditive]
jun.densities.additive    <- ind.motif.counts[inds.additive,    jun.motif.name] / peak.sizes[inds.additive]
jun.densities.superadditive    <- ind.motif.counts[inds.superadditive, jun.motif.name] / peak.sizes[inds.superadditive]
print("jun motif density t test, subadditive vs. additive")
print(t.test(jun.densities.subadditive, jun.densities.additive))

fos.motif.name     <- "ENSG00000170345_LINE476_FOS_D_N4"
fos.densities.subadditive <- ind.motif.counts[inds.subadditive, fos.motif.name] / peak.sizes[inds.subadditive]
fos.densities.additive    <- ind.motif.counts[inds.additive,    fos.motif.name] / peak.sizes[inds.additive]
fos.densities.superadditive    <- ind.motif.counts[inds.superadditive, fos.motif.name] / peak.sizes[inds.superadditive]
print("fos motif density t test, subadditive vs. additive")
print(t.test(fos.densities.subadditive, fos.densities.additive))

ctcf.motif.name    <- "ENSG00000102974_LINE747_CTCF_D_N67"
ctcf.densities.subadditive <- ind.motif.counts[inds.subadditive, ctcf.motif.name] / peak.sizes[inds.subadditive]
ctcf.densities.additive    <- ind.motif.counts[inds.additive,    ctcf.motif.name] / peak.sizes[inds.additive]
ctcf.densities.superadditive    <- ind.motif.counts[inds.superadditive, ctcf.motif.name] / peak.sizes[inds.superadditive]
print("ctcf motif density t test, subadditive vs. additive")
print(t.test(ctcf.densities.subadditive, ctcf.densities.additive))
print("ctcf motif density t test, subadditive vs. superadditive")
print(t.test(ctcf.densities.subadditive, ctcf.densities.superadditive))

smad3.motif.name   <- "ENSG00000166949_LINE3324_SMAD3_D_N1"
smad3.densities.subadditive <- ind.motif.counts[inds.subadditive, smad3.motif.name] / peak.sizes[inds.subadditive]
smad3.densities.additive    <- ind.motif.counts[inds.additive,    smad3.motif.name] / peak.sizes[inds.additive]
smad3.densities.superadditive    <- ind.motif.counts[inds.superadditive, smad3.motif.name] / peak.sizes[inds.superadditive]
print("smad3 motif density t test, subadditive vs. superadditive")
print(t.test(smad3.densities.additive, smad3.densities.superadditive))

smad4.motif.name   <- "ENSG00000141646_LINE3314_SMAD4_D"
smad4.densities.subadditive <- ind.motif.counts[inds.subadditive, smad4.motif.name] / peak.sizes[inds.subadditive]
smad4.densities.additive    <- ind.motif.counts[inds.additive,    smad4.motif.name] / peak.sizes[inds.additive]
smad4.densities.superadditive    <- ind.motif.counts[inds.superadditive, smad4.motif.name] / peak.sizes[inds.superadditive]
print("smad4 motif density t test, subadditive vs. superadditive")
print(t.test(smad4.densities.additive, smad4.densities.superadditive))

smad9.motif.name   <- "ENSG00000120693_LINE19375_SMAD9_I_N1"
smad9.densities.subadditive <- ind.motif.counts[inds.subadditive, smad9.motif.name] / peak.sizes[inds.subadditive]
smad9.densities.additive    <- ind.motif.counts[inds.additive,    smad9.motif.name] / peak.sizes[inds.additive]
smad9.densities.superadditive    <- ind.motif.counts[inds.superadditive, smad9.motif.name] / peak.sizes[inds.superadditive]
print("smad9 motif density t test, subadditive vs. superadditive")
print(t.test(smad9.densities.additive, smad9.densities.superadditive))

elf1.motif.name    <- "ENSG00000120690_LINE1828_ELF1_D_N8"
elf1.densities.subadditive <- ind.motif.counts[inds.subadditive, elf1.motif.name] / peak.sizes[inds.subadditive]
elf1.densities.additive    <- ind.motif.counts[inds.additive,    elf1.motif.name] / peak.sizes[inds.additive]
elf1.densities.superadditive    <- ind.motif.counts[inds.superadditive, elf1.motif.name] / peak.sizes[inds.superadditive]
print("elf1 motif density t test, subadditive vs. superadditive")
print(t.test(elf1.densities.additive, elf1.densities.superadditive))


