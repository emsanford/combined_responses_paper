library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)
library(chromVARmotifs)
library(BiocParallel)
library(motifmatchr)
library(SummarizedExperiment)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", version = "3.8")
library(BSgenome.Hsapiens.UCSC.hg38)
register(MulticoreParam(4, progressbar = TRUE))

# manually_include_these_motifs <- c("ENSG00000126778_LINE2315_SIX1_I", "ENSG00000102974_LINE747_CTCF_D_N67")  # from homer de novo motif analysis on upregulated peaks
# manually_include_these_motifs <- c("SIX1", "CTCF")  # from homer de novo motif analysis on upregulated peaks
manually_include_these_motifs <- c("CTCF")

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  fragmentCounts <- read_rds("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/final_diffPeaks_fragment_counts_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
  outputFile <- here("extractedData", "most_variable_motifs_Robject.rds")
} else {
  fragmentCounts <- read_rds(cmdargs[1])
  outputFile <- cmdargs[2]
}

# select the most variable motifs from the data set (selected by finding "natural looking cutoff" in the full set of curated motifs)
fragmentCounts <- addGCBias(fragmentCounts, 
                            genome = BSgenome.Hsapiens.UCSC.hg38)

set.seed(2019)
data("human_pwms_v2") #loads the curated cisBP motif set from the chromVar paper
motifSet <- human_pwms_v2 #loads the curated cisBP motif set from the chromVar paper
motif_ix <- matchMotifs(motifSet, fragmentCounts,
                        genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object = fragmentCounts, annotations = motif_ix)
variability <- computeVariability(dev)
p <- plotVariability(variability, use_plotly = FALSE)
top_n_motifs_to_keep <- 50  # this was chosen by manually looking at the plot below to see where the natural break was
p <- p + geom_vline(xintercept = top_n_motifs_to_keep)
# print(p)
tvar <- as.tibble(variability)

TF.names.variabilitycutoff <- arrange(tvar, -variability)$name[1:top_n_motifs_to_keep]
TF.names.manual.additions  <- tvar %>% filter(name %in% manually_include_these_motifs) %>% pull("name")
TF.names.variabilitycutoff.plus.manual.additions <- union(TF.names.variabilitycutoff, TF.names.manual.additions)
PWM.object.indices.bool <- sapply(strsplit(names(motifSet), "_"), function (x) x[3]) %in% TF.names.variabilitycutoff.plus.manual.additions
PWM.object.indices.num  <- which(PWM.object.indices.bool)
selected.PWM.objects <- motifSet[PWM.object.indices.num]
saveRDS(selected.PWM.objects, file = outputFile)
################################################################################################################  

