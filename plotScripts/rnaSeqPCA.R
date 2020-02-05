library(here)
library(tidyverse)
sample.metadata <- read_tsv(here('extractedData', 'sampleMetadata_SI2-SI4.txt'))
rna.count.matrix <- readRDS(here('extractedData', 'rnaSeqMatrixFormatted', 'counts.RNA-seq-matrix.min-count-filtered.rds'))
rpm.matrix <- readRDS(here('extractedData', 'rnaSeqMatrixFormatted', 'rpm.RNA-seq-matrix.min-count-filtered.rds'))
tpm.matrix <- readRDS(here('extractedData', 'rnaSeqMatrixFormatted', 'tpm.RNA-seq-matrix.min-count-filtered.rds'))

# # test effect of including relevant HD3 samples; old code, need to generate TPM matrix files first...
# hd3_rnaseq_tib <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/HistoryDependence/SequencingResults/Analysis_HD3_to_HD6/extractedData/rnaSeqSampleGeneCounts.tsv", col_names = T)
# relevant_fields <- hd3_rnaseq_tib %>% filter(experiment=='HD3_RNA-Seq', sampleID %in% c('04-TGFb-3d-rep1', '06-RA-3d-rep1'))
# filt_genes <- relevant_fields %>% filter(gene_id %in% rownames(rna.count.matrix))
# filt_gene_sprd <- spread(filt_genes, key=sampleID, value=counts)
# filt_gene_matx <- as.matrix(filt_gene_sprd[, c(3,4)])
# rownames(filt_gene_matx) <- filt_gene_sprd$gene_id
# rna.count.matrix <- cbind(rna.count.matrix, filt_gene_matx)
# sample.metadata <- rbind(sample.metadata, c('04-TGFb-3d-rep1', '04-TGFb-3d-rep1', '04-TGFb-3d-rep1', '04-TGFb-3d-rep1', '04-TGFb-3d-rep1'), c('06-RA-3d-rep1', '06-RA-3d-rep1', '06-RA-3d-rep1', '06-RA-3d-rep1', '06-RA-3d-rep1'))

sample.coverage <- colSums(rna.count.matrix)
# #### this block of code creates a new "summed" sample from the SI3 rerun for SI3 original samples 17 and 18 ####################################
# rna.count.matrix <- cbind(rna.count.matrix, rna.count.matrix[, '49-RA-med-second-RNA-extraction'] + rna.count.matrix[, '47-RA-med'])
# rna.count.matrix <- cbind(rna.count.matrix, rna.count.matrix[, '50-TGFb-and-RA-med-second-RNA-extraction'] + rna.count.matrix[, '48-TGFb-and-RA-med'])
# colnames(rna.count.matrix)[40:41] <- c('47+49-RA-med', '48+50-TGFb-and-RA-med')
# sample.metadata <- rbind(sample.metadata, c('47+49-RA-med', NA, NA, NA, NA, NA, NA, NA))
# sample.metadata <- rbind(sample.metadata, c('48+50-TGFb-and-RA-med', NA, NA, NA, NA, NA, NA, NA))
# ################################################################################################################################################

####### select a single condition to look at for PCA #######
## conditionSelect <- c('RA-low', 'RA-med', 'RA-high')
## conditionSelect <- c('TGFb-low', 'TGFb-med', 'TGFb-high')
# conditionSelect <- c('TGFb-and-RA-low', 'TGFb-and-RA-med', 'TGFb-and-RA-high')
# sampleSelectionIndexVector <- sample.metadata$condition %in% conditionSelect
# sample.metadata <- sample.metadata[sampleSelectionIndexVector,]
# rpm.matrix <- rpm.matrix[, sampleSelectionIndexVector]
###############################################################

#unscaled PCA, PC1 vs. PC2, using the rpm input matrix
pcres <- prcomp(t((rpm.matrix)), scale = F)
pcscores <- pcres$x
pcpctvars <- pcres$sdev^2 / (sum(pcres$sdev^2))
sample.metadata.matched <- filter(sample.metadata, SampleID %in% colnames(rpm.matrix))
restib <- sample.metadata.matched %>% mutate(pc1score = pcscores[,"PC1"], pc2score = pcscores[,"PC2"], pc3score = pcscores[,"PC3"], pc4score = pcscores[,"PC4"])
restib %>% ggplot(aes(pc1score, pc2score, color=condition, label=SampleID)) + geom_text() + 
  ggtitle('Dose response RNA-seq PCA (centered, unscaled), on rpm values of minCountFiltered genes') + 
  xlab(sprintf('PC1 (%.02f pct variance)', pcpctvars[1] * 100)) +
  ylab(sprintf('PC2 (%.02f pct variance)', pcpctvars[2] * 100)) + expand_limits(x = c(-115, 140)) + theme_light(base_size = 20)

#unscaled PCA, PC3 vs. PC4, using the counts per million input matrix
restib %>% ggplot(aes(pc3score, pc4score, color=sample.coverage, label=SampleID)) + geom_text() + 
  ggtitle('Dose response RNA-seq PCA (centered, unscaled), on rpm values of minCountFiltered genes') + 
  xlab(sprintf('PC3 (%.02f pct variance)', pcpctvars[3] * 100)) +
  ylab(sprintf('PC4 (%.02f pct variance)', pcpctvars[4] * 100))


#unscaled PCA, PC1 vs. PC2, using the fold change matrix (matched to each individual control)
findTheMatchedControl <- function(sampleID) {
  sampleNumberString <- substr(sampleID,1,2)
  sampleNumber <- as.numeric(sampleNumberString)
  if (sampleNumber %>% between(1,12)) {
    return("05-EtOH-nlDensity")
  }
  if (sampleNumber %>% between(13,24)) {
    return("15-EtOH-nlDensity")
  }
  if (sampleNumber %>% between(25,36)) {
    return("27-EtOH-nlDensity")
  }
  if (sampleNumber %>% between(46,52)) {
    return("46-EtOH-nlDensity")
  } else {
    return(NA)
  }
}
matchedControlVector <- c()
for (sid in sample.metadata.matched$SampleID) {
  matchedControlVector <- c(matchedControlVector, findTheMatchedControl(sid))
}
sample.metadata.matched <- sample.metadata.matched %>% mutate(controlSampleID = matchedControlVector)

foldchange.matrix <- rpm.matrix / rpm.matrix[, sample.metadata.matched$controlSampleID]
foldchange.matrix <- foldchange.matrix[, !colnames(foldchange.matrix) %in% c("05-EtOH-nlDensity", "15-EtOH-nlDensity", "27-EtOH-nlDensity", "46-EtOH-nlDensity")]
sample.metadata.matched.without.nlcontrols <- filter(sample.metadata.matched, !SampleID %in% c("05-EtOH-nlDensity", "15-EtOH-nlDensity", "27-EtOH-nlDensity", "46-EtOH-nlDensity"))
foldchange.matrix.infinite.values.removed <- foldchange.matrix[is.finite(rowMeans(foldchange.matrix)),]
logfoldchange.matrix.infinite.values.removed <- log2(foldchange.matrix.infinite.values.removed)
logfoldchange.matrix.infinite.values.removed <- logfoldchange.matrix.infinite.values.removed[is.finite(rowMeans(logfoldchange.matrix.infinite.values.removed)),]

pcres <- prcomp(t((logfoldchange.matrix.infinite.values.removed)), scale = F)
pcscores <- pcres$x
pcpctvars <- pcres$sdev^2 / (sum(pcres$sdev^2))
restib <- sample.metadata.matched.without.nlcontrols %>% mutate(pc1score = pcscores[,"PC1"], pc2score = pcscores[,"PC2"], pc3score = pcscores[,"PC3"], pc4score = pcscores[,"PC4"])
restib %>% ggplot(aes(pc1score, pc2score, color=condition, label=SampleID)) + geom_text() + 
  ggtitle('Dose response RNA-seq PCA (centered, unscaled), on log-fold changes of minCountFiltered genes') + 
  xlab(sprintf('PC1 (%.02f pct variance)', pcpctvars[1] * 100)) +
  ylab(sprintf('PC2 (%.02f pct variance)', pcpctvars[2] * 100))




#### do PCA on fold change of TPM, normalize each to its paired control

#### differentially expressed genes only
deSeqTib <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.tsv'))
isDeGeneCols <- grepl('isDeGene', colnames(deSeqTib))
isDeGeneColumnOnly <- deSeqTib[, isDeGeneCols]
numDeCondsPerGene  <- rowSums(isDeGeneColumnOnly)
isDeGeneAnyCond <- numDeCondsPerGene > 0
pcacols <- grepl('log2fc', colnames(deSeqTib))
# pcacols <- grepl('avgTPM', colnames(deSeqTib))
pcainput <- as.matrix(deSeqTib[isDeGeneAnyCond, pcacols])
# pcainput <- 2 ^ pcainput - 1 # uncomment this line to run PCA on fold change rather than log fold change
pcares <- prcomp(t(pcainput))
pcpctvars <- pcares$sdev^2 / (sum(pcares$sdev^2))
pcscores <- pcares$x
forvis <- tibble(pc1 = pcscores[, 1], pc2 = pcscores[, 2], pc3 = pcscores[,3], pc4 = pcscores[, 4], condition = rownames(pcscores))
cond_names <- map(rownames(pcscores), function(x) strsplit(x, '_')[[1]][[1]])
ggplot(forvis, aes(pc1, pc2, label=cond_names)) + geom_text(size=8) + theme_light(base_size = 30) +
  xlab(sprintf('PC1 (%.02f pct variance)', pcpctvars[1] * 100)) +
  ylab(sprintf('PC2 (%.02f pct variance)', pcpctvars[2] * 100))
pcares$sdev^2 / sum(pcares$sdev ^2)