library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  full.rna.seq.table <- read_tsv(here('extractedData', 'si2-si4_RNA-seq-pipeline-output-normalized.tsv'))
  counts.output.matrix <- here('extractedData', 'rnaSeqMatrixFormatted', 'counts.RNA-seq-matrix.min-count-filtered.rds')
  rpm.output.matrix    <- here('extractedData', 'rnaSeqMatrixFormatted', 'rpm.RNA-seq-matrix.min-count-filtered.rds')
  tpm.output.matrix    <- here('extractedData', 'rnaSeqMatrixFormatted',  'tpm.RNA-seq-matrix.min-count-filtered.rds')
} else {
  full.rna.seq.table <- read_tsv(cmdargs[1])
  outputFolder   <- cmdargs[2]
  
  counts.output.matrix <- paste0(outputFolder, '/counts.RNA-seq-matrix.min-count-filtered.rds')
  rpm.output.matrix    <- paste0(outputFolder, '/rpm.RNA-seq-matrix.min-count-filtered.rds')
  tpm.output.matrix    <- paste0(outputFolder, '/tpm.RNA-seq-matrix.min-count-filtered.rds')
}

N.SAMPLES.TOTAL <- 37 # 36 "core" samples plus one additional technical replicate (sample 46; EtOH control)
MINIMUM.COUNTS.ACROSS.ALL.SAMPLES <- N.SAMPLES.TOTAL * 2 #this requires at least 2 counts per sample on average for a gene to be "saved" in the final matrix
# use hist(rowSums(count.matrix), breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,24,30,36, 48, 60 , 72, 100, 5000000), xlim = c(0, 100), ylim = c(0, 0.05)) to look at histograms

wide.rna.seq.counts <- full.rna.seq.table %>% dplyr::select(sampleID, counts, gene_id) %>% spread(key=sampleID, value=counts)
wide.rna.seq.rpm    <- full.rna.seq.table %>% dplyr::select(sampleID, rpm, gene_id) %>% spread(key=sampleID, value=rpm)
wide.rna.seq.tpm    <- full.rna.seq.table %>% dplyr::select(sampleID, tpm, gene_id) %>% spread(key=sampleID, value=tpm)

names.removed.rna.seq.counts <- wide.rna.seq.counts %>% dplyr::select(-gene_id)
names.removed.rna.seq.rpm    <- wide.rna.seq.rpm    %>% dplyr::select(-gene_id)
names.removed.rna.seq.tpm    <- wide.rna.seq.tpm    %>% dplyr::select(-gene_id)
# convert fields of interest to matrix
count.matrix           <- as.matrix(names.removed.rna.seq.counts)
rownames(count.matrix) <- wide.rna.seq.counts$gene_id
rpm.matrix             <- as.matrix(names.removed.rna.seq.rpm)
rownames(rpm.matrix)   <- wide.rna.seq.rpm$gene_id
tpm.matrix             <- as.matrix(names.removed.rna.seq.tpm)
rownames(tpm.matrix)   <- wide.rna.seq.tpm$gene_id
# filter to minimum row sum
gene.indices.meeting.count.threshold <- rowSums(count.matrix) >= MINIMUM.COUNTS.ACROSS.ALL.SAMPLES
mrna.count.matrix <- count.matrix[gene.indices.meeting.count.threshold, ]
rpm.matrix        <- rpm.matrix[gene.indices.meeting.count.threshold, ]
tpm.matrix        <- tpm.matrix[gene.indices.meeting.count.threshold, ]
saveRDS(mrna.count.matrix, file = counts.output.matrix)
saveRDS(rpm.matrix, file = rpm.output.matrix)
saveRDS(tpm.matrix, file = tpm.output.matrix)