library("tidyverse")
library("DESeq2")
library(GenomicRanges)

controlCondition <- "EtOH-nlDensity"
sampleMetadata <- read_tsv(here('sampleMetadata_SI2-SI4.txt'), col_names=TRUE)
qval <- 0.05
min.abs.log2fc <- 0.5
outputFile <- here('extractedData', 'DeSeqOutputAllConds.tsv')

blacklistedSamples <- c('46-EtOH-nlDensity')  # '46-EtOH-nlDensity' is a technical replicate of sample 15. remove for differential expression analysis as including it could artificially decrease variance

count.matrix <- readRDS(here('extractedData', 'rnaSeqMatrixFormatted', 'counts.RNA-seq-matrix.min-count-filtered.rds'))
count.matrix <- count.matrix[, !colnames(count.matrix) %in% blacklistedSamples]
relevantMetadata <- filter(sampleMetadata, !sampleID %in% blacklistedSamples) %>% filter(sampleID %in% colnames(count.matrix))
dds <- DESeqDataSetFromMatrix(countData = count.matrix,
                              colData = relevantMetadata,
                              design= ~ replicate + condition)
dds <- DESeq(dds)


deSeqTibble <- -1
for (experimentalCondition in unique(relevantMetadata$condition)) {
  if (experimentalCondition != controlCondition) {
    res <- results(dds, contrast=c("condition", experimentalCondition, controlCondition))
    isDeGene <- (!is.na(res$padj) & (abs(res$log2FoldChange) >= min.abs.log2fc) & (res$padj <= qval))  # note that NA values only show up for padj values in genes with very low expression
    condition <- rep(experimentalCondition, length(isDeGene))
    restib <- tibble(condition=condition, baseMean=res$baseMean, log2fc=res$log2FoldChange, pval=res$pvalue, padj=res$padj, ensg=rownames(res), isDeGene=isDeGene)

    print(sprintf("num differentially expressed genes in %s with current parameters: %d", experimentalCondition, sum(isDeGene)))
    
    if (deSeqTibble == -1) {
      deSeqTibble <- restib
    } else {
      deSeqTibble <- rbind(deSeqTibble, restib)
    }
  }
}

### make the "long" deSeqTibble "wide" (want one entry per gene)
myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

wideDeSeqTib <- deSeqTibble %>% dplyr::select(-baseMean, -pval) %>% myspread(condition, c(padj, log2fc, isDeGene))

# we want to get the average TPM across replicates for a condition
pipeline.output.table <- read_tsv(here('extractedData', 'si2-si4_RNA-seq-pipeline-output-normalized.tsv'), col_names=T)

pipeline.output.table.with.condition <- relevantMetadata %>% 
  dplyr::select(sampleID, condition) %>%
  inner_join(pipeline.output.table) %>%
  ungroup()

pipeline.output.table.with.avgTPM <- pipeline.output.table.with.condition %>%
  group_by(condition, gene_id) %>%
  mutate(avgTPM = mean(tpm), avgRPM = mean(rpm), avgCounts = mean(counts)) %>%
  ungroup()

pipeline.output.table.wide <- pipeline.output.table.with.avgTPM %>% 
                            dplyr::select(avgTPM, avgRPM, avgCounts, gene_id, condition) %>% 
                            unique() %>% mutate(ensg = gene_id) %>% dplyr::select(-gene_id) %>%
                            myspread(condition, c(avgTPM, avgRPM, avgCounts))

# now join normalized count values back to table with differential gene expression
wideDeSeqTibWithNlCountData <- left_join(wideDeSeqTib, pipeline.output.table.wide, by = 'ensg')

# now add replicate-level info
pipeline.output.table.rmexpt <- pipeline.output.table %>% mutate(ensg=gene_id) %>% dplyr::select(-gene_id, -experiment, -rpkm)
pipeline.output.table.rmexpt.wide <- myspread(pipeline.output.table.rmexpt, sampleID, c(counts, rpm, tpm))
wideDeSeqTibWithNlCountData <- left_join(wideDeSeqTibWithNlCountData, pipeline.output.table.rmexpt.wide, by='ensg')
# for (sid in unique(pipeline.output.table$sampleID)) {
#   tibToJoin <- pipeline.output.table %>% filter(sampleID == sid) %>% mutate(ensg = gene_id) %>% dplyr::select(-gene_id)
#   tibToJoin
#   
#   wideDeSeqTibWithNlCountData <- left_join(wideDeSeqTibWithNlCountData, )
# }

write_tsv(wideDeSeqTibWithNlCountData, outputFile, col_names = T)





