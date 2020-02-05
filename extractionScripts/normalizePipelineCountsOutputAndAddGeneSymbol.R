library(tidyverse)
library(GenomicFeatures)
library(here)
library(refGenome)
# note: if you get errors from any of the above "library" lines, e.g. "there is no package called ‘here’, 
#       enter the following command in the RStudio console: install.packages(<name_of_library_that_failed_to_load>)

######### here beginneth user-defined parameters #########

meltedDataInputFile   <- here('extractedData', 'si2-si4_RNA-seq-pipeline-output-counts.tsv')
gtfFileUsedByPipeline <- here('refs', 'hg38.gtf')
ENSGtoGeneSymbolTable <- here('refs', 'EnsgHgncSymbolMapping.tsv')
geneCanonicalTssTable <- read_tsv(here('refs', 'EnsgToHg38TssMapping.tsv'))
outputFile            <- here('extractedData', 'si2-si4_RNA-seq-pipeline-output-normalized.tsv')
techRepPair1toSum <- c('47-RA-med', '49-RA-med-second-RNA-extraction')
techRepPair2toSum <- c('48-TGFb-and-RA-med', '50-TGFb-and-RA-med-second-RNA-extraction')
sampleIdBlacklist <- c('17-RA-med', '18-TGFb-and-RA-med') # failed post-sequencing QC due to high rRNA contamination
gene.types.to.keep <- c("protein_coding", "lincRNA")
######### here endeth user-defined parameters ############

# parse GTF file into an R object and use a union model to calculate gene lengths
txdb <- makeTxDbFromGFF(file = gtfFileUsedByPipeline, format="gtf")
lengthsPergeneid <- sum(width(IRanges::reduce(exonsBy(txdb, by = "gene"))))
lengthtbl <- tibble(gene_id = names(lengthsPergeneid), length = lengthsPergeneid)

# read in the "melted" HTSeq count table, sum technical replicates, remove blacklisted samples
htseq.table.original <- read_tsv(meltedDataInputFile, col_names = T)
htseq.table.rm.blacklisted <- filter(htseq.table.original, !sampleID %in% sampleIdBlacklist)
htseq.table.techRepPair1Sum <- htseq.table.rm.blacklisted %>%
  filter(sampleID %in% techRepPair1toSum) %>%
  group_by(gene_id) %>%
  mutate(countSum = sum(counts)) %>%
  dplyr::select(experiment, gene_id, countSum) %>%
  mutate(counts = countSum) %>% dplyr::select(-countSum) %>%
  unique() %>%
  mutate(sampleID = '51-RA-med')

htseq.table.techRepPair2Sum <- htseq.table.rm.blacklisted %>%
  filter(sampleID %in% techRepPair2toSum) %>%
  group_by(gene_id) %>%
  mutate(countSum = sum(counts)) %>%
  dplyr::select(experiment, gene_id, countSum) %>%
  mutate(counts = countSum) %>% dplyr::select(-countSum) %>%
  unique() %>%
  mutate(sampleID = '52-TGFb-and-RA-med')

htseq.table.rm.summed.samples <- htseq.table.rm.blacklisted %>% filter(!sampleID %in% c(techRepPair1toSum, techRepPair2toSum))
htseq.table <- bind_rows(htseq.table.rm.summed.samples, htseq.table.techRepPair1Sum, htseq.table.techRepPair2Sum)

# remove genes that aren't coding or lincRNAs
ens <- ensemblGenome()
setwd(here("refs"))
read.gtf(ens, "hg38.gtf")
my_gene <- getGenePositions(ens)
gene_names_to_keep <- filter(my_gene, gene_biotype %in% gene.types.to.keep)$gene_id
htseq.table <- filter(htseq.table, gene_id %in% gene_names_to_keep)

# add the following normalized values for each sample: RPM, RPKM, and TPM
htseq.table.withLength <- inner_join(htseq.table, lengthtbl, by = 'gene_id') # this step also removes non_gene features from the table, e.g. "__no_feature"
htseq.table.withRPM <- htseq.table.withLength %>% 
  group_by(experiment, sampleID) %>%
  mutate(totalMappedReads = sum(counts), rpm = 1000000*counts / totalMappedReads)
htseq.table.withRPKM <- htseq.table.withRPM %>%
  mutate(rpkm = 1000*rpm/length)
# note this TPM calc method doesn't take average fragment size into account (some papers use it, some papers don't)
htseq.table.withTPM <- htseq.table.withRPKM %>%
  mutate(ctsOverLength = counts/length, denominator = sum(ctsOverLength), tpm = ctsOverLength/denominator * 1000000) %>%
  dplyr::select(-c(ctsOverLength, denominator))

# now add the gene symbol to the final table and write the output file
HGNC.symbol.table <- read_tsv(ENSGtoGeneSymbolTable, col_names = T)
HGNC.symbol.table <- HGNC.symbol.table %>% mutate(gene_id = ensg) %>% dplyr::select(-ensg)
htseq.table.with.HGNCsymbol <- plyr::join(data.frame(htseq.table.withTPM), data.frame(HGNC.symbol.table), by = 'gene_id', type="left", match="first")  #uses plyr join on a data frame instead of dplyr left_join on a tibble due to first match option
htseq.table.almostfinal <- as_tibble(htseq.table.with.HGNCsymbol)

# add hg38.gtf names to the genes that are missing hgnc symbols
indices.genes.missing.hgnc.symbols <- which(is.na(htseq.table.almostfinal$hgnc_symbol))
gene.ids.missing.symbols <- htseq.table.almostfinal$gene_id[indices.genes.missing.hgnc.symbols]
temp.table1 <- data.frame(tibble(gene_id = gene.ids.missing.symbols))
temp.table2 <- data.frame(tibble(gene_id = my_gene$gene_id, name = my_gene$gene_name))
temp.table <- plyr::join(temp.table1, temp.table2, by = 'gene_id', type="left", match="first") 
htseq.table.almostfinal$hgnc_symbol[indices.genes.missing.hgnc.symbols] <- temp.table$name
stopifnot(sum(is.na(htseq.table.almostfinal$hgnc_symbol)) == 0)
# now rename HGNC column to "gene_name" column. Gene names are the HGNC symbol if one exists and the gene_name in the gtf file otherwise.
htseq.table.almostfinal[["gene_name"]] <- htseq.table.almostfinal$hgnc_symbol
htseq.table.almostfinal <- dplyr::select(htseq.table.almostfinal, -hgnc_symbol)

# add TSS from canonical file then add 
htseq.table.almostfinal <- left_join(htseq.table.almostfinal, dplyr::select(geneCanonicalTssTable, -transcript_id), by = "gene_id")
indices.genes.missing.tss <- which(is.na(htseq.table.almostfinal$tx_start))
gene.ids.missing.tss <- htseq.table.almostfinal$gene_id[indices.genes.missing.tss]
temp.table1 <- data.frame(tibble(gene_id = gene.ids.missing.tss))
temp.table2 <- data.frame(tibble(gene_id  = my_gene$gene_id,
                                 chrom    = my_gene$seqid,
                                 tx_start = my_gene$start,
                                 tx_end   = my_gene$end,
                                 strand   = my_gene$strand))
temp.table <- plyr::join(temp.table1, temp.table2, by = 'gene_id', type="left", match="first") 

htseq.table.almostfinal$chrom[indices.genes.missing.tss] <- temp.table$chrom
htseq.table.almostfinal$tx_start[indices.genes.missing.tss] <- temp.table$tx_start
htseq.table.almostfinal$tx_end[indices.genes.missing.tss] <- temp.table$tx_end
htseq.table.almostfinal$strand[indices.genes.missing.tss] <- temp.table$strand
# add TSS coordinate column, which depends on strand, to TSS location table
htseq.table.almostfinal[["TSS_loc"]] <- rep(0, nrow(htseq.table.almostfinal))
htseq.table.almostfinal[["TSS_loc"]][which(htseq.table.almostfinal$strand == "+")] <- htseq.table.almostfinal[["tx_start"]][which(htseq.table.almostfinal$strand == "+")]
htseq.table.almostfinal[["TSS_loc"]][which(htseq.table.almostfinal$strand == "-")] <- htseq.table.almostfinal[["tx_end"]][which(htseq.table.almostfinal$strand == "-")]

stopifnot(sum(is.na(htseq.table.almostfinal$chrom)) == 0)


htseq.table.final <- htseq.table.almostfinal %>% dplyr::select(-length, -totalMappedReads)
write_tsv(htseq.table.final, outputFile, col_names = T)


####### if you need to make a new ENSGtoGeneSymbolTable, uncomment, edit, and run this block of code ###############
# library(biomaRt)
# library(tidyverse)
# library(GenomicFeatures)
# library(here)
# gtfFileUsedByPipeline <- '/Users/emsanford/Downloads/Homo_sapiens.GRCh38.97.gtf'
# txdb <- makeTxDbFromGFF(file = gtfFileUsedByPipeline, format="gtf")
# vectorOfENSG_IDs_to_find_HGNC_symbols_for <- names(sum(width(IRanges::reduce(exonsBy(txdb2, by = "gene")))))  ## this is just a vector of "ENSG0000XXXX" strings. This line of code gives you this vector from a parsed GTF file that you make below
# outputTableWithGeneNameSymbols <- here('newGeneNameSymbolTable.tsv') 
# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "uswest.ensembl.org", ensemblRedirect = FALSE)
# hgnc.name.table <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
#                          filters = 'ensembl_gene_id', 
#                          values = vectorOfENSG_IDs_to_find_HGNC_symbols_for, 
#                          mart = ensembl)
# 
# hgnc.name.tibble <- tibble(ensg = hgnc.name.table$ensembl_gene_id, hgnc_symbol = hgnc.name.table$hgnc_symbol)
# write_tsv(hgnc.name.tibble, outputTableWithGeneNameSymbols, col_names = TRUE)
####################################################################################################################