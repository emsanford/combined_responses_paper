library(tidyverse)
library(chromVAR)
library(GenomicRanges)
library(here)

mappingtib <- read_tsv(here('refs', 'EnsgHgncSymbolMapping.tsv'))
deseqtib <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.tsv'))
tpmtib <- select(deseqtib, matches('(med|nlDensity)_avgTPM|ensg|gene_name'))


viewTPM <- function(genesymbol, lookuptable = tpmtib) {
  filt.table <- filter(lookuptable, gene_name == genesymbol)
  stopifnot(nrow(filt.table) == 1)
  etohTPM <- filt.table[1, "EtOH-nlDensity_avgTPM"]
  raTPM   <- filt.table[1, "TGFb-med_avgTPM"]
  tgfbTPM <- filt.table[1, "RA-med_avgTPM"]
  bothTPM <- filt.table[1, "TGFb-and-RA-med_avgTPM"]
  print(genesymbol)
  # print(filt.table[1, "ensg"])
  print(sprintf("EtOH TPM: %0.2f", etohTPM))
  print(sprintf("RA   TPM: %0.2f", raTPM))
  print(sprintf("TGFb TPM: %0.2f", tgfbTPM))
  print(sprintf("Both TPM: %0.2f", bothTPM))
  
  return(list(etohTPM, raTPM, tgfbTPM, bothTPM))
}

# viewTPM("GRHL1")
# viewTPM("GRHL2")
# viewTPM("HOXC9")
# viewTPM("KLF7")
hox_inds <- which(grepl("HOX", deseqtib$gene_name))
hox_gene_names <- deseqtib$gene_name[hox_inds]
for (hoxgene in hox_gene_names) {
  viewTPM(hoxgene)
}
