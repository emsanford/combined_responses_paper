deseqTibAnno <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))

deseqTibAnno %>% 
  filter(gene_name == 'MAP3K1') %>%
  dplyr::select(matches("integrationConstant"))

deseqTibAnno %>% 
  filter(gene_name == 'GPRC5A') %>%
  dplyr::select(matches("integrationConstant"))

deseqTibAnno %>% 
  filter(gene_name == 'EPHB2') %>%
  dplyr::select(matches("integrationConstant"))

deseqTibAnno %>% 
  filter(gene_name == 'RIPK4') %>%
  dplyr::select(matches("integrationConstant"))

deseqTibAnno %>% 
  filter(gene_name == 'ZNF469') %>%
  dplyr::select(matches("integrationConstant"))