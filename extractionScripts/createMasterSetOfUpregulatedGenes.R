library(tidyverse)
library(here)

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  geneTibAnno <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
  outputLoc   <- here("extractedData", "DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
} else {
  geneTibAnno  <- read_tsv(cmdargs[1])
  outputLoc <- cmdargs[2]
}



# this definition is "upregulated by combined treatment OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals
#    requires that the individual RA and TGFb effects be positive at all doses
siUpregGenes <- geneTibAnno %>% filter(`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene` | ((`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
                                       `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
                                       `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
                                       `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0)

write_tsv(siUpregGenes, outputLoc, col_names = T)


# how does the requirement of indiviudal increases in expression affect the number of genes in the master set?
mastersetgeneIDs <- siUpregGenes$ensg

ra.diff.gene.ids <- filter(geneTibAnno, `RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) %>% pull("ensg")
ra.upreg.gene.ids <- filter(geneTibAnno, `RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`, `RA-high_log2fc` > 0) %>% pull("ensg") 

tgfb.diff.gene.ids <- filter(geneTibAnno, `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`) %>% pull("ensg")
tgfb.upreg.gene.ids <- filter(geneTibAnno, `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`, `TGFb-high_log2fc` > 0) %>% pull("ensg") 

both.diff.gene.ids <- filter(geneTibAnno, `TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`) %>% pull("ensg")
both.upreg.gene.ids <- filter(geneTibAnno, `TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`, `TGFb-and-RA-high_log2fc` > 0) %>% pull("ensg") 
downreg.both.gene.ids <- setdiff(both.diff.gene.ids, both.upreg.gene.ids)
intersect(mastersetgeneIDs, both.upreg.gene.ids) %>% length
intersect(mastersetgeneIDs, downreg.both.gene.ids) %>% length

tgfb.and.ra.up <- intersect(tgfb.upreg.gene.ids, intersect(ra.upreg.gene.ids, both.upreg.gene.ids))
tgfb.and.both.up <- setdiff(intersect(tgfb.upreg.gene.ids, both.upreg.gene.ids), ra.upreg.gene.ids)
ra.and.both.up <- setdiff(intersect(ra.upreg.gene.ids, both.upreg.gene.ids), tgfb.upreg.gene.ids)
both.up.only <- setdiff(both.upreg.gene.ids, union_all(tgfb.and.ra.up, tgfb.and.both.up, ra.and.both.up))

num.genes.in.upreg.venn.diagram <- length(tgfb.and.ra.up) + length(tgfb.and.both.up) + length(ra.and.both.up) + length(both.up.only)
num.genes.in.master.set <- length(mastersetgeneIDs)
percent.genes.removed <- 100 * (1 - (num.genes.in.master.set)/num.genes.in.upreg.venn.diagram)
print(sprintf("percent genes removed = %.2f", percent.genes.removed))

# #number in each category removed
# setdiff(tgfb.and.ra.up, mastersetgeneIDs) %>% length
# setdiff(tgfb.and.both.up, mastersetgeneIDs) %>% length
# setdiff(ra.and.both.up, mastersetgeneIDs) %>% length
# setdiff(both.up.only, mastersetgeneIDs) %>% length

######### uncomment this section to see how different definitions of upregulated genes affects the final set of upregulated genes
# # filter out diff peaks for which to do signal integration analysis on
# siUpregGenes1 <- geneTibAnno %>% filter(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`,
#                                         `TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`,
#                                         `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
#                                         `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
#                                         `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0)
# 
# # this definition is "upregulated by both individually OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals
# siUpregGenes2 <- geneTibAnno %>% filter(`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene` | ((`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
#                                         `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
#                                         `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
#                                         `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0)
# 
# # this definition is "upregulated by both individually OR upregulated in the both treatment with positive effects (that may not be "significant") from individual signals AND not upregulated in a single individual condition
# siUpregGenes3 <- geneTibAnno %>% filter(`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene` | ((`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
#                                         !((`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`) &  (`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) & !(`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
#                                         !((`TGFb-and-RA-low_isDeGene` | `TGFb-and-RA-med_isDeGene` | `TGFb-and-RA-high_isDeGene`) & !(`RA-low_isDeGene` | `RA-med_isDeGene` | `RA-high_isDeGene`) &  (`TGFb-low_isDeGene` | `TGFb-med_isDeGene` | `TGFb-high_isDeGene`)),
#                                         `RA-low_log2fc` > 0, `TGFb-low_log2fc` > 0,
#                                         `RA-med_log2fc` > 0, `TGFb-med_log2fc` > 0,
#                                         `RA-high_log2fc` > 0, `TGFb-high_log2fc` > 0,
#                                         `TGFb-and-RA-low_log2fc` > 0, `TGFb-and-RA-med_log2fc` > 0, `TGFb-and-RA-high_log2fc` > 0)
# table(siUpregGenes1$`integrationCategory-med-dose`)
# table(siUpregGenes2$`integrationCategory-med-dose`)
# table(siUpregGenes3$`integrationCategory-med-dose`)
# 
# genes.of.interest <- setdiff(siUpregGenes2$gene_name, siUpregGenes3$gene_name)
# tib.genes.of.interest <- siUpregGenes2 %>% filter(gene_name %in% genes.of.interest)
# new.supermultiplicative.geneIDs <- tib.genes.of.interest %>% filter(`integrationCategory-med-dose` == "super-multiplicative") %>% pull("gene_name")
# table(tib.genes.of.interest$`integrationCategory-med-dose`)
