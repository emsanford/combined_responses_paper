library(tidyverse)
library(motifmatchr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
motif.pwms <- read_rds("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/mostVariableMotifs_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.rds")
#now do motif matching on this region to display RA & TGFb matches
# range2 <- chr1:12534290   12566459      #
# range1 <- chr3:189564184  189578559     #
locs.to.annotate <- GRanges(seqnames  = c("chr1", "chr3"),
                            ranges    = IRanges(start = c(12534290, 189564184),
                            end   = c(12566459, 189578559)))

motif.locs.res <- matchMotifs(motif.pwms, locs.to.annotate, genome = BSgenome.Hsapiens.UCSC.hg38, out = "positions")
rara.granges.res  <- motif.locs.res$ENSG00000131759_LINE2987_RARA_D_N3
smad3.granges.res <- motif.locs.res$ENSG00000166949_LINE3324_SMAD3_D_N1
smad4.granges.res <- motif.locs.res$ENSG00000141646_LINE3314_SMAD4_D
smad9.granges.res <- motif.locs.res$ENSG00000120693_LINE19375_SMAD9_I_N1

rara.chr <- as.vector(seqnames(rara.granges.res))
smad3.chr <- as.vector(seqnames(smad3.granges.res))
smad4.chr <- as.vector(seqnames(smad4.granges.res))
smad9.chr <- as.vector(seqnames(smad9.granges.res))

rara.starts <- start(rara.granges.res)
smad3.starts <- start(smad3.granges.res)
smad4.starts <- start(smad4.granges.res)
smad9.starts <- start(smad9.granges.res)

rara.ends  <- end(rara.granges.res)
smad3.ends <- end(smad3.granges.res)
smad4.ends <- end(smad4.granges.res)
smad9.ends <- end(smad9.granges.res)

motif.out.tib <- tibble(chr = c(rara.chr, smad3.chr, smad4.chr, smad9.chr),
                        start = c(rara.starts, smad3.starts, smad4.starts, smad9.starts),
                        end   = c(rara.ends, smad3.ends, smad4.ends, smad9.ends),
                        motif = c(rep('RARA', length(rara.chr)), rep('SMAD3', length(smad3.chr)), rep('SMAD4', length(smad4.chr)), rep('SMAD9', length(smad9.chr))))

write_tsv(motif.out.tib %>% filter(motif == "RARA"),  here("extractedData", "motif_tracks", "rara_motifs_fig2_intervals.bed"),  col_names = F)
write_tsv(motif.out.tib %>% filter(motif == "SMAD3"), here("extractedData", "motif_tracks", "smad3_motifs_fig2_intervals.bed"), col_names = F)
write_tsv(motif.out.tib %>% filter(motif == "SMAD4"), here("extractedData", "motif_tracks", "smad4_motifs_fig2_intervals.bed"), col_names = F)
write_tsv(motif.out.tib %>% filter(motif == "SMAD9"), here("extractedData", "motif_tracks", "smad9_motifs_fig2_intervals.bed"), col_names = F)


# retrieve diff peaks in these locations
diffpeaks <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.tsv")

reg1.peaks <- diffpeaks %>%
  filter(chrom == 'chr1',
         startLocs >= 12534290,
         endLocs <= 12566459)

reg2.peaks <- diffpeaks %>%
  filter(chrom == 'chr3',
         startLocs >= 189564184,
         endLocs <= 189578559)
