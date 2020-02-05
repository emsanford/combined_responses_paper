library(tidyverse)
library(here)

geneTib        <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
samplemetadata <- read_tsv(here('sampleMetadata_SI2-SI4.txt'))
relevantmetadata <- filter(samplemetadata, as.integer(substr(sampleID,1,2)) %in% c(1:16, 19:36, 51, 52))  # discard RNA-seq samples that failed library prep (17, 18) or extra samples from second run (46-50) 


makeCompositeBeeswarmPlot <- function(geneToPlot, geneTib, tpmLongFormatTib) {
  filt.tib <- filter(geneTib, gene_name == geneToPlot)

  addpred_low  <- filt.tib[["addPred-low"]][1] 
  addpred_med  <- filt.tib[["addPred-med"]][1] 
  addpred_high <- filt.tib[["addPred-high"]][1] 
  multpred_low  <- filt.tib[["multPred-low"]][1] 
  multpred_med  <- filt.tib[["multPred-med"]][1] 
  multpred_high <- filt.tib[["multPred-high"]][1] 
  lineradius <- 0.35
  linewidth <- 1
  xposns_bothSignals <- c(10, 11, 12)

  tfp <- tpmLongFormatTib %>% filter(gene_name == geneToPlot)
  
  p <- ggplot()  +
    geom_point(tfp, mapping = aes(x = order, y = TPM, color = replicate)) + 
    # scale_color_grey(start=0.65, end = 0.05) + 
    scale_color_manual(c("mult pred", "add pred", "rep1", "rep2", "rep3"), 
                       values = c("#009292", "#DB6D00", "#9C9C9C", "#707070", "#000000"), name = NULL) + 
    scale_x_continuous(breaks=seq(1,12,1), minor_breaks = NULL, labels = ordertib$condition) +
    labs(title = geneToPlot) + 
    xlab("condition and dosage") +
    ylab("transcripts per million (TPM)") + expand_limits(y = 0) +
    theme_classic(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

  # add minimalist error bars, without "chartjunk caps"
  tfpSubsetForErrorBars <- tfp %>% filter(replicate == "rep2") %>% dplyr::select(condition, order)
  errorBarLineWidth <- 0.5
  errorBarColor <- "rep1"
  x_vec <- c()
  y_lci_vec <- c()
  y_uci_vec <- c()
  for (ii in 1:nrow(tfpSubsetForErrorBars) ) {
    condName  <- tfpSubsetForErrorBars$condition[ii]
    condOrder <- tfpSubsetForErrorBars$order[ii]
    lci <- filt.tib[[paste0(condName, "_pseudoConfintLowerBoundaryTPM")]][1]
    uci <- filt.tib[[paste0(condName, "_pseudoConfintUpperBoundaryTPM")]][1]
    x_vec <- c(x_vec, condOrder)
    y_lci_vec <- c(y_lci_vec, lci)
    y_uci_vec <- c(y_uci_vec, uci)
  }
  p <- p + geom_segment(aes(x = x_vec, 
                            y = y_lci_vec, 
                            xend = x_vec, 
                            yend = y_uci_vec, color = errorBarColor),  size = errorBarLineWidth) 
  
  # add the "notch" in the center of the error bars to depict the mean value
  mean_notch_lineWidth <- 0.5
  mean_notch_line_radius <- 0.35
  x_vec <- c()
  y_vec <- c()
  for (ii in 1:nrow(tfpSubsetForErrorBars) ) {
    condName  <- tfpSubsetForErrorBars$condition[ii]
    condOrder <- tfpSubsetForErrorBars$order[ii]
    avgTpm <- filt.tib[[paste0(condName, "_avgTPM")]][1]
    x_vec <- c(x_vec, condOrder)
    y_vec <- c(y_vec, avgTpm)
  }
  p <- p + geom_segment(aes(x = x_vec - mean_notch_line_radius, 
                            y = y_vec, 
                            xend = x_vec + mean_notch_line_radius, 
                            yend = y_vec, color = errorBarColor),  size = mean_notch_lineWidth) 
  
  # add additive and multiplicative predictions
  p <- p +
    geom_segment(mapping = aes(x = xposns_bothSignals - lineradius, 
                               y = c(addpred_low, addpred_med, addpred_high), color = "add pred",),  
                 xend = xposns_bothSignals + lineradius, 
                 yend = c(addpred_low, addpred_med, addpred_high), size = linewidth) +
    geom_segment(mapping = aes(x = xposns_bothSignals - lineradius, 
                               y = c(multpred_low, multpred_med, multpred_high), color = "mult pred"), 
                 xend = xposns_bothSignals + lineradius, 
                 yend = c(multpred_low, multpred_med, multpred_high), size = linewidth)
  
  return(p)
}

# this "data massaging" block of code helps create a tibble needed to make the beeswarm plot with the correct x and y positions and x axis ordering
ind.tpm.values <- dplyr::select(geneTib, matches("(gene_name)|(_tpm)") ) %>% dplyr::select(-"46-EtOH-nlDensity_tpm")
ind.tpm.values.long <- gather(ind.tpm.values, key = "sampleID", value = "TPM", -gene_name)
ind.tpm.values.long[['sampleID']] <- sapply(ind.tpm.values.long$sampleID, function(x) strsplit(x, "_")[[1]][1])

tibforplot <- left_join(ind.tpm.values.long, dplyr::select(relevantmetadata, sampleID, condition), by="sampleID")
ordertib <- tibble(condition = c("EtOH-halfDensity", "EtOH-nlDensity", "EtOH-highDensity",
                                 "RA-low", "RA-med", "RA-high", 
                                 "TGFb-low", "TGFb-med", "TGFb-high", 
                                 "TGFb-and-RA-low", "TGFb-and-RA-med", "TGFb-and-RA-high"), 
                   order = 1:12)
tibforplot <- left_join(tibforplot, ordertib, by="condition")
tibforplot <- left_join(tibforplot, dplyr::select(relevantmetadata, sampleID, replicate), by = "sampleID")
replicate.spacing <- 0.15
tibforplot[tibforplot[["replicate"]] == "rep1", "order"] <- tibforplot[tibforplot[["replicate"]] == "rep1", "order"] - replicate.spacing
tibforplot[tibforplot[["replicate"]] == "rep3", "order"] <- tibforplot[tibforplot[["replicate"]] == "rep3", "order"] + replicate.spacing

genes.to.plot <- c('EPHB2', 'MAP3K1', 'GPRC5A', 'RIPK4') 
for (gene.to.plot in genes.to.plot) {
  outputLoc <- here("plots", "beeSwarmPlotsForPaperFigs", paste0(gene.to.plot, "_beeswarm.svg"))
  p <- makeCompositeBeeswarmPlot(gene.to.plot, geneTib, tibforplot)
  # optional: remove things from the plots that we will edit in illustrator. (re-add them later in illustrator)
  # p <- p + theme(axis.text.y = element_blank()) + ylab("") + xlab("")
  ggsave(outputLoc, plot = p, width = 5.25, height = 4)
}

genes.to.plot <- genes_of_interest  # i cheat here and define genes_of_interest in another script in Rstudio to explore random categories of genes
for (gene.to.plot in genes.to.plot) {
  outputLoc <- here("plots", "beeswarm_playground", paste0(gene.to.plot, "_beeswarm.svg"))
  p <- makeCompositeBeeswarmPlot(gene.to.plot, geneTib, tibforplot)
  # optional: remove things from the plots that we will edit in illustrator. (re-add them later in illustrator)
  # p <- p + theme(axis.text.y = element_blank()) + ylab("") + xlab("")
  ggsave(outputLoc, plot = p, width = 5.25, height = 4)
}

