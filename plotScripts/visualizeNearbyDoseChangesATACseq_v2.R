library(tidyverse)
library(GenomicRanges)
library(here)
library(gtable)
library(grid)

#########
samplemetadata   <- read_tsv(here('sampleMetadata_SI2-SI4.txt'))

TSS.window.radius <- 100000  #  include peaks within this many base pairs of a gene's transcription start site
  
diffpeakstib <- read_tsv(here('extractedData', 'differentialAtacPeaks.tsv'))
deseqtib     <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.tsv'))

## to do-- make this table "wide" and calc average values for each condition in order to plot nearby gene expr changes
rnaseqPipelineOut <- filter(read_tsv(here('extractedData', 'si2-si4_RNA-seq-pipeline-output-normalized.tsv')), sampleID == "01-TGFb-low")
TSS.loc.tib  <- read_tsv(here('refs', 'EnsgToTssMapping.tsv'))
TSS.loc.tib.extended <- inner_join(rnaseqPipelineOut, TSS.loc.tib, by = "gene_id") %>% mutate(ensg=gene_id) %>% dplyr::select(-gene_id)

genes_dose_responsive_up <- filter(deseqtib, 
                   `RA-med_isDeGene` == 1, 
                   `RA-low_log2fc` > 0, 
                   (`RA-low_avgTPM`) * 1.75 < (`RA-high_avgTPM`))$ensg
genes_dose_responsive_down <- filter(deseqtib, 
                   `RA-med_isDeGene` == 1, 
                   `RA-low_log2fc` < 0, 
                   (`RA-low_avgTPM`) > (`RA-high_avgTPM`) * 1.75)$ensg
gene_set <- genes_dose_responsive_up

outputfolder <- here("plots", "peaksNearGenesTrackView", "dose_responsive_up")
# outputfolder <- here("plots", "peaksNearGenesTrackView", "dose_responsive_down")
#########



# add TSS coordinate column, which depends on strand, to TSS location table
TSS.loc.tib.extended[["tss"]] <- rep(0, nrow(TSS.loc.tib))
TSS.loc.tib.extended[["tss"]][which(TSS.loc.tib.extended$strand == "+")] <- TSS.loc.tib.extended[["tx_start"]][which(TSS.loc.tib.extended$strand == "+")]
TSS.loc.tib.extended[["tss"]][which(TSS.loc.tib.extended$strand == "-")] <- TSS.loc.tib.extended[["tx_end"]][which(TSS.loc.tib.extended$strand == "-")]


# gene.to.plot <- sample(deseqtib$ensg, 1)
# gene.to.plot <- "ENSG00000087916" # this gene causes an error in the script

makeBeeswarmPlot <- function(geneToPlot, tibforplot, deseq.tib, ordertib) {
  
  gene.symbol <- deseq.tib$hgnc.symbol[deseq.tib$ensg == geneToPlot]
  
  tfp <- tibforplot %>% filter(ensg == geneToPlot, 
                               condition != "TGFb-low", 
                               condition != "TGFb-med",
                               condition != "TGFb-high")
  
  p <- ggplot() +
    geom_point(tfp, mapping = aes(x = order, y = yaxis, color = condition)) + 
    # scale_color_grey(start=0.65, end = 0.05) + 
    scale_color_manual(c("EtOH-halfDensity", "EtOH-nlDensity", "EtOH-highDensity",
                         "RA-low", "RA-med", "RA-high",
                         "TGFb-and-RA-low", "TGFb-and-RA-med", "TGFb-and-RA-high"), 
                       values = c("black", "black", "black",
                                  "red", "red", "red",
                                  "purple", "purple", "purple"), name = NULL) + 
    scale_x_continuous(breaks=seq(1,9,1), minor_breaks = NULL, labels = ordertib$condition) +
    xlab("") +
    ylab("") + expand_limits(y = 0) +
    theme_classic(base_size = 6) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          legend.position = "none",
          axis.line.x=element_blank(),
          axis.line.y=element_blank()) 
  
  return(p)
}

makeTrackOfBeeswarmPlots <- function(nearbyGeneInfo, deseqtib, yaxis = "tpm") {
  counter <- 1
  savedplots <- list()
  plotcenters <- c()
  plotymaxs <- c()
  plotymins <- c()
  for (i in 1:nrow(nearbyGeneInfo)) {
    ensg <- nearbyGeneInfo$ensg[i]
    tssloc <- nearbyGeneInfo$tss[i]
    if (ensg %in% deseqtib$ensg) {
      if (yaxis == "tpm") {
        # yaxis.only.subset <- deseqtib %>% dplyr::select(matches("^(TGFb|RA)-(low|med|high)_avgTPM$"))
        # yaxis.only.subset[["EtOH-nlDensity_avgTPM"]] <- deseqtib[["EtOH-nlDensity_avgTPM"]]
        ind.yaxis.values <- dplyr::select(deseqtib, matches("(ensg)|(_tpm)") ) %>% dplyr::select(-"46-EtOH-nlDensity_tpm")
      } else if (yaxis == "log2fc") {
        # yaxis.only.subset <- deseqtib %>% dplyr::select(matches("^(TGFb|RA)-(low|med|high)_avgTPM$"))
        # yaxis.only.subset[["EtOH-nlDensity_avgTPM"]] <- deseqtib[["EtOH-nlDensity_avgTPM"]]
        ind.yaxis.tpm.values <- dplyr::select(deseqtib, matches("(ensg)|(_tpm)") ) %>% dplyr::select(-"46-EtOH-nlDensity_tpm")
        # this is going to be a messy block of code but should work. the column names will say "tpm" but the tpm will later be removed.
        # this will change the tpm value to a fold-change value in three stepwise operations, one for each replicate.
        # we will "mutate" the data frame to replace tpm values with foldchange values, then change all foldchange values to log fold change values.
        # the index numbers are manually taken from the samplemetadata spreadsheet, with +1 added to account for the offset due to the ensg column.
        indicesrep1 <- 1:12 + 1
        indicesrep2 <- c(13:22, 35:36) + 1
        indicesrep3 <- 23:34 + 1
        ind.yaxis.foldchange.values <- ind.yaxis.tpm.values
        ind.yaxis.foldchange.values[, indicesrep1] <- ind.yaxis.tpm.values[, indicesrep1] / pull(ind.yaxis.tpm.values, "05-EtOH-nlDensity_tpm")
        ind.yaxis.foldchange.values[, indicesrep2] <- ind.yaxis.tpm.values[, indicesrep2] / pull(ind.yaxis.tpm.values, "15-EtOH-nlDensity_tpm")
        ind.yaxis.foldchange.values[, indicesrep3] <- ind.yaxis.tpm.values[, indicesrep3] / pull(ind.yaxis.tpm.values, "27-EtOH-nlDensity_tpm")
        ind.yaxis.log2fc.values <- ind.yaxis.foldchange.values
        ind.yaxis.log2fc.values[, 1:36 + 1] <- log2(ind.yaxis.foldchange.values[, 1:36 + 1])
        ind.yaxis.values <- ind.yaxis.log2fc.values
      } else {
        stop("yaxis value in makeTrackOfBeeswarmPlots must be 'tpm' or 'log2fc'")
      }
      ind.yaxis.values.long <- gather(ind.yaxis.values, key = "sampleID", value = "yaxis", -"ensg") 
      ind.yaxis.values.long[['sampleID']] <- sapply(ind.yaxis.values.long$sampleID, function(x) strsplit(x, "_")[[1]][1])  # remove _tpm or _log2fc 
      
      tibforplot <- left_join(ind.yaxis.values.long, dplyr::select(samplemetadata, sampleID, condition), by="sampleID")
      ordertib <- tibble(condition = c("EtOH-halfDensity", "EtOH-nlDensity", "EtOH-highDensity",
                                       "RA-low", "RA-med", "RA-high", 
                                       "TGFb-and-RA-low", "TGFb-and-RA-med", "TGFb-and-RA-high"), 
                         order = 1:9)
      tibforplot <- left_join(tibforplot, ordertib, by="condition")
      tibforplot <- left_join(tibforplot, dplyr::select(samplemetadata, sampleID, replicate), by = "sampleID")
      replicate.spacing <- 0.15
      tibforplot[tibforplot[["replicate"]] == "rep1", "order"] <- tibforplot[tibforplot[["replicate"]] == "rep1", "order"] - replicate.spacing
      tibforplot[tibforplot[["replicate"]] == "rep3", "order"] <- tibforplot[tibforplot[["replicate"]] == "rep3", "order"] + replicate.spacing
      
      p.this <- makeBeeswarmPlot(ensg, tibforplot, deseqtib, ordertib)
      savedplots[[counter]] <- p.this
      
      plotcenters <- c(plotcenters, tssloc - this.tsstib$tss[1])
      builtplot <- ggplot_build(p.this)
      plotymaxs <- c(plotymaxs, builtplot$layout$panel_scales_y[[1]]$range$range[2])
      plotymins <- c(plotymins, builtplot$layout$panel_scales_y[[1]]$range$range[1])
      
      counter <- counter + 1
    }
  }
  max.yaxis.val <- max(plotymaxs)
  min.yaxis.val <- min(plotymins)
  for (i in 1:length(savedplots)) {
    savedplots[[i]] <- savedplots[[i]] + ylim(min.yaxis.val, max.yaxis.val)
  }
  
  # now mush these plots together as "annotations"
  yaxis_radius_p1 <- max.yaxis.val
  p1 <- ggplot(NULL) + 
    xlim(-1 * TSS.window.radius, TSS.window.radius) + 
    ylim(0, yaxis_radius_p1) +
    theme_classic() +
    ylab(paste0("rna-seq ", yaxis)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank()) 
  
  width_x_topplots <- 20000
  left_offset <- width_x_topplots * 0.085
  for (i in 1:length(savedplots)) {
    if (plotcenters[i] == 0) {
      plotmelaterindex <- i # skip the "central gene" to plot last so that it doesn't get plotted over
    } else {
      g <- ggplotGrob(savedplots[[i]])
      p1 <- p1 + annotation_custom(grob = g,
                                   xmin = plotcenters[i] - width_x_topplots/2 - left_offset,
                                   xmax = plotcenters[i] + width_x_topplots/2 - left_offset,
                                   ymin = 0,
                                   ymax = yaxis_radius_p1)
    }
  }
  # finally, plot the "central gene"
  g <- ggplotGrob(savedplots[[plotmelaterindex]])
  

  p1 <- p1 + annotation_custom(grob = g,
                               xmin = plotcenters[plotmelaterindex] - width_x_topplots/2 - left_offset,
                               xmax = plotcenters[plotmelaterindex] + width_x_topplots/2 - left_offset,
                               ymin = 0,
                               ymax = yaxis_radius_p1)
  return(p1)
} 


for (gene.to.plot in gene_set) {
  gene.to.plot <- "ENSG00000204264"
  gene.to.plot <- "ENSG00000240065"
  this.tsstib <- filter(TSS.loc.tib.extended, ensg == gene.to.plot)
  if (nrow(this.tsstib) == 0) {
    print(paste0("no TSS info available for ", gene.to.plot))
    next
  }
  this.hgncsymbol <- this.tsstib$hgnc_symbol[1]
  grTssWindowThisGene <- GRanges(seqnames = c(this.tsstib$chrom),
                                 ranges   = IRanges(start = this.tsstib$tss - TSS.window.radius,
                                                    end   = this.tsstib$tss + TSS.window.radius),
                                 strand   = this.tsstib$strand,
                                 gene     = this.tsstib$ensg)
  
  grTssAllGenes <- GRanges(seqnames = c(TSS.loc.tib.extended$chrom),
                           ranges   = IRanges(start = TSS.loc.tib.extended$tss,
                                              end   = TSS.loc.tib.extended$tss),
                           strand   = TSS.loc.tib.extended$strand,
                           gene     = TSS.loc.tib.extended$ensg)
  
  nearbyGeneInfo <- TSS.loc.tib.extended[subjectHits(findOverlaps(grTssWindowThisGene, grTssAllGenes, ignore.strand=TRUE)), ]
  nearbyGeneInfo$hgnc_symbol[is.na(nearbyGeneInfo$hgnc_symbol)] <- nearbyGeneInfo$ensg[is.na(nearbyGeneInfo$hgnc_symbol)]
   
  grDiffPeaks   <-  GRanges(seqnames = diffpeakstib$chrom,
                            ranges = IRanges(start = diffpeakstib$startLocs,
                                             end   = diffpeakstib$endLocs) ,
                            fragCtsControl = diffpeakstib$`EtOH-nlDensity-avgNormFragmentCounts`)
  
  TssWindowPeakOverlaps <- findOverlaps(grTssWindowThisGene, grDiffPeaks, ignore.strand=TRUE)
  nearbyPeakInfo <- diffpeakstib[subjectHits(TssWindowPeakOverlaps), ] %>% 
    mutate(peakmidpoint = (startLocs + endLocs) / 2) %>%
    mutate(relativepeakmidpoint = peakmidpoint - this.tsstib$tss[1])
  
  # Plots 0 and 1: show the gene expression changes associated to the central and neighboring genes, in tpm + log2fc 
  p0 <- makeTrackOfBeeswarmPlots(nearbyGeneInfo, deseqtib, yaxis = "log2fc")
  p1 <- makeTrackOfBeeswarmPlots(nearbyGeneInfo, deseqtib, yaxis = "tpm")
  grid.arrange(p0, p1, ncol = 1)

  # Plot 2: show the relative locations of differential peaks and nearby genes in the chosen radius from the central gene's TSS
  p2tib <- tibble(x = nearbyPeakInfo$relativepeakmidpoint, y = rep(0, nrow(nearbyPeakInfo)))
  p2 <- ggplot(p2tib, aes(x = x, y = y)) + 
    geom_line(data = tibble(x = c(-1 * TSS.window.radius, TSS.window.radius), y = c(0, 0)), mapping = aes(x=x, y=y)) +
    geom_linerange(ymin = -1, ymax = 1) +
    xlim(-1 * TSS.window.radius, TSS.window.radius) + ylim(-10, 10) + 
    geom_text(data = nearbyGeneInfo, mapping = aes(x = nearbyGeneInfo$tss - this.tsstib$tss[1], 
                                                   y = rep(-2, nrow(nearbyGeneInfo)), 
                                                   label = nearbyGeneInfo$hgnc_symbol)) +
    geom_point(data = nearbyGeneInfo, mapping = aes(x = nearbyGeneInfo$tss - this.tsstib$tss[1], 
                                                   y = rep(0, nrow(nearbyGeneInfo))), size = I(2))
  p2 <- p2 + theme_void()
  
  # Plot 3: show how the accessibility of each peak changes according to dosage or cell density
  # just use a for loop and keep adding new plots to the same graphics object
  p3 <- ggplot(data = NULL)
  for (condition.set in list(c("RA-low", "RA-med", "RA-high", "red"),
                             c("TGFb-and-RA-low", "TGFb-and-RA-med", "TGFb-and-RA-high", "purple"),
                             c("EtOH-halfDensity", "EtOH-nlDensity", "EtOH-highDensity", "black"))) {
    for (i in 1:nrow(nearbyPeakInfo)) {
      y1 <- pull(nearbyPeakInfo, paste0(condition.set[1], "-avgNormFragmentCounts"))[i]
      y2 <- pull(nearbyPeakInfo, paste0(condition.set[2], "-avgNormFragmentCounts"))[i]
      y3 <- pull(nearbyPeakInfo, paste0(condition.set[3], "-avgNormFragmentCounts"))[i]
      condcolor <- condition.set[4]
      xbase <- pull(nearbyPeakInfo, "relativepeakmidpoint")[i]
      xbaseincrement <- 1000
      thislineplottib <- tibble(x = c(xbase - xbaseincrement, xbase, xbase + xbaseincrement),
                                y = c(y1, y2, y3))
      p3 <- p3 + geom_line(data = thislineplottib, mapping = aes(x = x, y = y), color = condcolor)
    }
  }
  p3 <- p3 + 
    ylab("avg. nl. fragment counts (atac)") +
    xlim(-1 * TSS.window.radius, TSS.window.radius) +
    theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank())
  
  # Plot 4: add histograms to the bottom to display the distribution of change in accessibility: fold change, log2fc, and raw change
  condForNearbyPeaks <- 'RA'
  nearbyPeakInfoEtOHAvgNormFragCts <- pull(nearbyPeakInfo, "EtOH-nlDensity-avgNormFragmentCounts")
  nearbyPeakInfoAvgNormFragCts <- select(nearbyPeakInfo, matches(paste0('^', condForNearbyPeaks, '-(low|med|high)', '-avgNormFragmentCounts')))
  nearbyPeakInfoAvgFoldChange  <- select(nearbyPeakInfo, matches(paste0('^', condForNearbyPeaks, '-(low|med|high)', '-avgFoldchange')))
  nearbyPeakInfoIsDiffPeakThisCond  <- select(nearbyPeakInfo, matches(paste0('^', condForNearbyPeaks, '-(low|med|high)', '-isDiffPeak')))
  isDiffPeakAnyDose <- c()

  for (i in 1:nrow(nearbyPeakInfo)) {
    isDiffPeakAnyDose <- c(isDiffPeakAnyDose, any(as_vector(nearbyPeakInfoIsDiffPeakThisCond[i,])))
  }
  if (nrow(nearbyPeakInfo) == 0) {
    print("no nearby peaks in any condition for")
    print(gene.to.plot)
    next()
  } else if (sum(isDiffPeakAnyDose) == 0) {
    print("no nearby peaks for")
    print(gene.to.plot)
    next()
  }
  nearbyPeakInfoAvgNormFragCts     <- nearbyPeakInfoAvgNormFragCts[isDiffPeakAnyDose,]
  nearbyPeakInfoEtOHAvgNormFragCts <- nearbyPeakInfoEtOHAvgNormFragCts[isDiffPeakAnyDose]
  nearbyPeakInfoRawAccDiff       <- nearbyPeakInfoAvgNormFragCts - nearbyPeakInfoEtOHAvgNormFragCts
  nearbyPeakInfoAvgFoldChange    <- nearbyPeakInfoAvgFoldChange[isDiffPeakAnyDose,]
  nearbyPeakInfoAvgLogFoldChange <- log2(nearbyPeakInfoAvgFoldChange)
  
  p4 <- nearbyPeakInfoRawAccDiff %>% 
    mutate(local_peak_id = 1:nrow(nearbyPeakInfoRawAccDiff)) %>%
    gather(key = "dose", value = "AvgNormFragCountDiff", -local_peak_id) %>%
    mutate(dose = sapply(dose, function(x) paste0(condForNearbyPeaks, " ", strsplit(x, "-")[[1]][2], " dose"))) %>%
    mutate(dose = factor(dose, levels = paste0(condForNearbyPeaks, " ", c("low", "med", "high"), " dose"))) %>%
    ggplot(aes(AvgNormFragCountDiff)) + 
    geom_histogram() +
    facet_wrap(vars(dose), nrow=3) +
    theme_minimal(base_size = 10)
  xaxisfactor <- max(abs(min(nearbyPeakInfoRawAccDiff)), abs(max(nearbyPeakInfoRawAccDiff))) * 1.05
  p4 <- p4 + xlim(-1 * xaxisfactor, xaxisfactor)
  
  p5 <- nearbyPeakInfoAvgFoldChange %>% 
    mutate(local_peak_id = 1:nrow(nearbyPeakInfoAvgFoldChange)) %>%
    gather(key = "dose", value = "AvgFoldChangeDiff", -local_peak_id) %>%
    mutate(dose = sapply(dose, function(x) paste0(condForNearbyPeaks, " ", strsplit(x, "-")[[1]][2], " dose"))) %>%
    mutate(dose = factor(dose, levels = paste0(condForNearbyPeaks, " ", c("low", "med", "high"), " dose"))) %>%
    ggplot(aes(AvgFoldChangeDiff)) + 
    geom_histogram() +
    facet_wrap(vars(dose), nrow=3) +
    theme_minimal(base_size = 10) +
    theme(axis.title.y = element_blank())
  xaxisfactor <- max(abs(min(nearbyPeakInfoAvgFoldChange)), abs(max(nearbyPeakInfoAvgFoldChange))) * 1.05
  p5 <- p5 + xlim(0, xaxisfactor)
  
  p6 <- nearbyPeakInfoAvgLogFoldChange %>% 
    mutate(local_peak_id = 1:nrow(nearbyPeakInfoAvgLogFoldChange)) %>%
    gather(key = "dose", value = "AvgLog2FoldChangeDiff", -local_peak_id) %>%
    mutate(dose = sapply(dose, function(x) paste0(condForNearbyPeaks, " ", strsplit(x, "-")[[1]][2], " dose"))) %>%
    mutate(dose = factor(dose, levels = paste0(condForNearbyPeaks, " ", c("low", "med", "high"), " dose"))) %>%
    ggplot(aes(AvgLog2FoldChangeDiff)) + 
    geom_histogram() +
    facet_wrap(vars(dose), nrow=3) +
    theme_minimal(base_size = 10) +
    theme(axis.title.y = element_blank())
  xaxisfactor <- max(abs(min(nearbyPeakInfoAvgLogFoldChange)), abs(max(nearbyPeakInfoAvgLogFoldChange))) * 1.05
  p6 <- p6 + xlim(-1 * xaxisfactor, xaxisfactor)

  g0 <- ggplotGrob(p0)
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g3 <- ggplotGrob(p3)
  g4 <- ggplotGrob(p4)
  g5 <- ggplotGrob(p5)
  g6 <- ggplotGrob(p6)
  g7 <- cbind(g4, g6, g5)
  
  
  g <- rbind(g0, g1, g2, g3, size = "first")
  g$widths <- unit.pmax(g2$widths, g3$widths)
  # grid.newpage()
  # grid.draw(g)
  
  layout.matx <- rbind(c(1), 
                       c(1),
                       c(1),
                       c(2))
  track.panel.final <- grid.arrange(g, g7, layout_matrix = layout.matx)

  ggsave(paste0(outputfolder, '/', this.hgncsymbol, '_', gene.to.plot, ".svg"), plot = track.panel.final, width = 16, height = 9)
}