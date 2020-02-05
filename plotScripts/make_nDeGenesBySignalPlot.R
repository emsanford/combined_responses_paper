library(tidyverse)
library(here)


deseqTibAnno  <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.annotated.tsv'))
peakTibAnno   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))


ndeg.ralow <- sum(deseqTibAnno$`RA-low_isDeGene`)
ndeg.ramed <- sum(deseqTibAnno$`RA-med_isDeGene`)
ndeg.rahigh <- sum(deseqTibAnno$`RA-high_isDeGene`)
ndeg.tgfblow <- sum(deseqTibAnno$`TGFb-low_isDeGene`)
ndeg.tgfbmed <- sum(deseqTibAnno$`TGFb-med_isDeGene`)
ndeg.tgfbhigh <- sum(deseqTibAnno$`TGFb-high_isDeGene`)
ndeg.bothlow <- sum(deseqTibAnno$`TGFb-and-RA-low_isDeGene`)
ndeg.bothmed <- sum(deseqTibAnno$`TGFb-and-RA-med_isDeGene`)
ndeg.bothhigh <- sum(deseqTibAnno$`TGFb-and-RA-high_isDeGene`)

ndeg.vec <- c(ndeg.ralow, ndeg.ramed, ndeg.rahigh, 
              ndeg.tgfblow, ndeg.tgfbmed, ndeg.tgfbhigh, 
              ndeg.bothlow, ndeg.bothmed, ndeg.bothhigh)
ord.factor <- c("RA, 50 nM", "RA, 200 nM", "RA, 400 nM", 
                "TGFb, 1.25 ng/mL", "TGFb, 5 ng/mL", "TGFb, 10ng/mL",
                "Both, 50 nM/1.25 ng/mL", "Both, 200 nM/5 ng/mL", "Both, 400 nM/10 ng/mL")
barplot.xlabel <- factor(ord.factor, levels = ord.factor)

bptib <- tibble(signal = barplot.xlabel, nDeGenes = ndeg.vec)
p1 <- ggplot(bptib, aes(x = barplot.xlabel, y = nDeGenes)) + 
  geom_bar(stat = "identity", width = .5, color = "black") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4))
p1 <- p1 + xlab("") + ylab("")
p1
ggsave(here('plots', 'nDeGenesBySignalPlot.svg'), plot = p1, width = 4, height = 3, useDingbats=FALSE)

svg(filename=here('plots', 'nDeGenesBySignalPlot2.svg'), 
    width=4, 
    height=3, 
    pointsize=6)
p1
dev.off()