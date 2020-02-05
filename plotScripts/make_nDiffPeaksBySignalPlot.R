library(tidyverse)
library(here)


peakTibAnno   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))


ndiffpeaks.ralow <- sum(peakTibAnno$`RA-low-isDiffPeak`)
ndiffpeaks.ramed <- sum(peakTibAnno$`RA-med-isDiffPeak`)
ndiffpeaks.rahigh <- sum(peakTibAnno$`RA-high-isDiffPeak`)
ndiffpeaks.tgfblow <- sum(peakTibAnno$`TGFb-low-isDiffPeak`)
ndiffpeaks.tgfbmed <- sum(peakTibAnno$`TGFb-med-isDiffPeak`)
ndiffpeaks.tgfbhigh <- sum(peakTibAnno$`TGFb-high-isDiffPeak`)
ndiffpeaks.bothlow <- sum(peakTibAnno$`TGFb-and-RA-low-isDiffPeak`)
ndiffpeaks.bothmed <- sum(peakTibAnno$`TGFb-and-RA-med-isDiffPeak`)
ndiffpeaks.bothhigh <- sum(peakTibAnno$`TGFb-and-RA-high-isDiffPeak`)

ndiffpeaks.vec <- c(ndiffpeaks.ralow, ndiffpeaks.ramed, ndiffpeaks.rahigh, 
              ndiffpeaks.tgfblow, ndiffpeaks.tgfbmed, ndiffpeaks.tgfbhigh, 
              ndiffpeaks.bothlow, ndiffpeaks.bothmed, ndiffpeaks.bothhigh)
ord.factor <- c("RA, 50 nM", "RA, 200 nM", "RA, 400 nM", 
                "TGFb, 1.25 ng/mL", "TGFb, 5 ng/mL", "TGFb, 10ng/mL",
                "Both, 50 nM/1.25 ng/mL", "Both, 200 nM/5 ng/mL", "Both, 400 nM/10 ng/mL")
barplot.xlabel <- factor(ord.factor, levels = ord.factor)

bptib <- tibble(signal = barplot.xlabel, nDiffPeaks = ndiffpeaks.vec)
p1 <- ggplot(bptib, aes(x = barplot.xlabel, y = nDiffPeaks)) + 
  geom_bar(stat = "identity", width = .5, color = "black") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.4))
p1 <- p1 + xlab("") + ylab("")
p1
ggsave(here('plots', 'nDiffPeaksBySignalPlot.svg'), plot = p1, width = 4, height = 3, useDingbats=FALSE)

svg(filename=here('plots', 'nDiffPeaksBySignalPlot2.svg'), 
    width=4, 
    height=3, 
    pointsize=6)
p1
dev.off()