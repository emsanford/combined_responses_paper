library(tidyverse)
library(here)

peakTibAnno   <- read_tsv(here('extractedData', 'differentialAtacPeaks.annotated.tsv'))

# filter medium dose upregulated peaks
filt_tib <- filter(peakTibAnno, `RA-med-isDiffPeak`, `TGFb-med-isDiffPeak`,
                   `RA-med-avgFoldchange` > 1, `TGFb-med-avgFoldchange` > 1)

ra.fc   <- filt_tib$`RA-med-avgFoldchange`
tgfb.fc <- filt_tib$`TGFb-med-avgFoldchange`
both.fc <- filt_tib$`TGFb-and-RA-med-avgFoldchange`
addpred.fc  <- ra.fc + tgfb.fc - 1
multpred.fc <- ra.fc * tgfb.fc

plot_tib <- tibble(addpred = addpred.fc, multpred = multpred.fc, measuredfc = both.fc)

p1 <- ggplot(plot_tib, aes(x = addpred.fc, y = both.fc)) + 
  geom_point(size = I(0.5)) + 
  xlim(0, 30) + ylim(0,30) +
  geom_smooth(method = "lm", se = F) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("fold-change of additive prediction") +
  ylab("measured fold change") +
  theme_classic(base_size = 12) + coord_fixed()
p1

p2 <- ggplot(plot_tib, aes(x = multpred.fc, y = both.fc)) + 
  geom_point(size = I(0.5)) + 
  xlim(0, 30) + ylim(0,30) +
  geom_smooth(method = "lm", se = F) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("fold-change of multiplicative prediction") +
  ylab("measured fold change") + 
  theme_classic(base_size = 12) + coord_fixed()
p2

ggsave(here('plots','f30plot2.svg'), p1)
ggsave(here('plots','f30plot3.svg'), p2)

