library(tidyverse)
library(here)
library(grid)
library(gridExtra)

deseq.tib   <- read_tsv(here("extractedData", "DeSeqOutputAllConds.tsv"))
confint.tib <- read_tsv(here("extractedData", "geneExprConfIntsAndAddMultPredictions.tsv"))
big.tib <- left_join(confint.tib, deseq.tib, by = "ensg")

################################### typical direction -- effects of adding signals
p.both.up <- big.tib %>%
  # filter(`TGFb-and-RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_isDeGene` == 1, `RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_log2fc` > 0, `RA-med_log2fc` > 0) %>%
  filter(`EtOH-nlDensity_avgTPM` > 2) %>%
  # mutate(fcdiffrankorder = rank(`multPred-med-dose`/`EtOH-nlDensity_avgTPM`)) %>%
  mutate(fcdiffrankorder = rank(`fcDiffMultVsAddPred-med-dose`)) %>%
  ggplot() +
    geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`multPred-med-dose` / `EtOH-nlDensity_avgTPM`)), color = "blue") +
    geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`addPred-med-dose` / `EtOH-nlDensity_avgTPM`)), color = "red") +
    geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`TGFb-and-RA-med_avgTPM` / `EtOH-nlDensity_avgTPM`)), color = "black") +
    geom_errorbar(mapping = aes(x = fcdiffrankorder, 
                                ymax = log2(`TGFb-and-RA-med_upperBoundary` / `EtOH-nlDensity_avgTPM`),
                                ymin = log2(`TGFb-and-RA-med_lowerBoundary` / `EtOH-nlDensity_avgTPM`)), color = "gray") +
    xlab("rank order of ab/x^2, both signals up") +
    ylab("log-fold change") +
    theme_classic(base_size = 18)

p.both.down <- big.tib %>%
  # filter(`TGFb-and-RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_isDeGene` == 1, `RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_log2fc` < 0, `RA-med_log2fc` < 0) %>%
  filter(`EtOH-nlDensity_avgTPM` > 2) %>%
  # mutate(fcdiffrankorder = rank(`multPred-med-dose`/`EtOH-nlDensity_avgTPM`)) %>%
  mutate(fcdiffrankorder = rank(`fcDiffMultVsAddPred-med-dose`)) %>%
  ggplot() +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`multPred-med-dose` / `EtOH-nlDensity_avgTPM`)), color = "blue") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`addPred-med-dose` / `EtOH-nlDensity_avgTPM`)), color = "red") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`TGFb-and-RA-med_avgTPM` / `EtOH-nlDensity_avgTPM`)), color = "black") +
  geom_errorbar(mapping = aes(x = fcdiffrankorder, 
                              ymax = log2(`TGFb-and-RA-med_upperBoundary` / `EtOH-nlDensity_avgTPM`),
                              ymin = log2(`TGFb-and-RA-med_lowerBoundary` / `EtOH-nlDensity_avgTPM`)), color = "gray") +
  xlab("rank order of ab/x^2, both signals down") +
  ylab("log-fold change") +
  theme_classic(base_size = 18)

p.opposing <- big.tib %>%
  # filter(`TGFb-and-RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_isDeGene` == 1, `RA-med_isDeGene` == 1) %>%
  filter(((`TGFb-med_log2fc` > 0) & (`RA-med_log2fc` < 0)) | ((`TGFb-med_log2fc` < 0) & (`RA-med_log2fc` > 0))) %>%
  filter(`EtOH-nlDensity_avgTPM` > 2) %>%
  # mutate(fcdiffrankorder = rank(`multPred-med-dose`/`EtOH-nlDensity_avgTPM`)) %>%
  mutate(fcdiffrankorder = rank(`fcDiffMultVsAddPred-med-dose`)) %>%
  ggplot() +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`multPred-med-dose` / `EtOH-nlDensity_avgTPM`)), color = "blue") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`addPred-med-dose` / `EtOH-nlDensity_avgTPM`)), color = "red") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`TGFb-and-RA-med_avgTPM` / `EtOH-nlDensity_avgTPM`)), color = "black") +
  geom_errorbar(mapping = aes(x = fcdiffrankorder, 
                              ymax = log2(`TGFb-and-RA-med_upperBoundary` / `EtOH-nlDensity_avgTPM`),
                              ymin = log2(`TGFb-and-RA-med_lowerBoundary` / `EtOH-nlDensity_avgTPM`)), color = "gray") +
  xlab("rank order of ab/x^2, opposing signals") +
  ylab("log-fold change") +
  theme_classic(base_size = 18)

grid.arrange(p.both.up, p.both.down, p.opposing, ncol = 1)

################################### reverse direction -- effects of removing signals
# both "up" here means going down from the original prediction
p.both.up.reverse <- big.tib %>%
  # filter(`TGFb-and-RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_isDeGene` == 1, `RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_log2fc` > 0, `RA-med_log2fc` > 0) %>%
  filter(`EtOH-nlDensity_avgTPM` > 2) %>%
  # mutate(fcdiffrankorder = rank(`multPred-med-dose`/`EtOH-nlDensity_avgTPM`)) %>%
  mutate(fcdiffrankorder = rank(`fcDiffMultVsAddPredBackwards-med-dose`)) %>%
  ggplot() +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`multPredBackwards-med-dose` / `TGFb-and-RA-med_avgTPM`)), color = "blue") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`addPredBackwards-med-dose` / `TGFb-and-RA-med_avgTPM`)), color = "red") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`EtOH-nlDensity_avgTPM` / `TGFb-and-RA-med_avgTPM`)), color = "black") +
  geom_errorbar(mapping = aes(x = fcdiffrankorder, 
                              ymax = log2(`EtOH-nlDensity_upperBoundary` / `TGFb-and-RA-med_avgTPM`),
                              ymin = log2(`EtOH-nlDensity_lowerBoundary` / `TGFb-and-RA-med_avgTPM`)), color = "gray") +
  xlab("rank order of ab/x^2, reverse of both signals up") +
  ylab("log-fold change") +
  theme_classic(base_size = 18)


p.both.down.reverse <- big.tib %>%
  # filter(`TGFb-and-RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_isDeGene` == 1, `RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_log2fc` < 0, `RA-med_log2fc` < 0) %>%
  filter(`EtOH-nlDensity_avgTPM` > 2) %>%
  # mutate(fcdiffrankorder = rank(`multPred-med-dose`/`EtOH-nlDensity_avgTPM`)) %>%
  mutate(fcdiffrankorder = rank(`fcDiffMultVsAddPredBackwards-med-dose`)) %>%
  ggplot() +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`multPredBackwards-med-dose` / `TGFb-and-RA-med_avgTPM`)), color = "blue") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`addPredBackwards-med-dose` / `TGFb-and-RA-med_avgTPM`)), color = "red") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`EtOH-nlDensity_avgTPM` / `TGFb-and-RA-med_avgTPM`)), color = "black") +
  geom_errorbar(mapping = aes(x = fcdiffrankorder, 
                              ymax = log2(`EtOH-nlDensity_upperBoundary` / `TGFb-and-RA-med_avgTPM`),
                              ymin = log2(`EtOH-nlDensity_lowerBoundary` / `TGFb-and-RA-med_avgTPM`)), color = "gray") +
  xlab("rank order of ab/x^2, reverse of both signals down") +
  ylab("log-fold change") +
  theme_classic(base_size = 18)


p.opposing.reverse <- big.tib %>%
  # filter(`TGFb-and-RA-med_isDeGene` == 1) %>%
  filter(`TGFb-med_isDeGene` == 1, `RA-med_isDeGene` == 1) %>%
  filter(((`TGFb-med_log2fc` > 0) & (`RA-med_log2fc` < 0)) | ((`TGFb-med_log2fc` < 0) & (`RA-med_log2fc` > 0))) %>%
  filter(`EtOH-nlDensity_avgTPM` > 2) %>%
  # mutate(fcdiffrankorder = rank(`multPred-med-dose`/`EtOH-nlDensity_avgTPM`)) %>%
  mutate(fcdiffrankorder = rank(`fcDiffMultVsAddPredBackwards-med-dose`)) %>%
  ggplot() +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`multPredBackwards-med-dose` / `TGFb-and-RA-med_avgTPM`)), color = "blue") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`addPredBackwards-med-dose` / `TGFb-and-RA-med_avgTPM`)), color = "red") +
  geom_point(mapping = aes(x = fcdiffrankorder, y = log2(`EtOH-nlDensity_avgTPM` / `TGFb-and-RA-med_avgTPM`)), color = "black") +
  geom_errorbar(mapping = aes(x = fcdiffrankorder, 
                              ymax = log2(`EtOH-nlDensity_upperBoundary` / `TGFb-and-RA-med_avgTPM`),
                              ymin = log2(`EtOH-nlDensity_lowerBoundary` / `TGFb-and-RA-med_avgTPM`)), color = "gray") +
  xlab("rank order of ab/x^2, reverse of opposing signals") +
  ylab("log-fold change") +
  theme_classic(base_size = 18)

grid.arrange(p.both.up.reverse, p.both.down.reverse, p.opposing.reverse, ncol = 1)

grid.arrange(p.both.up, p.both.up.reverse, p.both.down, p.both.down.reverse,
             p.opposing, p.opposing.reverse, ncol = 2, nrow = 3)

grid.arrange(p.both.up, p.both.down.reverse, p.both.down, p.both.up.reverse,
             p.opposing, p.opposing.reverse, ncol = 2, nrow = 3)

grid.arrange(p.both.down, p.both.down.reverse, ncol = 1)
grid.arrange(p.both.up, p.both.up.reverse, ncol = 1)
