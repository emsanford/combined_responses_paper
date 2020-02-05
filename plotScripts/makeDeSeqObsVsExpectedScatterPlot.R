library("tidyverse")
library("DESeq2")
library(GenomicRanges)
library(here)
library(grid)
library(gridExtra)

deSeqTibble <- read_tsv(here('extractedData', 'DeSeqOutputAllConds.tsv'))

tib.subset.OR.logic <- dplyr::select(deSeqTibble, ensg, `EtOH-nlDensity_avgCounts`, matches("log2fc"), matches("avgCounts"), matches("isDeGene")) %>%
  dplyr::select(ensg, `EtOH-nlDensity_avgCounts`, matches("TGFb-med_|RA-med")) %>%
  filter(`RA-med_isDeGene` == 1 | `TGFb-med_isDeGene` == 1, `EtOH-nlDensity_avgCounts` >= 5)

tib.subset.AND.logic <- dplyr::select(deSeqTibble, ensg, `EtOH-nlDensity_avgCounts`, matches("log2fc"), matches("avgCounts"), matches("isDeGene")) %>%
  dplyr::select(ensg, `EtOH-nlDensity_avgCounts`, matches("TGFb-med|RA-med")) %>%
  filter(`RA-med_isDeGene` == 1, `TGFb-med_isDeGene` == 1, `EtOH-nlDensity_avgCounts` >= 5)

# observed vs. expected, multiplicative model, log scale, differentially expressed under TGFb OR RA
axis.tick.locs <- seq(-6, 9, by = 3)
p1 <- ggplot(tib.subset.OR.logic, aes(x=(`TGFb-med_log2fc` + `RA-med_log2fc`), y=`TGFb-and-RA-med_log2fc`)) +
  geom_point(size=I(0.30)) +
  xlab("Predicted response: log-fold change in TGFb\nplus log-fold change in RA") +
  ylab("Measured response: log-fold change\nin combined TGFb + RA treatment") +
  coord_fixed() +
  scale_x_continuous(breaks = axis.tick.locs) +
  scale_y_continuous(breaks = axis.tick.locs) +
  theme_light(base_size = 20)
ggsave(here('plots', 'obsVsExpectedGeneExpr_logScale_multModel_OR_logic.svg'), plot = p1)
p1


# observed vs. expected, multiplicative model, log scale, differentially expressed under single TGFb and single RA tx
axis.tick.locs <- seq(-6, 9, by = 3)
p2 <- ggplot(tib.subset.AND.logic, aes(x=(`TGFb-med_log2fc` + `RA-med_log2fc`), y=`TGFb-and-RA-med_log2fc`)) +
  geom_point(size=I(0.6)) +
  xlab("Predicted response: log-fold change in TGFb\nplus log-fold change in RA") +
  ylab("Measured response: log-fold change\nin combined TGFb + RA treatment") +
  coord_fixed() +
  scale_x_continuous(breaks = axis.tick.locs) +
  scale_y_continuous(breaks = axis.tick.locs) +
  theme_light(base_size = 20)
ggsave(here('plots', 'obsVsExpectedGeneExpr_logScale_multModel_AND_logic.svg'), plot = p2)
p2

# and logic: show distribution of single-signal effects for these genes
p3 <- ggplot(tib.subset.AND.logic, aes(x=`TGFb-med_log2fc`, y=`RA-med_log2fc`)) +
  geom_point(size=I(0.6)) +
  xlab("log-fold change in TGFb") +
  ylab("log-fold change in RA") +
  ggtitle("genes differentially expressed in RA and TGFb") +
  coord_fixed() +
  theme_light(base_size = 12)
ggsave(here('plots', 'tgfbEffect_vs_raEffect_logScale_AND_logic_filters.svg'), plot = p3)
p3

# 
# 
# ##### below::: legacy
# 
# # plot for genes that are differentially expressed in at least one condition
# filtTib <- filter(deSeqTibble, isDeGeneRA3d | isDeGeneTGFb3d, deSeqTibble$avgTPM_EtOHnlDensity > 2.5)
# filtTib %>%
#   ggplot(aes(log2fcTGFb, log2fcRA)) + geom_point(size=I(0.75)) +
#   xlab('log-fold change in TGFb') + 
#   ylab('log-fold change in RA') + 
#  # ggtitle('The log-fold-change of differentially expressed \ngenes is weakly correlated between RA and TGFb treatments') + 
#   coord_fixed() + xlim(-4.5, 4.5) + ylim(-4.5, 4.5) +
#   theme(axis.text=element_text(size=16), axis.title=element_text(size=50), title = element_text(size=20)) +
#   theme_light(base_size=20)
# cor(filtTib$log2fcTGFb, filtTib$log2fcRA)
# 
# filtTib <- filter(deSeqTibble, isDeGeneRA3d | isDeGeneTGFb3d, deSeqTibble$avgTPM_EtOHnlDensity > 2.5)
# filtTib %>%
#   ggplot(aes(log2fcTGFb + log2fcRA, log2fcBoth)) + geom_point(size=I(0.75)) +
#   xlab('Multiplicative model prediction \n (log-fold change in TGFb + log-fold change in RA)') + 
#   ylab('measured log-fold change in TGFb + RA treatment') + 
#   #ggtitle('The log-fold-change of the integrated response genes \nis well-predicted by the individual signal responses') + 
#   coord_fixed() + xlim(-6, 6) + ylim(-6, 6) +
#   theme(axis.text=element_text(size=16), axis.title=element_text(size=20), title = element_text(size=20)) +
#   theme_light(base_size=17)
# cor(filtTib$log2fcTGFb + filtTib$log2fcRA, filtTib$log2fcBoth)
# 
# 
# filtTib <- deSeqTibble %>%
#   mutate(a = 2^log2fcRA - 1, b = 2^log2fcTGFb - 1, addPredFC = a + b + 1, multPredFC = 2^(log2fcRA + log2fcTGFb)) %>%
#   mutate(addMultPredDiffFC = multPredFC - addPredFC, addMultPredAbsDiffFC = abs(addMultPredDiffFC))
# 
# filtTib %>%
#   filter(isDeGeneRA3d & isDeGeneTGFb3d, deSeqTibble$avgTPM_EtOHnlDensity > 2.5) %>%
#   ggplot(aes(x=multPredFC, y=addPredFC)) + geom_point(size=I(0.75)) +
#   xlab('Multiplicative model prediction (fold change difference)') + 
#   ylab('Additive model prediction (fold change difference)') + 
#   #ggtitle('The log-fold-change of the integrated response genes \nis well-predicted by the individual signal responses') + 
#   coord_fixed() + xlim(-2, 18) + ylim(-2, 18) + geom_abline(slope=1, intercept=0) +
#   theme(axis.text=element_text(size=16), axis.title=element_text(size=20), title = element_text(size=20)) +
#   theme_light(base_size=17)
# 
# #calculate percent of genes where add and mult predictions are close to each other
# filtTib %>%
#   filter(isDeGeneRA3d & isDeGeneTGFb3d, deSeqTibble$avgTPM_EtOHnlDensity > 2.5) %>%
#   ggplot(aes(abs(multPredFC-addPredFC))) + geom_histogram(binwidth = 1) + xlim(0, 10) + ylim(0, 200) + theme_light(base_size=18) +
#   xlab('Difference in fold-change prediction between\nadditive and multiplicative models') + 
#   ylab('number of genes')
# tibForPctSimilarPred <- filtTib %>%
#   filter(isDeGeneRA3d & isDeGeneTGFb3d, deSeqTibble$avgTPM_EtOHnlDensity > 2.5)
# frac.within.one.fc <- sum(tibForPctSimilarPred$addMultPredAbsDiffFC <= 1) / length(tibForPctSimilarPred$addMultPredAbsDiffFC)
# frac.within.two.fc <- sum(tibForPctSimilarPred$addMultPredAbsDiffFC <= 2) / length(tibForPctSimilarPred$addMultPredAbsDiffFC)
# 
# filtTib2 <- filtTib %>%
#   filter(addMultPredAbsDiffFC > 1.0) %>%
#   filter(avgCounts_EtOHnlDensity >= 12) %>%
#   filter(avgTPM_EtOHnlDensity > 2.5)
#   
# p1 <- ggplot(filtTib2, aes(addPredFC, 2^log2fcBoth)) + geom_point() + coord_fixed() + xlab('add pred linear space') + geom_abline(slope=1, intercept=0) + xlim(0, 50) + ylim(0, 50)
# p2 <- ggplot(filtTib2, aes(multPredFC, 2^log2fcBoth)) + geom_point() + coord_fixed() + xlab('mult pred linear space') + geom_abline(slope=1, intercept=0) + xlim(0, 50) + ylim(0, 50)
# p3 <- ggplot(filtTib2, aes(log2(addPredFC), log2fcBoth)) + geom_point() + coord_fixed() + xlab('add pred log space') + geom_abline(slope=1, intercept=0)
# p4 <- ggplot(filtTib2, aes(log2(multPredFC), log2fcBoth)) + geom_point() + coord_fixed() + xlab('mult pred log space') + geom_abline(slope=1, intercept=0)
# grid.arrange(p1,p2,p3,p4, top='differentially expressed genes where additive + multiplicative model \ndiffer by at least a fold-change of 3 and the EtOH TPM is > 2.5')
# 
# 
# # remake original graph (fold changes >4 in both RA and TGFb singleton tx)
# filtTib3 <- filtTib %>%
#   filter(log2fcRA > 2, log2fcTGFb > 2)  %>%
#   filter(isDeGeneRA3d, isDeGeneRA3d)
# this.tpm.cutoff <- 1
# p1 <- ggplot(filtTib3, aes(addPredFC, 2^log2fcBoth, color=avgTPM_EtOHnlDensity > 2)) + geom_point() + coord_fixed() + xlab('add pred linear space') + geom_abline(slope=1, intercept=0) + xlim(0, 100) + ylim(0, 250)
# p2 <- ggplot(filtTib3, aes(multPredFC, 2^log2fcBoth, color=avgTPM_EtOHnlDensity > 2)) + geom_point() + coord_fixed() + xlab('mult pred linear space') + geom_abline(slope=1, intercept=0) + xlim(0, 625) #+ ylim(0, 50)
# p3 <- ggplot(filtTib3, aes(log2(addPredFC), log2fcBoth, color=avgTPM_EtOHnlDensity > 2)) + geom_point() + coord_fixed() + xlab('add pred log space') + geom_abline(slope=1, intercept=0)
# p4 <- ggplot(filtTib3, aes(log2(multPredFC), log2fcBoth, color=avgTPM_EtOHnlDensity > 2)) + geom_point() + coord_fixed() + xlab('mult pred log space') + geom_abline(slope=1, intercept=0)
# grid.arrange(p1,p2,p3,p4, top='differentially expressed genes where log2fcRA > 2 and log2fcTGFb > 2')
# 
# 
# ############### figure generation for thesis committee meeting #2 ###############
# 
# filtTib <- deSeqTibble %>%
#   mutate(a = 2^log2fcRA - 1, b = 2^log2fcTGFb - 1, addPredFC = a + b + 1, multPredFC = 2^(log2fcRA + log2fcTGFb)) %>%
#   mutate(addMultPredDiffFC = multPredFC - addPredFC, addMultPredAbsDiffFC = abs(addMultPredDiffFC))
# 
# minimum.fc.diff <- 1
# min.counts.denominator <- 12
# min.TPM.denominator <- 2.5
# 
# filtTib2 <- filtTib %>%
#   filter(addMultPredAbsDiffFC > minimum.fc.diff) %>%
#   filter(avgCounts_EtOHnlDensity >= min.counts.denominator) %>%
#   filter(avgTPM_EtOHnlDensity > min.TPM.denominator)
# 
# p1 <- ggplot(filtTib2, aes(addPredFC, 2^log2fcBoth)) + geom_point() + coord_fixed() + xlab('additive prediction, fold-change') + ylab('measured fold-change') + geom_abline(slope=1, intercept=0) + xlim(0, 50) + ylim(0, 50)
# p2 <- ggplot(filtTib2, aes(multPredFC, 2^log2fcBoth)) + geom_point() + coord_fixed() + xlab('multiplicative prediction, fold-change') + ylab('measured fold-change') + geom_abline(slope=1, intercept=0) + xlim(0, 50) + ylim(0, 50)
# p3 <- ggplot(filtTib2, aes(log2(addPredFC), log2fcBoth)) + geom_point() + coord_fixed() + xlab('additive prediction, log2-fold-change') + ylab('measured log2-fold-change') + geom_abline(slope=1, intercept=0)
# p4 <- ggplot(filtTib2, aes(log2(multPredFC), log2fcBoth)) + geom_point() + coord_fixed() + xlab('multiplicative prediction, log2-fold-change') + ylab('measured log2-fold-change') + geom_abline(slope=1, intercept=0)
# grid.arrange(p1,p2,p3,p4, top=sprintf('differentially expressed genes where additive + multiplicative model \ndiffer by at least a fold-change of %d and the EtOH TPM is > %.1f', minimum.fc.diff, min.TPM.denominator))
# grid.arrange(p3,p4, top=sprintf('differentially expressed genes where additive + multiplicative model \ndiffer by at least a fold-change of %d and the EtOH TPM is > %.1f', minimum.fc.diff, min.TPM.denominator), ncol=2)

## filtTib2 %>% select(ensg, hgnc.symbol, avgTPM_EtOHnlDensity, avgTPM_TGFb3d, avgTPM_RA3d, avgTPM_Both3d, addMultPredDiffFC) %>% View

# 
# ## Post lab-meeting analysis, use TPM as well
# minTPM <- 10
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc > 1 | ra3d$log2fc > 1) & (tgfb3d$avgTPM_EtOHnlDensity > minTPM)
# qplot(tgfb3d$log2fc[bvec], ra3d$log2fc[bvec]) + coord_fixed() + xlab("log fold change in TGFb") + ylab("log fold change in RA") + ggtitle('Genes with fold-change >4 in both RA and TGFb')
# p1 <- qplot(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec], both3d$log2fc[bvec]) +
#   coord_fixed() + xlab('Multiplicative model prediction log-fold-change') + ylab('measured log-fold change') +
#   xlim(0, 11.4) + ylim(0, 10)
# 
# foldchangeminus1TGFb3d <- 2^(tgfb3d$log2fc[bvec]) - 1
# foldchangeminus1RA3d <- 2^(ra3d$log2fc[bvec]) - 1
# additivePredDiffFromBaseline <- foldchangeminus1TGFb3d + foldchangeminus1RA3d
# additivePredFoldChange <- 1 + additivePredDiffFromBaseline
# additivePredLogFoldChange <- log2(additivePredFoldChange)
# p2 <- qplot(additivePredLogFoldChange, both3d$log2fc[bvec]) +
#   coord_fixed() + xlab('Additive model prediction log-fold-change') + ylab('measured log-fold change') +
#   xlim(0, 7.5) + ylim(0, 10)
# 
# grid.arrange(p1, p2)
# 
# p3 <-  qplot(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]), 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Multiplicative model prediction fold-change') + ylab('measured fold change')
# p4 <-  qplot(2^additivePredLogFoldChange, 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('measured fold change')
# 
# grid.arrange(p1, p2, p3, p4, top='Genes with fold-change >4 in both RA and TGFb')


# 
################################ legacy code from old tibble format style
# # tgfb3d <- filter(deSeqTibble, condition == 'TGFb 3d')
# # ra3d   <- filter(deSeqTibble, condition == 'RA 3d')
# # both3d <- filter(deSeqTibble, condition == 'Both 3d')
# 
# # 
# # # plot for all genes above minimum expression level to be included in deseqtibble
# # qplot(tgfb3d$log2fc, ra3d$log2fc, size=I(0.5))
# # qplot(tgfb3d$log2fc + ra3d$log2fc, both3d$log2fc, size=I(0.5))
# 
# # plot Uri Alon-style "fold change - 1" instaed of log-fold change
# qplot(2^tgfb3d$log2fc[bvec] - 1, 2^ra3d$log2fc[bvec] - 1, size=I(0.5))
# qplot(2^tgfb3d$log2fc[bvec] - 1 + 2^ra3d$log2fc[bvec] - 1, 2^both3d$log2fc[bvec] - 1, size=I(0.5)) + xlim(-2, 2) + ylim(-2, 2)
# 
# 
# # plot for genes that are differentially expressed in both conditions
# bvec <- tgfb3d$isDeGene & ra3d$isDeGene
# qplot(tgfb3d$log2fc[bvec], ra3d$log2fc[bvec], size=I(0.5))
# qplot(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec], both3d$log2fc[bvec], size=I(0.5))
# 
# foldchangeminus1TGFb3d <- 2^(tgfb3d$log2fc[bvec]) - 1
# foldchangeminus1RA3d <- 2^(ra3d$log2fc[bvec]) - 1
# additivePredDiffFromBaseline <- foldchangeminus1TGFb3d + foldchangeminus1RA3d
# additivePredFoldChange <- 1 + additivePredDiffFromBaseline
# additivePredLogFoldChange <- log2(additivePredFoldChange)
# qplot(additivePredInLogSpace, both3d$log2fc[bvec])
# qplot(2^additivePredInLogSpace, 2^both3d$log2fc[bvec])
# 
# # plot for genes that are differentially expressed in both conditions AND increase expression by at least 4 fold
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc > 1.584963 | ra3d$log2fc > 1.584963)
# #bvec <- (tgfb3d$isDeGene & ra3d$isDeGene) & (tgfb3d$log2fc > 2 & ra3d$log2fc > 2)
# qplot(tgfb3d$log2fc[bvec], ra3d$log2fc[bvec]) + coord_fixed() + xlab("log fold change in TGFb") + ylab("log fold change in RA") + ggtitle('Genes with fold-change >4 in both RA and TGFb')
# p1 <- qplot(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec], both3d$log2fc[bvec]) + 
#   coord_fixed() + xlab('Multiplicative model prediction log-fold-change') + ylab('measured log-fold change') +
#   xlim(0, 11.4) + ylim(0, 10)
# 
# foldchangeminus1TGFb3d <- 2^(tgfb3d$log2fc[bvec]) - 1
# foldchangeminus1RA3d <- 2^(ra3d$log2fc[bvec]) - 1
# additivePredDiffFromBaseline <- foldchangeminus1TGFb3d + foldchangeminus1RA3d
# additivePredFoldChange <- 1 + additivePredDiffFromBaseline
# additivePredLogFoldChange <- log2(additivePredFoldChange)
# p2 <- qplot(additivePredLogFoldChange, both3d$log2fc[bvec]) + 
#   coord_fixed() + xlab('Additive model prediction log-fold-change') + ylab('measured log-fold change') +
#   xlim(0, 7.5) + ylim(0, 10)
# 
# grid.arrange(p1, p2)
# 
# p3 <-  qplot(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]), 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Multiplicative model prediction fold-change') + ylab('measured fold change') 
# p4 <-  qplot(2^additivePredLogFoldChange, 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('measured fold change') 
# 
# grid.arrange(p1, p2, p3, p4, top='Genes with fold-change >4 in both RA and TGFb')
# 
# minTPM <- 4
# p5 <-  qplot(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]), 2^both3d$log2fc[bvec], color=(both3d$avgTPM_EtOHnlDensity[bvec]) > minTPM) + coord_fixed() + xlab('Multiplicative model prediction fold-change') + ylab('measured fold change') +
#   ylim(0, 250) + xlim(0, 700) + ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# p6 <-  qplot(2^additivePredLogFoldChange, 2^both3d$log2fc[bvec], color=(both3d$avgTPM_EtOHnlDensity[bvec]) > minTPM) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('measured fold change') +
#   xlim(0, 100) + ylim(0,225) + ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# grid.arrange(p5, p6, top='Genes with fold-change >4 in both RA and TGFb', ncol=2)
# 
# # tgfb3d[bvec,] %>% mutate(multpred=(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec])), addpred=2^additivePredLogFoldChange) %>% View()
# ra3d[bvec,] %>% mutate(multpred=(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec])), addpred=2^additivePredLogFoldChange, measured=2^both3d$log2fc[bvec]) %>% View()
# 
# 
# #p2 <- grid.arrange(grobs=plotsToGrid, top='individual element responses')
# 
# 
# # how different is log fold change from fold-change - 1 when fold-change -1 is in a narrow range?
# rangeToCheck <- seq(0.5, 3, by=0.025) - 1
# qplot(rangeToCheck, log2(rangeToCheck + 1)) + geom_path(mapping=aes(x=rangeToCheck, y=rangeToCheck)) + xlab('fc - 1') + theme(axis.text=element_text(size=20), axis.title=element_text(size=24))
# 
# rangeToCheck <- seq(1, 4.5, by=0.025) - 1
# qplot(rangeToCheck, log2(rangeToCheck + 1)) + geom_path(mapping=aes(x=rangeToCheck, y=rangeToCheck)) + xlab('fc - 1') + theme(axis.text=element_text(size=20), axis.title=element_text(size=24))
# 
# 
# # get linear model coefficients when both genes differentially expressed
# bvec <- tgfb3d$isDeGene & ra3d$isDeGene
# sum(bvec)
# lm(both ~ ra + tgfb +0, tibble(ra=ra3d$log2fc[bvec], tgfb=tgfb3d$log2fc[bvec], both=both3d$log2fc[bvec]))
# 
# # "" when one or the other is differntially expressed
# bvec <- tgfb3d$isDeGene | ra3d$isDeGene
# sum(bvec)
# lm(both ~ ra + tgfb +0, tibble(ra=ra3d$log2fc[bvec], tgfb=tgfb3d$log2fc[bvec], both=both3d$log2fc[bvec]))
# 
# # "" when the "Both" condition has differential expression
# bvec <- both3d$isDeGene
# sum(bvec)
# lm(both ~ ra + tgfb +0, tibble(ra=ra3d$log2fc[bvec], tgfb=tgfb3d$log2fc[bvec], both=both3d$log2fc[bvec]))
# 
# # "" for all genes regardless of diff expression
# lm(both ~ ra + tgfb +0, tibble(ra=ra3d$log2fc, tgfb=tgfb3d$log2fc, both=both3d$log2fc))
# 
# 
# 
# 
# ## Post lab-meeting analysis, use TPM as well
# minTPM <- 10
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc > 1 | ra3d$log2fc > 1) & (tgfb3d$avgTPM_EtOHnlDensity > minTPM)
# qplot(tgfb3d$log2fc[bvec], ra3d$log2fc[bvec]) + coord_fixed() + xlab("log fold change in TGFb") + ylab("log fold change in RA") + ggtitle('Genes with fold-change >4 in both RA and TGFb')
# p1 <- qplot(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec], both3d$log2fc[bvec]) + 
#   coord_fixed() + xlab('Multiplicative model prediction log-fold-change') + ylab('measured log-fold change') +
#   xlim(0, 11.4) + ylim(0, 10)
# 
# foldchangeminus1TGFb3d <- 2^(tgfb3d$log2fc[bvec]) - 1
# foldchangeminus1RA3d <- 2^(ra3d$log2fc[bvec]) - 1
# additivePredDiffFromBaseline <- foldchangeminus1TGFb3d + foldchangeminus1RA3d
# additivePredFoldChange <- 1 + additivePredDiffFromBaseline
# additivePredLogFoldChange <- log2(additivePredFoldChange)
# p2 <- qplot(additivePredLogFoldChange, both3d$log2fc[bvec]) + 
#   coord_fixed() + xlab('Additive model prediction log-fold-change') + ylab('measured log-fold change') +
#   xlim(0, 7.5) + ylim(0, 10)
# 
# grid.arrange(p1, p2)
# 
# p3 <-  qplot(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]), 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Multiplicative model prediction fold-change') + ylab('measured fold change') 
# p4 <-  qplot(2^additivePredLogFoldChange, 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('measured fold change') 
# 
# grid.arrange(p1, p2, p3, p4, top='Genes with fold-change >4 in both RA and TGFb')
# 
# p5 <-  qplot(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]), 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Multiplicative model prediction fold-change') + ylab('measured fold change') +
#   ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# p6 <-  qplot(2^additivePredLogFoldChange, 2^both3d$log2fc[bvec]) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('measured fold change') +
#   ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# grid.arrange(p5 + xlim(0,45) + ylim(0,45), p6 + xlim(0,45) + ylim(0,45), ncol=2)
# 
# 
# #plot additive vs. multiplicative prediction
# qplot(2^additivePredLogFoldChange, 2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec])) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('Multiplicative model prediction') +
#   xlim(0, 100) + ylim(0,225) + ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# 
# qplot(2^additivePredLogFoldChange, 2^(both3d$log2fc[bvec])) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('Multiplicative model prediction') +
#   xlim(0, 100) + ylim(0,225) + ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# 
# qplot(2^(both3d$log2fc[bvec]), 2^(both3d$log2fc[bvec]), color=(both3d$avgTPM_EtOHnlDensity[bvec] > 4)) + coord_fixed() + xlab('Additive model prediction fold-change') + ylab('Multiplicative model prediction') +
#   xlim(0, 100) + ylim(0,225) + ggtitle('Genes with fold-change >4 in both RA and TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# 
# qplot(additivePredLogFoldChange, (tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec])) + coord_fixed() + xlab('Additive model prediction log-fold-change') + ylab('Multiplicative model prediction') + 
#   ggtitle('Genes with fold-change >1 in either RA or TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# 
# qplot(additivePredFoldChange, 2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec])) + coord_fixed() + xlab('Additive model prediction log-fold-change') + ylab('Multiplicative model prediction') + 
#   ggtitle('Genes with fold-change >1 in either RA or TGFb') + theme(axis.text = element_text(size=16), axis.title = element_text(size=16)) +
#   xlim(0, 100) + ylim (0, 200)
# 
# qplot(additivePredFoldChange, 2^(both3d$log2fc[bvec]), color=(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]) > 49)) + 
#   coord_fixed() + theme(axis.text = element_text(size=16), axis.title = element_text(size=16))
# 
# qplot(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]), 2^(both3d$log2fc[bvec]),  color=(2^(tgfb3d$log2fc[bvec] + ra3d$log2fc[bvec]) > 49)) + coord_fixed() + 
#   theme(axis.text = element_text(size=16), axis.title = element_text(size=16)) 
# 
# 
# ## Post lab-meeting analysis, use TPM as well -- use linear model coefficients too.
# ## Post lab-meeting analysis, use TPM as well
# minTPM <- 5
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$avgTPM_EtOHnlDensity > minTPM)
# measuredFoldChange <- 2^(both3d$log2fc[bvec])
# 
# foldchangeminus1TGFb3d <- 2^(tgfb3d$log2fc[bvec]) - 1
# foldchangeminus1RA3d <- 2^(ra3d$log2fc[bvec]) - 1
# foldchangeminus1Measured <- measuredFoldChange - 1
# 
# addLm <- lm(foldchangeminus1Measured ~ foldchangeminus1TGFb3d + foldchangeminus1RA3d + 0, tibble(foldchangeminus1TGFb3d, foldchangeminus1RA3d, foldchangeminus1Measured))
# additivePredDiffFromBaseline <- coefficients(addLm)[1]*foldchangeminus1TGFb3d + coefficients(addLm)[2]*foldchangeminus1RA3d
# additivePredFoldChange <- 1 + additivePredDiffFromBaseline
# additivePredLogFoldChange <- log2(additivePredFoldChange)
# p1 <- qplot(additivePredFoldChange, 2^both3d$log2fc[bvec]) + 
#   coord_fixed() + xlab('Additive model prediction fold-change') + ylab('measured log-fold change')
# 
# logfcTGFb <- tgfb3d$log2fc[bvec]
# logfcRA   <- ra3d$log2fc[bvec]
# measuredLogFc <- both3d$log2fc[bvec]
# 
# multLm <- lm(measuredLogFc ~ logfcTGFb + logfcRA + 0, tibble(measuredLogFc, logfcRA, logfcTGFb))
# multiplicativePredLogFoldChange <- coefficients(addLm)[1]*logfcTGFb + coefficients(addLm)[2]*logfcRA
# multiplicativePredFoldChange <- 2 ^ multiplicativePredLogFoldChange
# 
# p2 <- qplot(multiplicativePredFoldChange, measuredFoldChange) + 
#   coord_fixed() + xlab('Multiplicative model prediction fold-change') + ylab('measured log-fold change')
# 
# p3 <- qplot(multiplicativePredFoldChange, additivePredFoldChange)
# 
# 
# diffIndsToTest <- which(abs(multiplicativePredFoldChange - additivePredFoldChange)/measuredFoldChange > 0.5)
# p4 <- qplot(multiplicativePredFoldChange[diffIndsToTest], measuredFoldChange[diffIndsToTest]) + coord_fixed()
# p5 <- qplot(additivePredFoldChange[diffIndsToTest], measuredFoldChange[diffIndsToTest]) + coord_fixed()
# 
# grid.arrange(p1 + xlim(-0.5,15) + ylim(-0.5,15), p2 + xlim(-0.5,15) + ylim(-0.5,15), ncol=2)
# grid.arrange(p1, p2, ncol=2)
# 
# 
# 
# ###
# # examine genes that are differentially expressed AND increase expression by at least 3 fold
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc > 1.584963 | ra3d$log2fc > 1.584963)
# both3d[bvec,] %>% 
#   mutate(a=(avgTPM_RA3d-avgTPM_EtOHnlDensity), b=(avgTPM_TGFb3d-avgTPM_EtOHnlDensity), predDiffMagn=abs(a*b/avgTPM_EtOHnlDensity), predDiff=(a*b/avgTPM_EtOHnlDensity), additivePredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b)), multPredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b + predDiff))) %>% 
#   dplyr::select(ensg, avgTPM_Both3d, avgTPM_RA3d, avgTPM_TGFb3d, avgTPM_EtOHnlDensity, a, b, predDiffMagn, additivePredDiff, multPredDiff) %>% 
#   View
# 
# both3d[bvec,] %>% 
#   mutate(a=(avgTPM_RA3d-avgTPM_EtOHnlDensity), b=(avgTPM_TGFb3d-avgTPM_EtOHnlDensity), predDiffMagn=abs(a*b/avgTPM_EtOHnlDensity), predDiff=(a*b/avgTPM_EtOHnlDensity), additivePredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b)), multPredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b + predDiff))) %>% 
#   dplyr::select(ensg, avgTPM_Both3d, avgTPM_RA3d, avgTPM_TGFb3d, avgTPM_EtOHnlDensity, a, b, predDiffMagn, additivePredDiff, multPredDiff) %>% 
#   filter(!is.infinite(predDiffMagn), avgTPM_EtOHnlDensity > 2, a * b < 0) %>%
#   ggplot() + 
#   geom_point(mapping=aes(additivePredDiff/avgTPM_Both3d, multPredDiff/avgTPM_Both3d, color=(avgTPM_EtOHnlDensity > 5))) +
#   xlim(-4, 4) + ylim(-4,4)
# 
# both3d[bvec,] %>% 
#   mutate(measuredFC= avgTPM_Both3d/avgTPM_EtOHnlDensity, a=(avgTPM_RA3d-avgTPM_EtOHnlDensity), b=(avgTPM_TGFb3d-avgTPM_EtOHnlDensity), predDiffMagn=abs(a*b/avgTPM_EtOHnlDensity), predDiff=(a*b/avgTPM_EtOHnlDensity), additivePredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b))/avgTPM_Both3d, multPredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b + predDiff))/avgTPM_Both3d) %>% 
#   dplyr::select(ensg, avgTPM_Both3d, avgTPM_RA3d, avgTPM_TGFb3d, avgTPM_EtOHnlDensity, a, b, measuredFC, predDiffMagn, additivePredDiff, multPredDiff) %>% 
#   filter(!is.infinite(predDiffMagn), avgTPM_EtOHnlDensity > 2) %>%
#   View
# 
# #when TGFb and RA have effects in opposing directions, multiplicative model fits better
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc > 1.584963 | ra3d$log2fc > 1.584963)
# both3d[bvec,] %>% 
#   mutate(a=(avgTPM_RA3d-avgTPM_EtOHnlDensity), b=(avgTPM_TGFb3d-avgTPM_EtOHnlDensity), predDiffMagn=abs(a*b/avgTPM_EtOHnlDensity), predDiff=(a*b/avgTPM_EtOHnlDensity), additivePredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b)), multPredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b + predDiff))) %>% 
#   dplyr::select(ensg, avgTPM_Both3d, avgTPM_RA3d, avgTPM_TGFb3d, avgTPM_EtOHnlDensity, a, b, predDiffMagn, additivePredDiff, multPredDiff) %>% 
#   filter(!is.infinite(predDiffMagn), avgTPM_EtOHnlDensity > 2, a * b < 0) %>%
#   ggplot() + 
#   geom_point(mapping=aes(additivePredDiff/avgTPM_Both3d, multPredDiff/avgTPM_Both3d, color=(avgTPM_EtOHnlDensity > 5))) +
#   xlim(-3, 3) + ylim(-3, 3)
# 
# #when TGFb and RA have effects in SAME directions, additive model fits better?
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc > 1.584963 | ra3d$log2fc > 1.584963)
# both3d[bvec,] %>% 
#   mutate(a=(avgTPM_RA3d-avgTPM_EtOHnlDensity), b=(avgTPM_TGFb3d-avgTPM_EtOHnlDensity), predDiffMagn=abs(a*b/avgTPM_EtOHnlDensity), predDiff=(a*b/avgTPM_EtOHnlDensity), additivePredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b)), multPredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b + predDiff))) %>% 
#   dplyr::select(ensg, avgTPM_Both3d, avgTPM_RA3d, avgTPM_TGFb3d, avgTPM_EtOHnlDensity, a, b, predDiffMagn, additivePredDiff, multPredDiff) %>% 
#   filter(!is.infinite(predDiffMagn), avgTPM_EtOHnlDensity > 2, a * b > 0) %>%
#   ggplot() + 
#   geom_point(mapping=aes(additivePredDiff/avgTPM_Both3d, multPredDiff/avgTPM_Both3d, color=(avgTPM_EtOHnlDensity > 5))) + coord_fixed()
#   #xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
# 
# #what about when both signals downregulate the gene?
# bvec <- (tgfb3d$isDeGene | ra3d$isDeGene) & (tgfb3d$log2fc < 1.584963 | ra3d$log2fc < 1.584963)
# both3d[bvec,] %>% 
#   mutate(a=(avgTPM_RA3d-avgTPM_EtOHnlDensity), b=(avgTPM_TGFb3d-avgTPM_EtOHnlDensity), predDiffMagn=abs(a*b/avgTPM_EtOHnlDensity), predDiff=(a*b/avgTPM_EtOHnlDensity), additivePredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b)), multPredDiff=(avgTPM_Both3d-(avgTPM_EtOHnlDensity + a + b + predDiff))) %>% 
#   dplyr::select(ensg, avgTPM_Both3d, avgTPM_RA3d, avgTPM_TGFb3d, avgTPM_EtOHnlDensity, a, b, predDiffMagn, additivePredDiff, multPredDiff) %>% 
#   filter(!is.infinite(predDiffMagn), avgTPM_EtOHnlDensity > 2, a * b > 0) %>%
#   ggplot() + 
#   geom_point(mapping=aes(additivePredDiff/avgTPM_Both3d, multPredDiff/avgTPM_Both3d, color=(avgTPM_EtOHnlDensity > 5))) + 
#   xlim(-5, 5) + ylim(-5,5)


