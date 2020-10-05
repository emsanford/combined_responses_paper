library(tidyverse)
library(here)

inputdatapath = here("extractedData", "immunofluorescence_image_data", "curated_image_nuclei_IF_statistics.tsv")
outputfolder  = here("plots", "immunofluorescence_plots")
nucleus_tib   = read_tsv(inputdatapath)
time.factor.order <- c("0hrPseudo", "40min", "2hr", "4hr", "12hr", "24hr", "72hr")
condition.factor.order <- c("control", "ra_treated", "tgfb_treated")

nucleus_tib$timepoint <- factor(nucleus_tib$timepoint, levels = time.factor.order)
nucleus_tib$condition <- factor(nucleus_tib$condition, levels = condition.factor.order)
nucleus_tib$replicate_ID <- factor(nucleus_tib$replicate_ID, levels = c("IFT3", "IFT1"))

outlier.range.pSMAD2.lower <- -1000
outlier.range.pSMAD2.upper <-  7000

nucleus_tib_pSMAD2 <- nucleus_tib %>% filter(antibody_stained == "pSMAD2", 
                                             annulus_subtracted_intensity >= outlier.range.pSMAD2.lower, 
                                             annulus_subtracted_intensity <= outlier.range.pSMAD2.upper)

activated.nucleus.threshold <- 1000

outlier.range.RARA.lower <- -2000
outlier.range.RARA.upper <- 10000

nucleus_tib_RARA <- nucleus_tib %>% filter(antibody_stained == "RARA", 
                                             annulus_subtracted_intensity >= outlier.range.RARA.lower, 
                                             annulus_subtracted_intensity <= outlier.range.RARA.upper)

saved.plots <- list()
# show distributions for pSMAD2 with dotty box plots
for (replicate in c("IFT1", "IFT3")) {
  this.plot.tibble <- nucleus_tib_pSMAD2 %>% filter(replicate_ID == replicate, timepoint != "0hrPseudo")
  this.plot <- this.plot.tibble %>%
    ggplot(aes(x = timepoint, y = annulus_subtracted_intensity, fill = condition)) + 
    #ggrastr::rasterise(geom_point(position = position_jitterdodge(), size = 0.3, alpha = 1, color = "dark grey"), dpi = 180) +
    geom_point(position = position_jitterdodge(), size = 0.4, alpha = 1, color = "dark grey") +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=c("#656565", "#39CF29", "#2396D6")) + 
    theme_classic(base_size = 12) + theme(legend.position="none") +
    ggtitle(sprintf("%s", replicate))
  
  # now add threshold lines to the plot
  num.ref.lines <- length(unique(this.plot.tibble$timepoint))
  for (ii in 1:num.ref.lines) {
    timepoint.for.threshold   <- sort(unique(this.plot.tibble$timepoint))[ii]
    median.val.this.timepoint <- median(this.plot.tibble %>% filter(timepoint == timepoint.for.threshold) %>% pull("annulus_subtracted_intensity"))
    antibody.on.threshold <- median.val.this.timepoint + activated.nucleus.threshold
    spacer_radius <- 0.4
    this.plot <- this.plot + geom_segment(x = ii - spacer_radius, xend = ii + spacer_radius, y = antibody.on.threshold, yend = antibody.on.threshold)
  }
  
  print(this.plot)
  saved.plots[[length(saved.plots) + 1]] <- this.plot
}

ggsave(paste0(outputfolder, '/boxplots_pSMAD2_rep1.svg'),  plot = saved.plots[[1]], width = 6, height = 3.65)
ggsave(paste0(outputfolder, '/boxplots_pSMAD2_rep2.svg'),  plot = saved.plots[[2]], width = 6, height = 3.65)
ggsave(paste0(outputfolder, '/boxplots_pSMAD2_rep1.png'),  plot = saved.plots[[1]] + ggtitle("") + theme(axis.title.x = element_blank(),
                                                                                                         axis.title.y = element_blank(), 
                                                                                                         axis.text.x = element_blank(),
                                                                                                         axis.text.y = element_blank()), width = 6.5, height = 3.65)
ggsave(paste0(outputfolder, '/boxplots_pSMAD2_rep2.png'),  plot = saved.plots[[2]] + ggtitle("") + theme(axis.title.x = element_blank(),
                                                                                                         axis.title.y = element_blank(), 
                                                                                                         axis.text.x = element_blank(),
                                                                                                         axis.text.y = element_blank()), width = 6.5, height = 3.65)

for (replicate in c("IFT1", "IFT3")) {
  this.plot.tibble <- nucleus_tib_RARA %>% filter(replicate_ID == replicate, timepoint != "0hrPseudo")
  this.plot <- this.plot.tibble %>%
    ggplot(aes(x = timepoint, y = annulus_subtracted_intensity, fill = condition)) + 
    #ggrastr::rasterise(geom_point(position = position_jitterdodge(), size = 0.3, alpha = 1, color = "dark grey"), dpi = 180) +
    geom_point(position = position_jitterdodge(), size = 0.4, alpha = 1, color = "dark grey") +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=c("#656565", "#39CF29", "#2396D6")) + 
    theme_classic(base_size = 12) + theme(legend.position="none") +
    ggtitle(sprintf("%s", replicate))
  
  
  print(this.plot)
  saved.plots[[length(saved.plots) + 1]] <- this.plot
}

ggsave(paste0(outputfolder, '/boxplots_RARA_rep1.svg'),  plot = saved.plots[[3]], width = 6, height = 3.65)
ggsave(paste0(outputfolder, '/boxplots_RARA_rep2.svg'),  plot = saved.plots[[4]], width = 6, height = 3.65)
ggsave(paste0(outputfolder, '/boxplots_RARA_rep1.png'),  plot = saved.plots[[3]] + ggtitle("") + theme(axis.title.x = element_blank(),
                                                                                                         axis.title.y = element_blank(), 
                                                                                                         axis.text.x = element_blank(),
                                                                                                         axis.text.y = element_blank()), width = 6.5, height = 3.65)
ggsave(paste0(outputfolder, '/boxplots_RARA_rep2.png'),  plot = saved.plots[[4]] + ggtitle("") + theme(axis.title.x = element_blank(),
                                                                                                         axis.title.y = element_blank(), 
                                                                                                         axis.text.x = element_blank(),
                                                                                                         axis.text.y = element_blank()), width = 6.5, height = 3.65)



# test percent above "on" threshold 
for (thisrep in c("IFT1", "IFT3")) {
  vector.of.frac.above.thresh <- c()
  time.vector <- c()
  condition.vector <- c()
  timepoints.this.rep <-sort(unique(filter(nucleus_tib_pSMAD2, replicate_ID == thisrep)$timepoint))
  for (time.point in timepoints.this.rep) {
    medianval = median(nucleus_tib_pSMAD2 %>% filter(replicate_ID == thisrep, timepoint == time.point) %>% pull("annulus_subtracted_intensity"))
    for (cond1 in  c("control", "ra_treated", "tgfb_treated")) {
      condsubset = filter(nucleus_tib_pSMAD2, replicate_ID == thisrep, timepoint == time.point, condition == cond1)
      num1 <- sum(condsubset$annulus_subtracted_intensity > (medianval + activated.nucleus.threshold))
      den1 <- length(condsubset$annulus_subtracted_intensity)
      # print(sprintf("%s %s %0.3f", time.point, cond1, num1 / den1))
      if (is.nan(num1/den1)) {
        next
      }
      vector.of.frac.above.thresh <- c(vector.of.frac.above.thresh, num1 / den1)
      time.vector <- c(time.vector, time.point)
      condition.vector <- c(condition.vector, cond1)
    }
  }
  
  time.factor.vals <- factor(time.vector, levels = time.factor.order)
  condition.as.factors <- factor(condition.vector, levels = condition.factor.order)
  tib.above.thresh.plot <- tibble(condition = condition.as.factors, 
                                  pct_above_thresh = vector.of.frac.above.thresh,
                                  timepoint = time.factor.vals)
  
  ggplot(tib.above.thresh.plot, aes(x = timepoint, y = pct_above_thresh, fill = condition, color = condition, group = condition)) +
    geom_line() +
    theme_classic() + 
    ggtitle(thisrep)
    
  # now do bootstrap sampling for fraction of nuclei above threshold
  set.seed(0)
  n.bootstrap.samples <- 1000
  bootstrap.vector.of.frac.above.thresh <- c()
  bootstrap.time.vector <- c()
  bootstrap.condition.vector <- c()
  for (ii in 1:n.bootstrap.samples) {
    nuc_tib_bootstrap_sample <- sample_n(nucleus_tib_pSMAD2, nrow(nucleus_tib_pSMAD2), replace = TRUE)
    
    for (time.point in timepoints.this.rep) {
      medianval = median(nuc_tib_bootstrap_sample %>% filter(replicate_ID == thisrep, timepoint == time.point) %>% pull("annulus_subtracted_intensity"))
      for (cond1 in  c("control", "ra_treated", "tgfb_treated")) {
        condsubset = filter(nuc_tib_bootstrap_sample, replicate_ID == thisrep, timepoint == time.point, condition == cond1)
        num1 <- sum(condsubset$annulus_subtracted_intensity > (medianval + activated.nucleus.threshold))
        den1 <- length(condsubset$annulus_subtracted_intensity)
        # print(sprintf("%s %s %0.3f", time.point, cond1, num1 / den1))
        this.pct.above.thresh <- num1/den1 
        if (is.nan(this.pct.above.thresh)) {
          next
        }
        bootstrap.vector.of.frac.above.thresh <- c(bootstrap.vector.of.frac.above.thresh, this.pct.above.thresh)
        bootstrap.time.vector <- c(bootstrap.time.vector, time.point)
        bootstrap.condition.vector <- c(bootstrap.condition.vector, cond1)
      }
    }
  }
  
  bootstrap.tib.above.thresh.plot <- tibble(condition = factor(bootstrap.condition.vector,  levels = condition.factor.order),
                                            pct_above_thresh = bootstrap.vector.of.frac.above.thresh,
                                            timepoint = factor(bootstrap.time.vector, levels = time.factor.order))
  bootstrap.upper.intervals <- c()
  bootstrap.lower.intervals <- c()
  bootstrap.ci <- .95
  for (time.point in timepoints.this.rep) {
    for (cond1 in  c("control", "ra_treated", "tgfb_treated")) {
      this.fracs.above.thresh <- bootstrap.tib.above.thresh.plot %>% filter(condition == cond1, timepoint == time.point) %>% pull("pct_above_thresh")
      this.measured.value <- pull(filter(tib.above.thresh.plot, condition == cond1, timepoint == time.point), "pct_above_thresh")[1]
      if (is.na(this.measured.value)) {
        next
      }
      bootstrap.upper.intervals <- c(bootstrap.upper.intervals, 2 * this.measured.value - quantile(this.fracs.above.thresh, .05))
      bootstrap.lower.intervals <- c(bootstrap.lower.intervals, 2 * this.measured.value - quantile(this.fracs.above.thresh, .95))
    }
  }
  
  tib.above.thresh.plot[["bootstrap_ci_upper"]] <- bootstrap.upper.intervals
  tib.above.thresh.plot[["bootstrap_ci_lower"]] <- bootstrap.lower.intervals
  
  p.traces <- ggplot(tib.above.thresh.plot, aes(x = timepoint, y = pct_above_thresh, fill = condition, color = condition, 
                                    group = condition, ymin = bootstrap_ci_lower, ymax = bootstrap_ci_upper)) +
    geom_line() +
    geom_point(size = I(1)) +
    geom_errorbar(color = "gray", width = 0.05) +
    theme_classic() + theme(legend.position="none") +
    ggtitle(thisrep)
  
  print(p.traces)
  saved.plots[[length(saved.plots) + 1]] <- p.traces
}

# install.packages('devtools')
# remotes::install_github('VPetukhov/ggrastr')

ggsave(paste0(outputfolder, '/lineplots_pSMAD2_rep1.svg'), plot = saved.plots[[5]], width = 6, height = 2.5)
ggsave(paste0(outputfolder, '/lineplots_pSMAD2_rep2.svg'), plot = saved.plots[[6]], width = 6, height = 2.5)


# old code graveyard

# # plot both replicates on same axes with different color
# p1 <- nucleus_tib %>%
#   #filter(replicate_ID == "IFT1") %>%
#   ggplot(aes(annulus_subtracted_intensity, fill = replicate_ID)) +
#   geom_histogram(aes(y = ..density..), binwidth = 100) + 
#   xlim(-1000, 7000) +
#   facet_grid(vars(condition), vars(timepoint)) + 
#   xlab("annulus-subtracted nuclear intensity histograms") +
#   theme_classic()
# 
# 
# nucleus_tib_rep1_only <- filter(nucleus_tib, replicate_ID == "IFT1",  antibody_stained == "pSMAD2")
# nucleus_tib_rep2_only <- filter(nucleus_tib, replicate_ID == "IFT3",  antibody_stained == "pSMAD2")
# 
# xlim.lower <- -1000
# xlim.upper <-  7000
# 
# p2 <- ggplot(NULL) +
#   geom_histogram(data = nucleus_tib_rep1_only, mapping = aes(annulus_subtracted_intensity, y = ..density..), color = "dark gray",  binwidth = 25) + 
#   xlim(xlim.lower, xlim.upper) +
#   facet_grid(vars(condition), vars(timepoint)) + 
#   theme_classic()
# 
# p3 <- ggplot(nucleus_tib_rep2_only, aes(annulus_subtracted_intensity, y = ..density..)) +
#   geom_histogram(color = "dark gray",  binwidth = 25) +
#   xlim(xlim.lower, xlim.upper) +
#   facet_grid(vars(condition), vars(timepoint))
#
# ggsave(paste0(outputfolder, '/test1.svg'), plot = p2)
# ggsave(paste0(outputfolder, '/test2.svg'), plot = p3)
