library(tidyverse)
library(here)
library(VGAM)  # required for rfoldnorm function

source(here('extractionScripts', 'util.R'))

n.experiment.replicates <- 3
max.n.sample.attempts.to.be.higher.than.control <- 100000

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  dist.for.peaks.vs.genes  <- "genes" 
  min.fc.diff.mult.add.for.c.histogram   <- 0
  min.raw.val.diff.for.c.histogram       <- 0
  # input files -- use upregulated peaks or genes
  # input.table  <- read_tsv(here("extractedData", "DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv"))
  input.table  <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
  # input.table      <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv")
  output.file.prefix <- paste0("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4/plots/gene_integration_summary_plots/null_distributions/", dist.for.peaks.vs.genes, "_")
} else {
  input.table                            <- read_tsv(cmdargs[1])
  min.fc.diff.mult.add.for.c.histogram   <- as.numeric(cmdargs[2])
  min.raw.val.diff.for.c.histogram       <- as.numeric(cmdargs[3])
  dist.for.peaks.vs.genes                <- cmdargs[4]  # choose "peaks" or "genes"
  output.file.prefix                     <- paste0(cmdargs[5], '/', sprintf("null_distribution_upreg_%s_", dist.for.peaks.vs.genes))
}

if (dist.for.peaks.vs.genes == "genes") {
  n.iterations.per.input.observation <- 50  # if 1, null distribution will have the same number of observations as the input data, if 2 it'll be 2x, etc...
} else if (dist.for.peaks.vs.genes == "peaks") {
  n.iterations.per.input.observation <- 10  
}


# using the estimates for the mean and CV from our data set, draw a new set of observations under the constraint that the values for 
# tgfb and ra must be greater than the etoh value (this is a model for upregulation only)
drawNewValueSetFromDistribution <- function(etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.est.CV, add.vs.mult.pred) {
  # 1. calc new EtOH value
  if (etoh.mean.value == 0) {  # edge case: use pseudocount when etoh mean value is zero (otherwise get an error in rfoldnorm for variance est.)
    new.etoh.mean.value <- 0.0001
  } else {
    new.etoh.mean.value <- mean(rfoldnorm(n.experiment.replicates, mean = etoh.mean.value, sd = etoh.mean.value * etoh.est.CV))
  }
  # 2. calc new RA and TGFb effects (they must be greater than the etoh value. if not, resample them)
  new.ra.mean.value   <- 0
  n.attempts          <- 0
  while (new.ra.mean.value < new.etoh.mean.value) {
    new.ra.mean.value <- mean(rfoldnorm(n.experiment.replicates, mean = ra.mean.value, sd = ra.mean.value * ra.est.CV))
    n.attempts <- n.attempts + 1
    if (n.attempts > max.n.sample.attempts.to.be.higher.than.control) {
      print("pathological edge case for ra? too many attempts to sample a value higher than the control value")
      print(c(new.etoh.mean.value, etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.est.CV, add.vs.mult.pred))
      print("doubling CV term for ra ...")
      ra.est.CV <- 2 * ra.est.CV
      n.attempts <- 0
    }
  }
  new.tgfb.mean.value   <- 0
  n.attempts <- 0
  while (new.tgfb.mean.value < new.etoh.mean.value) {
    new.tgfb.mean.value <- mean(rfoldnorm(n.experiment.replicates, mean = tgfb.mean.value, sd = tgfb.mean.value * tgfb.est.CV))
    n.attempts <- n.attempts + 1
    if (n.attempts > max.n.sample.attempts.to.be.higher.than.control) {
      print("pathological edge case for tgfb? too many attempts to sample a value higher than the control value")
      print(c(new.etoh.mean.value, etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.est.CV, add.vs.mult.pred))
      print("doubling CV term for tgfb ...")
      tgfb.est.CV <- 2 * tgfb.est.CV
      n.attempts <- 0
    }
  }
  # 3. select add vs. mult model prediction and calculate new expected mean value
  stopifnot(add.vs.mult.pred %in% c("add", "mult"))
  if (add.vs.mult.pred == "add") {
    both.signals.mean.value <- makeAddModelPrediction(new.etoh.mean.value, new.ra.mean.value, new.tgfb.mean.value)
  } else if (add.vs.mult.pred == "mult") {
    both.signals.mean.value <- makeMultModelPrediction(new.etoh.mean.value, new.ra.mean.value, new.tgfb.mean.value)
  }
  # 4. calculate new add / mult prediction
  new.both.signals.mean.value   <- 0
  n.attempts <- 0
  while (new.both.signals.mean.value < new.etoh.mean.value) {
    new.both.signals.mean.value <- mean(rfoldnorm(n.experiment.replicates, mean = both.signals.mean.value, sd = both.signals.mean.value * both.est.CV))
    n.attempts <- n.attempts + 1
    if (n.attempts > max.n.sample.attempts.to.be.higher.than.control) {
      print("pathological edge case for both signals? too many attempts to sample a value higher than the control value")
      print(c(new.etoh.mean.value, etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.signals.mean.value, both.est.CV, add.vs.mult.pred))
      stopifnot(F)
    }
  }
  # 5. return new c and d values
  new.add.prediction  <- makeAddModelPrediction(new.etoh.mean.value, new.ra.mean.value, new.tgfb.mean.value)
  new.mult.prediction <- makeMultModelPrediction(new.etoh.mean.value, new.ra.mean.value, new.tgfb.mean.value)
  # the c value reflects the "integration spectrum" of the final result. when c = 0, it's additive, when c = 1, it's multiplicative
  ra.effect   <- new.ra.mean.value - new.etoh.mean.value
  tgfb.effect <- new.tgfb.mean.value - new.etoh.mean.value
  new.c.value <- new.etoh.mean.value * (new.both.signals.mean.value - new.add.prediction) / (ra.effect * tgfb.effect)
  # the d value is the fold-change difference between the measured expression and the additive model prediction
  new.d.value <- (new.both.signals.mean.value / new.etoh.mean.value) - (new.add.prediction / new.etoh.mean.value)
  
  raw.val.mult.add.diff    <- new.mult.prediction - new.add.prediction
  foldchange.mult.add.diff <- new.mult.prediction/new.etoh.mean.value - new.add.prediction/new.etoh.mean.value
  
  return(list(new.c.value, new.d.value, raw.val.mult.add.diff, foldchange.mult.add.diff))
}

makeAddModelPrediction <- function(etoh.value, ra.value, tgfb.value) {
  ra.effect   <- ra.value - etoh.value
  tgfb.effect <- tgfb.value - etoh.value
  return(etoh.value + ra.effect + tgfb.effect)
}

makeMultModelPrediction <- function(etoh.value, ra.value, tgfb.value) {
  ra.foldchange <- ra.value / etoh.value
  tgfb.foldchange <- tgfb.value / etoh.value
  return(etoh.value * ra.foldchange * tgfb.foldchange)
}


grand.tib.for.cval.plot.grid <- NULL
grand.tib.for.dval.plot.grid <- NULL
for (add.vs.mult.null.model in c("additive", "multiplicative", "mixture")) {
  if (add.vs.mult.null.model == "additive") {
    add.mult.mixture.frac.add <- 1
  } else if (add.vs.mult.null.model == "multiplicative") {
    add.mult.mixture.frac.add <- 0
  } 
  
  for (dose in c("low", "med", "high")) {
    if (add.vs.mult.null.model == "mixture") {
      if (dist.for.peaks.vs.genes == "genes") {
        n.add.genes  <- sum(input.table[[paste0("integrationCategory-", dose, "-dose")]] == "additive")
        n.mult.genes <- sum(input.table[[paste0("integrationCategory-", dose, "-dose")]] == "multiplicative")
        add.mult.mixture.frac.add <- n.add.genes / (n.add.genes + n.mult.genes)
      } else if (dist.for.peaks.vs.genes == "peaks") {
        n.add.peaks  <- sum(input.table[[paste0("peak_integrationCategory-", dose, "-dose")]] == "additive")
        n.mult.peaks <- sum(input.table[[paste0("peak_integrationCategory-", dose, "-dose")]] == "multiplicative")
        add.mult.mixture.frac.add <- n.add.peaks / (n.add.peaks + n.mult.peaks)
      }
    }
    set.seed(0)
    new.cvals.all              <- c()
    new.dvals.all              <- c()
    new.add.mult.raw.diffs.all <- c()
    new.add.mult.fc.diffs.all  <- c()
    selected.null.dists.all    <- c()
    for (jj in 1:n.iterations.per.input.observation) {
      if (dist.for.peaks.vs.genes == "peaks") {
        # define the column names in peaks we are using normalized fragment counts
        etoh.mean.values   <- input.table[["EtOH-nlDensity-avgNormFragmentCounts"]]
        etoh.est.CVs       <- input.table[["EtOH_peak_cvEstPooled"]]
        ra.mean.values     <- input.table[[paste0("RA-", dose, "-avgNormFragmentCounts")]]
        ra.est.CVs         <- input.table[["RA_peak_cvEstPooled"]]
        tgfb.mean.values   <- input.table[[paste0("TGFb-", dose, "-avgNormFragmentCounts")]]
        tgfb.est.CVs       <- input.table[["TGFb_peak_cvEstPooled"]]
        both.est.CVs       <- input.table[["TGFb-and-RA_peak_cvEstPooled"]]  # use the measured CV for the both condition. the mean will be determined by the definition of the null model and the individual mean signal effects
      } else if (dist.for.peaks.vs.genes == "genes") {
        # define the column names in peaks we are using transcripts per million
        etoh.mean.values   <- input.table[["EtOH-nlDensity_avgTPM"]]
        etoh.est.CVs       <- input.table[["EtOH_cvEstPooled"]]
        ra.mean.values     <- input.table[[paste0("RA-", dose, "_avgTPM")]]
        ra.est.CVs         <- input.table[["RA_cvEstPooled"]]
        tgfb.mean.values   <- input.table[[paste0("TGFb-", dose, "_avgTPM")]]
        tgfb.est.CVs       <- input.table[["TGFb_cvEstPooled"]]
        both.est.CVs       <- input.table[["TGFb-and-RA_cvEstPooled"]] # use the measured CV for the both condition. the mean will be determined by the definition of the null model and the individual mean signal effects
      }
  
      n.observations <- length(etoh.mean.values)
      # iterate over the rows of the input table, use a weighted coin model to select rows one-by-one, and call drawNewValueSetFromDistribution
      new.sampled.val.list <- list()
      selected.null.dists  <- c()
      for (ii in 1:n.observations) {
        use.add.null.model <- runif(1) <= add.mult.mixture.frac.add
        if (use.add.null.model) {
          this.null.model <- "add"
        } else {
          this.null.model <- "mult"
        }
        new.sampled.values <- drawNewValueSetFromDistribution(etoh.mean.values[ii], etoh.est.CVs[ii], 
                                                              ra.mean.values[ii], ra.est.CVs[ii], 
                                                              tgfb.mean.values[ii], tgfb.est.CVs[ii], 
                                                              both.est.CVs[ii], 
                                                              this.null.model)
        new.sampled.val.list[[ii]] <- new.sampled.values
        selected.null.dists <- c(selected.null.dists, this.null.model)
      }
      new.cvals              <- sapply(new.sampled.val.list, function(x) x[[1]])
      new.dvals              <- sapply(new.sampled.val.list, function(x) x[[2]])
      new.add.mult.raw.diffs <- sapply(new.sampled.val.list, function(x) x[[3]])
      new.add.mult.fc.diffs  <- sapply(new.sampled.val.list, function(x) x[[4]])
      
      new.cvals.all <- c(new.cvals.all, new.cvals)
      new.dvals.all <- c(new.dvals.all, new.dvals)
      new.add.mult.raw.diffs.all <- c(new.add.mult.raw.diffs.all, new.add.mult.raw.diffs)
      new.add.mult.fc.diffs.all  <- c(new.add.mult.fc.diffs.all, new.add.mult.fc.diffs)
      selected.null.dists.all <- c(selected.null.dists.all, selected.null.dists)
    }
  
    # make plot for c values
    bin.leftmost  <- -3
    bin.rightmost <-  5
    bin.step.size <-  0.125
    plot.width    <-  8
    plot.height   <-  5
     
    categorical.values <- selected.null.dists.all  # to do--add these if interested. could show which distribution they came from (add vs. mult)
    hist.values        <- new.cvals.all
    cval_hist_output   <- makeHistogramOfValues(hist.values, categorical.values, bin.leftmost, bin.rightmost,
                                                   bin.step.size, paste0(add.vs.mult.null.model, ", ", dose, " dose, c-values, addmixfrac = ", add.mult.mixture.frac.add), 
                                                   xlabel = "c-value", ylabel = "counts", color.by.category = T, y.axis.units = "counts")
    stackedBarHistTibCvals <- cval_hist_output[[1]]
    cval.hist.tib          <- cval_hist_output[[2]]
    cval.hist.tib[["dose"]] <- dose
    cval.hist.tib[["null_model"]] <- add.vs.mult.null.model
    grand.tib.for.cval.plot.grid <- rbind(grand.tib.for.cval.plot.grid, cval.hist.tib)
    
    ggsave(paste0(output.file.prefix, "cval_plot_", add.vs.mult.null.model, "_model_", dose, "_dose.svg"), plot = stackedBarHistTibCvals, width = plot.width, height = plot.height)
    
    # make plot for d values
    bin.step.size   <-  0.05
    bin.leftmost    <- -2
    bin.rightmost   <-  2
    plot.width    <-  8
    plot.height   <-  5
    
    categorical.values <- selected.null.dists.all  # to do--add these if interested. could show which distribution they came from (add vs. mult)
    hist.values        <- new.dvals.all
    dval_hist_output <- makeHistogramOfValues(hist.values, categorical.values, bin.leftmost, bin.rightmost,
                                                    bin.step.size, paste0(add.vs.mult.null.model, ", ", dose, " dose, d-values, addmixfrac = ", add.mult.mixture.frac.add), 
                                                    xlabel = "d-value", ylabel = "counts", color.by.category = T, y.axis.units = "counts")
    
    stackedBarHistTibDvals <- dval_hist_output[[1]]
    dval.hist.tib          <- dval_hist_output[[2]]
    dval.hist.tib[["dose"]] <- dose
    dval.hist.tib[["null_model"]] <- add.vs.mult.null.model
    grand.tib.for.dval.plot.grid <- rbind(grand.tib.for.dval.plot.grid, dval.hist.tib)
    
    ggsave(paste0(output.file.prefix, "dval_plot_", add.vs.mult.null.model, "_model_", dose, "_dose.svg"), plot = stackedBarHistTibDvals, width = plot.width, height = plot.height)
  }
}

grand.tib.for.cval.plot.grid[["dose"]]       <- factor(grand.tib.for.cval.plot.grid[["dose"]], levels = c("low", "med", "high"))
grand.tib.for.cval.plot.grid[["null_model"]] <- factor(grand.tib.for.cval.plot.grid[["null_model"]], levels = c("additive", "multiplicative", "mixture"))
n.dose.nullmodel.pairs = 9
reduced.grand.tib.cvals <- grand.tib.for.cval.plot.grid %>% 
  group_by(intConstantHhistBin, intCategory, dose,  null_model) %>% 
  mutate(n_this_bin = n(), freq_this_bin = n_this_bin / (nrow(grand.tib.for.cval.plot.grid) / n.dose.nullmodel.pairs)) %>%
  ungroup() %>%
  unique()
  
gridhistogramplot1 <- grand.tib.for.cval.plot.grid %>%
  filter(null_model == "mixture") %>%
  ggplot(aes(x = intConstantHhistBin, fill = intCategory)) +
  geom_bar(stat="count", width = bin.step.size * .90) +
  facet_grid(~dose) +
  theme_classic()

gridhistogramplot2 <- ggplot(reduced.grand.tib.cvals, aes(x = intConstantHhistBin, y = freq_this_bin, fill = intCategory)) +
  geom_bar(stat="identity", width = bin.step.size * .90) +
  facet_grid(~dose) +
  theme_classic()

ggsave(paste0(output.file.prefix, "cval_grid_counts.svg"), plot = gridhistogramplot1, width = 18, height = 8.15)
ggsave(paste0(output.file.prefix, "cval_grid_frequencies.svg"), plot = gridhistogramplot2, width = 18, height = 8.15)
saveRDS(grand.tib.for.cval.plot.grid, file = paste0(output.file.prefix, "cval_grid_grand_tibble.rds"))

# to do: print percent above c = 2 for each null distribution:
c.cutoff <- 4
for (dosage in c("low", "med", "high")) {
  for (model in c("additive", "multiplicative", "mixture")) {
    freq.above.c.2 <- reduced.grand.tib.cvals %>%
      filter(null_model == model, dose == dosage) %>%
      filter(intConstantHhistBin >= c.cutoff) %>%
      pull("freq_this_bin") %>%
      sum()
    
    print(sprintf("%s %s %0.3f", model, dosage, freq.above.c.2))
  }
}




grand.tib.for.dval.plot.grid[["dose"]]       <- factor(grand.tib.for.dval.plot.grid[["dose"]], levels = c("low", "med", "high"))
grand.tib.for.dval.plot.grid[["null_model"]] <- factor(grand.tib.for.dval.plot.grid[["null_model"]], levels = c("additive", "multiplicative", "mixture"))
n.dose.nullmodel.pairs = 9
reduced.grand.tib <- grand.tib.for.dval.plot.grid %>% 
  group_by(intConstantHhistBin, intCategory, dose,  null_model) %>% 
  mutate(n_this_bin = n(), freq_this_bin = n_this_bin / (nrow(grand.tib.for.dval.plot.grid) / n.dose.nullmodel.pairs)) %>%
  ungroup() %>%
  unique()

gridhistogramplot1 <- ggplot(grand.tib.for.dval.plot.grid, aes(x = intConstantHhistBin, fill = intCategory)) +
  geom_bar(stat="count", width = bin.step.size * .90) +
  facet_grid(dose ~ null_model) +
  theme_classic()

gridhistogramplot2 <- ggplot(reduced.grand.tib, aes(x = intConstantHhistBin, y = freq_this_bin, fill = intCategory)) +
  geom_bar(stat="identity", width = bin.step.size * .90) +
  facet_grid(dose ~ null_model) +
  theme_classic()

histogramplot3 <- ggplot(filter(reduced.grand.tib, null_model == "additive", dose == "low"), aes(x = intConstantHhistBin, y = freq_this_bin, fill = intCategory)) +
  geom_bar(stat="identity", width = bin.step.size * .90) +
  theme_classic()

histogramplot4 <- ggplot(filter(reduced.grand.tib, null_model == "additive", dose == "med"), aes(x = intConstantHhistBin, y = freq_this_bin, fill = intCategory)) +
  geom_bar(stat="identity", width = bin.step.size * .90) +
  theme_classic()

histogramplot5 <- ggplot(filter(reduced.grand.tib, null_model == "additive", dose == "high"), aes(x = intConstantHhistBin, y = freq_this_bin, fill = intCategory)) +
  geom_bar(stat="identity", width = bin.step.size * .90) +
  theme_classic()

library(patchwork)
rowhistogramplot1 <- histogramplot3 + histogramplot4 + histogramplot5
plot.width      <-  16
plot.height     <-  8
ggsave(paste0(output.file.prefix, "dval_row_frequencies.svg"), plot = rowhistogramplot1, width = plot.width * 3, height = plot.height)

ggsave(paste0(output.file.prefix, "dval_grid_counts.svg"), plot = gridhistogramplot1, width = 18, height = 8.15)
ggsave(paste0(output.file.prefix, "dval_grid_frequencies.svg"), plot = gridhistogramplot2, width = 18, height = 8.15)
saveRDS(grand.tib.for.dval.plot.grid, file = paste0(output.file.prefix, "dval_grid_grand_tibble.rds"))

# to do: print percent above c = 2 for each null distribution:
d.cutoff <- 1.5
for (dosage in c("low", "med", "high")) {
  for (model in c("additive", "multiplicative", "mixture")) {
    freq.above.c.2 <- reduced.grand.tib %>%
      filter(null_model == model, dose == dosage) %>%
      filter(intConstantHhistBin >= d.cutoff) %>%
      pull("freq_this_bin") %>%
      sum()
    
    print(sprintf("dval: %s %s %0.5f", model, dosage, freq.above.c.2))
  }
}

