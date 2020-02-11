library(tidyverse)
library(here)
library(VGAM)  # required for rfoldnorm function

n.iterations.per.input.observation <- 5  # if 1, null distribution will have the same number of observations as the input data, if 2 it'll be 2x, etc...
n.experiment.replicates <- 3
max.n.sample.attempts.to.be.higher.than.control <- 100000

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  dist.for.peaks.vs.genes  <- "genes" 
  min.fc.diff.mult.add.for.c.histogram <- 0.2
  min.raw.val.diff.for.c.histogram       <- 1
  add.vs.mult.null.model   <- "mixture" # choose "additive" or "multiplicative" or "mixture"
  add.mult.mixture.frac.add <- 0.60
  # input files -- use upregulated peaks or genes
  input.table  <- read_tsv(here("extractedData", "DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv"))
  # input.table  <- read_tsv(here("extractedData", "differentialAtacPeaks_mergedist250_peakwidth150_minNormFrags30_minFoldChange1.5.annotated.upregulated.tsv"))
  output.file.prefix <- here("plots", "null_distributions",
                             sprintf("null_distribution_upreg_%s_%s_model_%0.2f_", dist.for.peaks.vs.genes, add.vs.mult.null.model, add.mult.mixture.frac.add))
} else {
  input.table                            <- read_tsv(cmdargs[1])
  min.fc.diff.mult.add.for.c.histogram   <- as.numeric(cmdargs[2])
  min.raw.val.diff.for.c.histogram       <- as.numeric(cmdargs[3])
  dist.for.peaks.vs.genes                <- cmdargs[4]  # choose "peaks" or "genes"
  add.vs.mult.null.model                 <- cmdargs[5]  # choose "additive" or "multiplicative" or "mixture"
  add.mult.mixture.frac.add              <- as.numeric(cmdargs[6])
  output.file.prefix                     <- paste0(cmdargs[7], '/', sprintf("null_distribution_upreg_%s_%s_model_%0.2f_", dist.for.peaks.vs.genes, add.vs.mult.null.model, add.mult.mixture.frac.add))
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

# for each gene, draw ?3 values from the folded normal distribution, average them, then report the measured c and d values 
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

if (add.vs.mult.null.model == "additive") {
  add.mult.mixture.frac.add <- 1
} else if (add.vs.mult.null.model == "multiplicative") {
  add.mult.mixture.frac.add <- 0
} else if (add.vs.mult.null.model == "mixture") {
  add.mult.mixture.frac.add <- add.mult.mixture.frac.add
}


for (dose in c("low", "med", "high")) {
  set.seed(0)
  new.cvals.all              <- c()
  new.dvals.all              <- c()
  new.add.mult.raw.diffs.all <- c()
  new.add.mult.fc.diffs.all  <- c()
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
    }
    new.cvals              <- sapply(new.sampled.val.list, function(x) x[[1]])
    new.dvals              <- sapply(new.sampled.val.list, function(x) x[[2]])
    new.add.mult.raw.diffs <- sapply(new.sampled.val.list, function(x) x[[3]])
    new.add.mult.fc.diffs  <- sapply(new.sampled.val.list, function(x) x[[4]])
    
    new.cvals.all <- c(new.cvals.all, new.cvals)
    new.dvals.all <- c(new.dvals.all, new.dvals)
    new.add.mult.raw.diffs.all <- c(new.add.mult.raw.diffs.all, new.add.mult.raw.diffs)
    new.add.mult.fc.diffs.all  <- c(new.add.mult.fc.diffs.all, new.add.mult.fc.diffs)
  }
  nulldist.result.tib <- tibble(cval = new.cvals.all, dval = new.dvals.all, 
                                raw.add.mult.diff = new.add.mult.raw.diffs.all, fc.add.mult.diff = new.add.mult.fc.diffs.all)
  
  p1 <- nulldist.result.tib %>% 
    filter(raw.add.mult.diff >= min.raw.val.diff.for.c.histogram,
           fc.add.mult.diff >= min.fc.diff.mult.add.for.c.histogram) %>% 
      ggplot(aes(cval)) +
        geom_histogram(bins = 100) + xlim(-2.5, 5) + 
      ggtitle(paste0("cvals from ", add.vs.mult.null.model, " null model, ", dose, " dose", ",\nminAddMultFcDiff = ", 
                     min.fc.diff.mult.add.for.c.histogram, ", min raw value diff ", min.raw.val.diff.for.c.histogram)) +
      theme_minimal()
  ggsave(paste0(output.file.prefix, "cval_plot_", dose, "_dose.svg"), plot = p1)
  
  p2 <- qplot(new.dvals.all, bins = 100) + xlim(-2.5, 5) + ggtitle(paste0("dvals from ", add.vs.mult.null.model, " null model, ", dose, " dose")) +
    theme_minimal()
  ggsave(paste0(output.file.prefix, "dval_plot_", dose, "_dose.svg"), plot = p2)
}
