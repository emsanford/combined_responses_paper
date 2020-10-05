# 1. make the original histogram: ( copy-pasted code from geneIntegrationSummaryPieChartsAndHistograms.R )

library(tidyverse)
library(here)
library(VGAM)  # required for rfoldnorm function
library(patchwork)


cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  siUpregGenes     <- read_tsv("/Users/emsanford/Dropbox (RajLab)/Shared_Eric/SIgnal_Integration/Analysis_SI2-SI4_github_testing/signal_integration_paper_scripts/extractedData/DeSeqOutputAllConds.annotated.upregulatedGeneSet.tsv")
  output.folder    <- here('plots', 'secondary_c_value_peak_analysis_plots')
  output.folder2   <- here('plots', 'gene_integration_summary_plots')
} else {
  siUpregGenes     <- read_tsv(cmdargs[1])
  output.folder    <- cmdargs[2]
  output.folder2   <- cmdargs[3]
}

control.tpm.zero.rows <- which(siUpregGenes$`EtOH-nlDensity_avgTPM` == 0)
siUpregGenes <- siUpregGenes[-control.tpm.zero.rows, ] # remove rows with control TPM = 0 because the erroneously estimate c to be zero
print(paste0("removed ", length(control.tpm.zero.rows), " genes with control TPM = 0"))

bin.step.size   <-  0.125
bin.leftmost    <- -4
bin.rightmost   <-  5
single.plot.width    <-   6
single.plot.height   <-   6
bin.radius    <- bin.step.size / 2
bin.midpoints <- seq(bin.leftmost + bin.step.size, bin.rightmost, by = bin.step.size) - bin.radius

plot.vertical.shrinkage.factor <- 0.60

n.experiment.replicates <- 3
max.n.sample.attempts.to.be.higher.than.control <- 1000

stackedBarHistogram.location.prefix <- paste0(output.folder, '/gene_integration_mode_stackedBarHistogram_AddSimDataSubtracted')

n_upreg_genes <- nrow(siUpregGenes)

####### make stacked histogram plot showing signal integration constant & colored by frequency of each category in the plot

# function: reassign weird categories that may come up when a specific dose integrates in a different direction to "uncategorized" 
mapCatsToReducesCatSet <- function(cat.values) {
  allowedCategories <- c("ambiguous", "sub-additive", "additive", "between-add-and-mult", "multiplicative", "super-multiplicative")
  cat.values[! cat.values %in% allowedCategories] <- "uncategorized"
  return(factor(cat.values, levels = rev(c("uncategorized", allowedCategories))))
}

sumTable <- function(contingency.table.one.dim) {
  cum.sum <- 0
  for (ii in contingency.table.one.dim) {
    cum.sum <- cum.sum + ii
  }
  return(cum.sum)
} 

# functions for generating simulated data distributions
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
      #print("pathological edge case for ra? too many attempts to sample a value higher than the control value")
      #print(c(new.etoh.mean.value, etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.est.CV, add.vs.mult.pred))
      #print("doubling CV term for ra ...")
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
      #print("pathological edge case for tgfb? too many attempts to sample a value higher than the control value")
      #print(c(new.etoh.mean.value, etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.est.CV, add.vs.mult.pred))
      #print("doubling CV term for tgfb ...")
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
      #print("pathological edge case for both signals? too many attempts to sample a value higher than the control value")
      #print(c(new.etoh.mean.value, etoh.mean.value, etoh.est.CV, ra.mean.value, ra.est.CV, tgfb.mean.value, tgfb.est.CV, both.signals.mean.value, both.est.CV, add.vs.mult.pred))
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


create.simulated.set.of.cvals <- function(geneTibble, add.mult.mixture.frac.add, dosage, n.iterations.per.input.observation = 10, seed = 0) {
  set.seed(seed)
  new.cvals.all      <- c()
  etoh.mean.values   <- geneTibble[["EtOH-nlDensity_avgTPM"]]
  etoh.est.CVs       <- geneTibble[["EtOH_cvEstPooled"]]
  ra.mean.values     <- geneTibble[[paste0("RA-", dosage, "_avgTPM")]]
  ra.est.CVs         <- geneTibble[["RA_cvEstPooled"]]
  tgfb.mean.values   <- geneTibble[[paste0("TGFb-", dosage, "_avgTPM")]]
  tgfb.est.CVs       <- geneTibble[["TGFb_cvEstPooled"]]
  both.est.CVs       <- geneTibble[["TGFb-and-RA_cvEstPooled"]]
  
  n.observations <- length(etoh.mean.values)
  new.sampled.val.list <- list()
  new.cvals.all        <- c()
  
  for (jj in 1:n.iterations.per.input.observation) {
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
    new.cvals     <- sapply(new.sampled.val.list, function(x) x[[1]])
    new.cvals.all <- c(new.cvals.all, new.cvals)
  }
  
  return(new.cvals.all)
}

create.histogram.y.values <- function(values, bin.midpoints) {
  n.vals <- length(values)
  n.bins <- length(bin.midpoints)
  bin.counts <- rep(0, n.bins)
  
  for (ii in 1:n.vals) {
    this.val <- values[ii]
    distvec <- abs(this.val - bin.midpoints)
    lowest.bin.distance <- min(distvec)
    bin.index <- which(distvec == lowest.bin.distance)[1]
    bin.counts[bin.index] <- bin.counts[bin.index] + 1
  }
  
  return(bin.counts)
}


## here endeth functions for generating simulated data distributions


# dosages.to.test <- c("high")

dosages.to.test <- c("low", "med", "high")
for (k in 1:3) {
  dosage <- dosages.to.test[k]
  
  additive.fraction.empirical.estimates       <- c(11.6/(11.3 + 11.6), 15.1/(11.7 + 15.1), 15.1/(11.9 + 15.1))
  multiplicative.fraction.empirical.estimates <- c(11.3/(11.3 + 11.6), 11.7/(11.7 + 15.1), 11.9/(11.9 + 15.1))
  fraction.add.mult.or.ambiguous.estimates    <- c(.473 + .116 + .113, .437 + .151 + .117, .414 + .151 + .119)
  
  additive.fraction.empirical.estimate       <- additive.fraction.empirical.estimates[k] 
  multiplicative.fraction.empirical.estimate <- multiplicative.fraction.empirical.estimates[k]
  fraction.add.mult.or.ambiguous.estimate    <- fraction.add.mult.or.ambiguous.estimates[k]
  
  # parameter determines if we subtract additive, multiplicative, or mixture simulated histogram from observed data
  # distributions.to.subtract <- c("additive", "multiplicative", "mixture")
  distributions.to.subtract <- c("additive")
  
  # regardless of which distribution we subtract, we need to calculate the simulated data mixture model to 
  # in case we use the method of scaling the mixture model to the peak of the data before subtracting it...
  n.iterations.per.input.observation <- 250
  # make simulated data for additive, multiplicative, and mixture model. we will use these in the main figure to show different simulated data distributions.
  # to outline the primary peak at c = 0 in the observed data before showing the secondary peak at/close to c = 1 after performing the subtraction.
  additive.model.c.vals <- create.simulated.set.of.cvals(siUpregGenes, 1, dosage, n.iterations.per.input.observation) 
  additive.model.c.vals.table <- create.histogram.y.values(additive.model.c.vals, bin.midpoints)
  additive.model.c.vals.table.as.freqs <- additive.model.c.vals.table / sum(additive.model.c.vals.table)
  multiplicative.model.c.vals <- create.simulated.set.of.cvals(siUpregGenes, 0, dosage, n.iterations.per.input.observation) 
  multiplicative.model.c.vals.table <- create.histogram.y.values(multiplicative.model.c.vals, bin.midpoints)
  multiplicative.model.c.vals.table.as.freqs <- multiplicative.model.c.vals.table / sum(multiplicative.model.c.vals.table)
  mixture.model.c.vals <- create.simulated.set.of.cvals(siUpregGenes, additive.fraction.empirical.estimate, dosage, n.iterations.per.input.observation)
  mixture.model.c.vals.table <- create.histogram.y.values(mixture.model.c.vals, bin.midpoints)
  mixture.model.c.vals.table.as.freqs <- mixture.model.c.vals.table / sum(mixture.model.c.vals.table)
  hist.values.measured <- pull(siUpregGenes, paste0("integrationConstant-", dosage))
  measured.cvalTable <- create.histogram.y.values(hist.values.measured, bin.midpoints)
  measured.cvalTable.as.freqs <- measured.cvalTable / sum(measured.cvalTable)
  
  tib.for.plotting.obs.data.and.single.models <- tibble(bar_locations = bin.midpoints,
                                                        add_simdata = additive.model.c.vals.table.as.freqs,
                                                        mult_simdata = multiplicative.model.c.vals.table.as.freqs,
                                                        obs_data = measured.cvalTable.as.freqs)
  
  p11 <- ggplot(tib.for.plotting.obs.data.and.single.models, aes(bar_locations, obs_data)) + 
    geom_bar(stat="identity") + 
    ylim(0, 0.04) + 
    geom_vline(xintercept = 1) + 
    geom_vline(xintercept = 0) + 
    ylab("frequency") + xlab("c value") + 
    ggtitle(paste0("measured combined response factor values\n", dosage, " dose")) + theme_classic()
  p12 <- ggplot(tib.for.plotting.obs.data.and.single.models, aes(bar_locations, add_simdata)) + 
    geom_bar(stat="identity") + 
    ylim(0, 0.10) + 
    geom_vline(xintercept = 1) + 
    geom_vline(xintercept = 0) + 
    ylab("frequency") + xlab("c value") + 
    ggtitle(paste0("additive simulated c vals\n", dosage, " dose")) + theme_classic()
  p13 <- ggplot(tib.for.plotting.obs.data.and.single.models, aes(bar_locations, mult_simdata)) + 
    geom_bar(stat="identity") + 
    ylim(0, 0.10) + 
    geom_vline(xintercept = 1) + 
    geom_vline(xintercept = 0) + 
    ylab("frequency") + xlab("c value") + 
    ggtitle(paste0("multiplicative simulated c vals\n", dosage, " dose")) + theme_classic()
    
  # need second tibble to show add vs. mult colors as separate (just subtract mult from add values w appropriate scaling)
  tib.for.two.color.mixture.model <- tibble(hist_bar_locations = rep(bin.midpoints, 2),
                                            mixture_data_split = c(mixture.model.c.vals.table.as.freqs - (multiplicative.fraction.empirical.estimate * multiplicative.model.c.vals.table.as.freqs), 
                                                                   mixture.model.c.vals.table.as.freqs - (additive.fraction.empirical.estimate * additive.model.c.vals.table.as.freqs)),
                                            simulated_data_model = factor(c(rep("add", length(bin.midpoints)), rep("mult", length(bin.midpoints))), levels = c("mult", "add")))
  
  p14 <- ggplot(tib.for.two.color.mixture.model, aes(x = hist_bar_locations, y = mixture_data_split, fill = simulated_data_model, 
                                                     order = simulated_data_model)) + 
    geom_bar(stat = "identity") + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + 
    ylim(0, 0.07) + ylab('frequency') + xlab("c value") + ggtitle(paste0("add/mult mixture model c vals", "\n", dosage, " dose")) + 
    theme_classic() + theme(legend.position = "none") 
  
  # plot panel of additive model, multiplicative model, mixture model, and observed data
  p15 <- p12 + p13 + p14 + p11 + plot_layout(nrow = 1)
  ggsave(paste0(output.folder2, '/composite_simulated_and_obs_cval_histograms_', dosage, '_dose.svg'), plot = p15, width = 16, height = 4)

  for (simulated.distribution.to.subtract in distributions.to.subtract) {
    
    if (simulated.distribution.to.subtract == "additive") {
      percent.of.combined.mixture.distribution.to.subtract <- additive.fraction.empirical.estimate
      simulated.cvals <- additive.model.c.vals
    } else if (simulated.distribution.to.subtract == "multiplicative") {
      percent.of.combined.mixture.distribution.to.subtract <- multiplicative.fraction.empirical.estimate
      simulated.cvals <- multiplicative.model.c.vals
    } else if (simulated.distribution.to.subtract == "mixture") {
      percent.of.combined.mixture.distribution.to.subtract <- 1
      simulated.cvals <- mixture.model.c.vals
    }

    closest.cval.bin.to.right.of.zero <- min(bin.midpoints[bin.midpoints >= 0])
    closest.cval.bin.to.right.of.one <- min((bin.midpoints - 1)[(bin.midpoints - 1) >= 0]) + 1

    index.of.closest.bin.to.right.of.zero <- which(bin.midpoints == closest.cval.bin.to.right.of.zero)
    index.of.closest.bin.to.left.of.zero <- index.of.closest.bin.to.right.of.zero - 1
    index.of.closest.bin.to.right.of.one <- which(bin.midpoints == closest.cval.bin.to.right.of.one)
    index.of.closest.bin.to.left.of.one <- index.of.closest.bin.to.right.of.one - 1
    indices.near.zero.and.one <- c(index.of.closest.bin.to.left.of.zero, index.of.closest.bin.to.right.of.zero, 
                                   index.of.closest.bin.to.left.of.one, index.of.closest.bin.to.right.of.one)
    
    measured.raw.val.near.zero <- measured.cvalTable.as.freqs[index.of.closest.bin.to.right.of.zero]
    
    ############### now generate plot for simulated data   ###############  
    # using the estimates for the mean and CV from our data set, draw a new set of observations under the constraint that the values for 
    # tgfb and ra must be greater than the etoh value (this is a model for upregulation only)
    
    simulated.cvalTable <- create.histogram.y.values(simulated.cvals, bin.midpoints) 
    simulated.cvalTable.as.freqs <- simulated.cvalTable / sum(simulated.cvalTable)
    simulated.cvalHistTibble <- tibble(yvals = simulated.cvalTable.as.freqs, xvals = bin.midpoints)
    simulated.raw.val.near.zero <- simulated.cvalTable.as.freqs[index.of.closest.bin.to.right.of.zero]
    
    # calculate a factor to scale the simulated data size to the measured data before subtracting it
    scalemethod <- "match mixture model peak height"
    if (scalemethod == "match mixture model peak height") {
      mixture.model.values.to.match <- mixture.model.c.vals.table.as.freqs[indices.near.zero.and.one]
      obs.values.to.match <- measured.cvalTable.as.freqs[indices.near.zero.and.one]
      
      mixture.scaling.factor.model <- lm(formula = obs.values.to.match ~ 0 + mixture.model.values.to.match)
      mixture.scaling.factor <- mixture.scaling.factor.model$coefficients
      
      SIMULATED.DATA.SCALING.FACTOR <- mixture.scaling.factor * percent.of.combined.mixture.distribution.to.subtract
    } else if (scalemethod == "use arbitrary scaling factor") {
      ARBITRARY.SCALING.FACTOR <- .60
      SIMULATED.DATA.SCALING.FACTOR  <- ARBITRARY.SCALING.FACTOR * fraction.add.mult.or.ambiguous.estimate * percent.of.combined.mixture.distribution.to.subtract # currently calculated as additive + multiplicative + ambiguous percentage
    }
    
    simulated.data.scaled.for.subtraction <- simulated.cvalTable.as.freqs * SIMULATED.DATA.SCALING.FACTOR
    observed.cval.hist.minus.simulated.cval.hist <- measured.cvalTable.as.freqs - simulated.data.scaled.for.subtraction

    # make tibble containing observed data, simulated data, and subtracted data
    addModelSubtractedTib <- tibble(diff_bar_values = observed.cval.hist.minus.simulated.cval.hist, 
                                    hist_bar_locations = bin.midpoints,
                                    measured_bar_values = measured.cvalTable.as.freqs,
                                    simulated_bar_values = simulated.cvalTable.as.freqs,
                                    scaled_simulated_bar_values = simulated.data.scaled.for.subtraction,
                                    scaled_mixture_model = mixture.model.c.vals.table.as.freqs * mixture.scaling.factor)
    
    addMultMixtureTib <- tibble(hist_bar_locations = rep(bin.midpoints, 2),
                                simulated_scaled_data = c(simulated.data.scaled.for.subtraction, mixture.model.c.vals.table.as.freqs * mixture.scaling.factor - simulated.data.scaled.for.subtraction),
                                simulated_data_component = c(rep("add", length(bin.midpoints)), rep("mult", length(bin.midpoints))))
    addMultMixtureTib$simulated_data_component <- factor(addMultMixtureTib$simulated_data_component, levels = c("mult", "add"))
    
    # fit gaussian to to subtracted data
    y <- addModelSubtractedTib$diff_bar_values
    x <- addModelSubtractedTib$hist_bar_locations
    mu_est    <- 0
    sigsq_est <- 1
    #scrub outliers
    y <- y[2:(length(y) - 1)]
    x <- x[2:(length(x) - 1)]
    # fit without vertical displacement constant
    testfit <- nls(y ~ coeff * dnorm(x, mu, sigsq), start = list(mu = mu_est, sigsq = sigsq_est, coeff = 1))
    # fit with vertical displacement constant
    # testfit <- nls(y ~ coeff * dnorm(x, mu, sigsq) + constant, start = list(mu = mu_est, sigsq = sigsq_est, coeff = 1, constant = 0))
    print(paste0(dosage, " dose params:"))
    print(testfit$m$getAllPars())
    peak.center.string = sprintf("%0.03f", as.numeric(as.character(testfit$m$getAllPars()[1])))
    x.range.for.curve.plot <- seq(bin.leftmost + bin.step.size, bin.rightmost, by = bin.step.size * 0.1) - bin.radius
    fitted.curve.tib <- tibble(fit.x = x.range.for.curve.plot, fit.y = predict(testfit, list(x = x.range.for.curve.plot)))
    
    # plot original data, simulated data, and subtracted data with gaussian fit overlaid
    ylim_lower <- -0.0005
    ylim_upper <- 0.041
    p1 <- ggplot(addModelSubtractedTib, aes(hist_bar_locations, measured_bar_values)) + geom_bar(stat="identity") + ylim(ylim_lower, ylim_upper) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + ylab("measured frequency") + xlab("c value") + ggtitle(paste0("measured combined response factor values\n", dosage, " dose")) + theme_classic()
    p2 <- ggplot(addModelSubtractedTib, aes(hist_bar_locations, scaled_simulated_bar_values)) + geom_bar(stat="identity") + ylim(ylim_lower, ylim_upper) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + ylab("simulated frequency (additive model, scaled)") + xlab("c value") + ggtitle(paste0("simulated c value distribution to subtract, ", simulated.distribution.to.subtract, "\n", dosage, " dose")) + theme_classic()
    p3 <- ggplot() + geom_bar(data = addModelSubtractedTib, aes(hist_bar_locations, diff_bar_values), stat="identity") + 
      geom_line(data = fitted.curve.tib, aes(x = fit.x, y = fit.y), color = "blue") +
      ylim(ylim_lower, ylim_upper) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + 
      ylab("measured - simulated frequency") + xlab("c value") + theme_classic(base_size = 16) + 
      ggtitle(paste0(dosage, " dose")) + ggtitle(paste0("after subtracting simulated dist.\n", "fitted peak center = ", peak.center.string)) +
      theme_classic()
    #p4 <- ggplot(addModelSubtractedTib, aes(hist_bar_locations, scaled_mixture_model)) + geom_bar(stat="identity") + ylim(ylim_lower, ylim_upper) + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + ylab("scaled simulated add/mult mixture model") + xlab("c value") + theme_classic(base_size = 16) + ggtitle(paste0("scaled simulated add/mult mixture model", "\n", dosage, " dose")) + theme_classic()
    p5 <- ggplot(addMultMixtureTib, aes(x = hist_bar_locations, y = simulated_scaled_data, fill = simulated_data_component, order = simulated_data_component)) + geom_bar(stat = "identity") + geom_vline(xintercept = 1) + geom_vline(xintercept = 0) + theme_classic() + theme(legend.position = "none") + ylim(ylim_lower, ylim_upper) + ggtitle(paste0("scaled simulated add/mult mixture model", "\n", dosage, " dose"))
    
    pcomposite <- p1 + p5 + p3 + plot_layout(nrow = 1)
    print(pcomposite)
    
    ggsave(paste0(output.folder, '/observed_cval_hist_', dosage, '_dose.svg'), plot = p1, width = single.plot.width, height = single.plot.height)
    ggsave(paste0(output.folder, '/simulated_', simulated.distribution.to.subtract, '_scaled_for_subtraction_', dosage, '_dose.svg'), plot = p2, width = single.plot.width, height = single.plot.height)
    ggsave(paste0(output.folder, '/subtracted_simulated_', simulated.distribution.to.subtract, '_from_obs_data_', dosage, '_dose.svg'), plot = p3, width = single.plot.width, height = single.plot.height)
    ggsave(paste0(output.folder, '/simulated_mixture_model_scaled_', dosage, '_dose.svg'), plot = p5, width = single.plot.width, height = single.plot.height)
    ggsave(paste0(output.folder, '/composite_subtracted_simulated_', simulated.distribution.to.subtract, '_from_obs_plots_', dosage, '_dose.svg'), plot = pcomposite, width = single.plot.width * 3, height = single.plot.height * plot.vertical.shrinkage.factor)
  }
}


saved.plot.list.logpvals <- list()
# make a bunch of histograms under additive model with N = 1 for sampling, first try for medium dose, 
for (mm in 1:3) {
  dosage <- c('low', 'med', 'high')[mm]

  bin.step.size <- 0.125
  bin.radius    <- bin.step.size / 2
  bin.midpoints <- seq(bin.leftmost + bin.step.size, bin.rightmost, by = bin.step.size) - bin.radius
  
  fraction.additive <- 1
  simulated.additive.cvalTable <- create.histogram.y.values(create.simulated.set.of.cvals(siUpregGenes, fraction.additive, dosage, n.iterations.per.input.observation = 50), bin.midpoints)
  simulated.additive.cvalTable.as.freqs <- simulated.additive.cvalTable / sum(simulated.additive.cvalTable)
  
  hist.values.measured <- pull(siUpregGenes, paste0("integrationConstant-", dosage))
  measured.cvalTable <- create.histogram.y.values(hist.values.measured, bin.midpoints)
  measured.cvalTable.as.freqs <- measured.cvalTable / sum(measured.cvalTable)
  
  index.of.closest.bin.to.right.of.zero <- which(bin.midpoints == closest.cval.bin.to.right.of.zero)
  index.of.closest.bin.to.left.of.zero <- index.of.closest.bin.to.right.of.zero - 1
  indices.near.zero <- c(index.of.closest.bin.to.right.of.zero, index.of.closest.bin.to.left.of.zero)
  
  additive.model.values.to.match <- simulated.additive.cvalTable.as.freqs[indices.near.zero]
  obs.values.to.match <- measured.cvalTable.as.freqs[indices.near.zero]
  
  sim.additive.scaling.factor.model <- lm(formula = obs.values.to.match ~ 0 + additive.model.values.to.match)
  sim.additive.scaling.factor <- sim.additive.scaling.factor.model$coefficients
  
  p8 <- ggplot(tibble(hist_bar_locations = bin.midpoints, simulated_bar_values = sim.additive.scaling.factor * simulated.additive.cvalTable.as.freqs), aes(hist_bar_locations, simulated_bar_values)) + geom_bar(stat="identity") + ylim(ylim_lower, ylim_upper) + geom_vline(xintercept = 0) + ylab("sim scaled add frequency") + xlab("c value") + ggtitle(paste0("best add model fit scaled\n", dosage, " dose")) + theme_classic()
  ggsave(paste0(output.folder, '/simulated_additive_model_peak_scaled_for_pval_calculations', dosage, '_dose.svg'), plot = p8, width = single.plot.width, height = single.plot.height)
  
  
  hist.values.measured <- pull(siUpregGenes, paste0("integrationConstant-", dosage))
  measured.cvalTable <- create.histogram.y.values(hist.values.measured, bin.midpoints)
  measured.cvalTable.as.freqs <- measured.cvalTable / sum(measured.cvalTable)
  p9 <- ggplot(tibble(hist_bar_locations = bin.midpoints, observed_bar_values = measured.cvalTable.as.freqs), aes(hist_bar_locations, observed_bar_values)) + geom_bar(stat="identity") + ylim(ylim_lower, ylim_upper) + geom_vline(xintercept = 0) + ylab("observed frequency") + xlab("c value") + ggtitle(paste0("best add model fit scaled\n", dosage, " dose")) + theme_classic()
  
  num.iterations.for.pval.calc <- 1000
  saved.cval.results <- list()
  for (ii in 1:num.iterations.for.pval.calc) {
    frac.add <- 1
    sim.cvals <- create.simulated.set.of.cvals(siUpregGenes, frac.add, dosage, n.iterations.per.input.observation = 1, seed = ii)
    saved.cval.results[[ii]] <- sim.cvals
    
    if ((ii %% 100) == 0) {
      print(sprintf("%d iterations complete", ii))
    }
  }
  
  # use larger bin size and do "sliding window" across data to generate p-value plot
  bin.step.size   <-  0.25
  bin.radius      <- bin.step.size / 2
  
  non.edge.pvals <- c()
  non.edge.bin.midpoints <- c()
  
  for (bin.adjustment in c(0, 0.25, 0.5, 0.75)) {
    bin.midpoints <- seq(bin.leftmost + bin.step.size, bin.rightmost, by = bin.step.size) - bin.radius + bin.adjustment * bin.step.size
    measured.cvalTable <- create.histogram.y.values(hist.values.measured, bin.midpoints)
    
    num.bins <- length(bin.midpoints)
    # n simulations x p bins
    histbin.matrix <- matrix(nrow = num.iterations.for.pval.calc, ncol = num.bins)
    # turn this into a matrix using the bin midpoints and histogram function
    for (jj in 1:num.iterations.for.pval.calc) {
      this.cval.result.vec   <- saved.cval.results[[jj]]
      this.binned.hist.data  <- create.histogram.y.values(this.cval.result.vec, bin.midpoints)
      # scale binned histogram data
      scaled.binned.hist.data <- this.binned.hist.data
      histbin.matrix[jj,] <- round(scaled.binned.hist.data)
    }
    
    #qplot(histbin.matrix[,10], binwidth = 1)
    
    p.values <- c()
    for (kk in 2:(num.bins - 1)) {
      mean.est <- mean(histbin.matrix[,kk])
      var.est <- var(histbin.matrix[,kk])
      # p for binomial distribution
      p.est <- mean.est / num.iterations.for.pval.calc
      #var.est <- var(histbin.matrix[,kk])
      measured.val.scaled <- round(measured.cvalTable[kk] * (1 / sim.additive.scaling.factor))
      #p.val.est <- 1 - pbinom(scaled.measured.val, num.iterations.for.pval.calc, p.est) 
      p.val.est <- min(ppois(measured.val.scaled, mean.est, lower.tail = F, log.p = T), 
                       ppois(measured.val.scaled, mean.est, lower.tail = T, log.p = T))
      
      p.values <- c(p.values, p.val.est)
    }
    
    non.edge.pvals <- c(non.edge.pvals, p.values)
    non.edge.bin.midpoints <- c(non.edge.bin.midpoints, bin.midpoints[2:(num.bins - 1)])
  }
  
  # see if mean and variance are about equal like in poisson
  bin_measured_means     <- c()
  bin_measured_variances <- c()
  for (ii in 1:ncol(histbin.matrix)) {
    print(sprintf("bin midpoint = %f, mean = %f, variance = %f", bin.midpoints[ii], mean(histbin.matrix[, ii]), var(histbin.matrix[, ii])))
    bin_measured_means <- c(bin_measured_means, mean(histbin.matrix[, ii]))
    bin_measured_variances <- c(bin_measured_variances, var(histbin.matrix[, ii]))
  }
  tib.for.meanvarplot <- tibble(bin.midpoints, bin_measured_means, bin_measured_variances)
  
  p1 <- ggplot(tib.for.meanvarplot) + geom_line(mapping = aes(x = bin.midpoints, y=  bin_measured_means), color = "blue") +
    geom_line(mapping = aes(x = bin.midpoints, y = bin_measured_variances), color = "green") + ggtitle(dosage)
  
  p2 <- qplot(non.edge.bin.midpoints, non.edge.pvals) + xlab("c value bin midpoint") + ylab("log p value estimate") + 
    geom_vline(xintercept = 0) + geom_vline(xintercept = 1) + theme_classic(base_size = 16) + ggtitle(paste0("log p value for finding the measured number of counts\nat each bin in an additive model,\n", dosage," dose, poisson model"))
  
  print(p1)
  print(p2)
  ggsave(paste0(output.folder, '/binned_log_p_values_', dosage, '_dose.svg'),       plot = p2, width = single.plot.width, height = single.plot.height)
  ggsave(paste0(output.folder, '/mean_var_plot_across_bins_', dosage, '_dose.svg'), plot = p1, width = single.plot.width, height = single.plot.height)
  
  p.composite2 <- p9 + p8 + p2
  ggsave(paste0(output.folder, '/composite_pval_plot_', dosage, '_dose.svg'), plot = p.composite2, width = single.plot.width * 3, height = single.plot.height * plot.vertical.shrinkage.factor)
  
  saved.plot.list.logpvals[[length(saved.plot.list.logpvals) + 1]] <- p1
  saved.plot.list.logpvals[[length(saved.plot.list.logpvals) + 1]] <- p2
}



