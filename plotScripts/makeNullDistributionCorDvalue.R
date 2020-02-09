library(tidyverse)
library(here)
library(VGAM)  # required for rfoldnorm function

cmdargs = commandArgs(trailingOnly=TRUE)
if (length(cmdargs) == 0) {
  dist.for.peaks.vs.genes  <- "peaks" 
  add.vs.mult.null.model   <- "additive" # choose "additive" or "multiplicative" or "mixture 
  add.mult.mixture.frac.add <- 0.50
  
  # input files -- use upregulated peaks or genes
} else {
  # siUpregGenes       <- read_tsv(cmdargs[1])
  # siUpregJoinedPeaks <- read_tsv(cmdargs[2])
  # outputloc.prefix   <- cmdargs[3]
  # selected.peak.category.arg <- cmdargs[4] 
}