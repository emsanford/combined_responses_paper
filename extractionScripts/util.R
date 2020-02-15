## functions used by multiple other extraction scripts

assignValuesToHistBin <- function(values, bin_midpoints, bin.radius) {
  n.vals <- length(values)
  outputVec <- c()
  for (ii in 1:n.vals) {
    this.val <- values[ii]
    distvec <- abs(this.val - bin_midpoints)
    lowest.bin.distance <- min(distvec)
    bin.index <- which(distvec == lowest.bin.distance)[1]
    outputVec <- c(outputVec, bin_midpoints[bin.index])
  }
  return(outputVec)
}

# optional todo: if we end up using the colored stacked histograms, we can specify the colors here  
makeHistogramOfValues <- function(data.vector, categories.vector, xlim_lower, xlim_upper,
                                  bin.step.size, plot.title, 
                                  xlabel = "", ylabel = "", color.by.category = T) {
  
  
  bin.radius      <- bin.step.size / 2
  bin.midpoints   <- seq(xlim_lower + bin.step.size, xlim_upper, by = bin.step.size) - bin.radius

  bin.values <- assignValuesToHistBin(hist.values, bin.midpoints, bin.radius)
  
  stackedBarHistTib <- tibble(intConstantHhistBin = bin.values, intCategory = categories.vector)
  
  if (color.by.category) {
    stackedBarHist <- ggplot(stackedBarHistTib, aes(x = intConstantHhistBin, fill = intCategory))
  } else {
    stackedBarHist <- ggplot(stackedBarHistTib, aes(x = intConstantHhistBin))
  }
  
  stackedBarHist <- stackedBarHist +
    geom_bar(stat="count", width = bin.radius * 2 * .90) + 
    theme_minimal(base_size = 12) + 
    # ylim(0, 45) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(paste0(plot.title, "\nleft-val ", table(bin.values)[1], ", right-val ", table(bin.values)[length(table(bin.values))])) +
    geom_vline(xintercept = 0) + geom_vline(xintercept = 1) +
    xlim(min(bin.midpoints) - bin.radius, max(bin.midpoints) + bin.radius)
  
  return(stackedBarHist)
}