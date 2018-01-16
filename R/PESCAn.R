#' PESCAn_covert
#'
#' From a ChIP bed file, a HiCpro matrix and a HiCpro bed file calculate a PE-scan like data structure. Run this per chromosome.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minComparables The minimal amount of bed-entries for a given chromosome. If this threshold is not reached, PE-SCAn will skip this chromosome.
#' @param minDist The minimal distance
#' @param add Add constant value to bed-start and -end.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor).
#' @return A score-matrix.
#' @import data.table
#' @examples
#' # Run PE-SCAn on a bed of super-enhancers, using WT Hi-C data.
#' SE <- read.delim(superEnhancers.bed, header = F)
#' SE_vs_WT <- PESCAn(experiment = WT, bed = SE, size = 500e3)
#'
#' # Do a circular perutation by adding 1Mb.
#' SE_vs_WT_perm <- PESCAn(experiment = WT, bed = SE, size = 500e3, add = 1e6)
#'
#' # Plot using, for example, persp
#' persp(SE_vs_WT/SE_vs_WT_perm, phi = 30, theta = 30, col = 'skyblue')
#'
PESCAn_covert <- function( experiment, bed, minComparables = 5, minDist = 5e6, size = 500e3, add = 0 ){
  #sorting the bed file is essential for the analysis
  bed <- bed[order(bed[,1],bed[,2]),]
  count = 0

  # there could be circumstances where there is only one bed-entry for a specific chromosome!
  chromsomesToLookAt <- names(which(table(bed[,1]) > 1))
  for( chr in chromsomesToLookAt ){
    message("Analyzing ", chr)
    BED <- bed[bed[,1]==chr,]
    if(nrow(BED) < minComparables){
      next()
      }
    pe.res <- cov2d(experiment, BED, minDist, size, add)
    if(exists("score.mat")){
      score.mat <- score.mat + pe.res$score
      count = count + pe.res$count
    }else{
      score.mat <- pe.res$score
      count = pe.res$count
    }
  }

  if(exists("score.mat")){
    return(score.mat/count)
  }
}

#' PE-SCAn
#'
#' From a ChIP bed file, a HiCpro matrix and a HiCpro bed file calculate a PE-scan like data structure. Run this per chromosome.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minDist The minimal distance
#' @param shift Set to X bp for circular permutation. Set to zero for just getting the signal-matrix.
#' @param size Size in bp of window.
#' @return A O/E score-matrix.
#' @import data.table
#' @examples
#' # Run PE-SCAn on a bed of super-enhancers, using WT Hi-C data.
#' SE <- read.delim(superEnhancers.bed, header = F)
#'
#' WT_PE = PESCAn_run(exp = WT_exp, bed = SE)
#' MUT_PE = PESCAn_run(exp = MUT_exp, bed = SE)
#'
#' list_PE = list("WT" = WT_PE, 'Mutant' = MUT_PE)
#' PESCAn_plot(list_PE, 4e4, title = 'PE-SCAn', zTop = c(1,1.5))
#'
#' @export
#'
PESCAn = function(exp, bed, shift = 1e6, mindist = 5e+06, size = 4e+05){

  # Get signal
  signal = suppressMessages(PESCAn_covert(experiment = exp, bed = bed, minDist = mindist, size = size))

  # Get background
  background = suppressMessages(PESCAn_covert(experiment = exp, bed = bed, add = shift, minDist = mindist, size = size))
  medianBackground = median(background)

  # Get O/E
  OE = NULL
  if(shift == 0){
    OE = signal
  } else {
    OE = signal/medianBackground
  }

  # Return OE-matrix
  return(OE)

}

#' Plot the PE-SCAn-results
#'
#' @param PESCAnlist A list of results from `PESCAn`.
#' @param title Text to plot
#' @param Focus Which sample does need to be the to-compare sample?
#' @param zTop The min and max colorscale-values for the first row of plots.
#' @param zBottom The min and max colorscale-values for the first row of plots.
#' @return A grid object, containing two ggplot-objects.
#' @export
visualise.PESCAn.ggplot = function (PESCAnlist, resolution, title = "PE-SCAn", zTop = NULL, zBottom = NULL, focus = 1, smooth = F, ...) {
  require(ggplot2)
  require(viridis)

  size <- dim(as.data.frame(PESCAnlist[[1]]))[1]
  size.banks <- (size - 1)/2
  # tickLabelDownstream <- as.character(1 * ((size.banks/2 *  resolution)/1000))
  # tickLabelUpstream <- as.character(-1 * ((size.banks/2 * resolution)/1000))

  allTicks = seq(-1*resolution*size.banks, resolution*size.banks, length.out = size)/1e3



  tickPosDownstream = median(1:size.banks)
  tickPosUpstream = median( (((size-1)/2)+2 ): size )
  tickLabelUpstream = mean(allTicks[c(tickPosDownstream-.5, tickPosDownstream+.5)])
  tickLabelDownstream = mean(allTicks[c(tickPosUpstream-.5, tickPosUpstream+.5)])



  abovePlots <- data.frame(Var1 = integer(), Var2 = integer(),
                           value = numeric(), sample = factor())
  belowPlots <- abovePlots
  list.len <- length(PESCAnlist)
  belownames <- vector()
  for (i in 1:list.len) {
    firstMat <- as.data.frame(PESCAnlist[[i]])

    firstMat <- t(apply(firstMat, 2, rev))
    colnames(firstMat) <- 1:size
    rownames(firstMat) <- 1:size
    secondMat <- as.data.frame(PESCAnlist[[focus]])

    secondMat <- t(apply(secondMat, 2, rev))

    colnames(secondMat) <- 1:size
    rownames(secondMat) <- 1:size
    a <- reshape2::melt(as.matrix(firstMat))
    a$sample <- factor(rep(names(PESCAnlist)[i], length(a[,1])))
    abovePlots <- rbind(abovePlots, a)
    e <- reshape2::melt(as.matrix(firstMat - secondMat))
    e$sample <- rep(paste0(names(PESCAnlist)[focus], " vs ",
                           names(PESCAnlist)[i]), length(e[, 1]))
    belowPlots <- rbind(belowPlots, e)
    belownames <- c(belownames, paste0(names(PESCAnlist)[focus], " vs ", names(PESCAnlist)[i]))
  }

  z <- NULL
  if (is.null(zTop)) {
    z <- c(quantile(na.exclude(abovePlots$value), 0.001),
           quantile(na.exclude(abovePlots$value), 0.999))
    z = (z[2] - 1)
    z = c(1-z , 1+z)

  } else {
    z <- zTop
  }

  #message(z)
  abovePlots$value[abovePlots$value > z[2]] <- z[2]
  abovePlots$value[abovePlots$value < z[1]] <- z[1]
  abovePlots$sample <- factor(abovePlots$sample, levels = c(levels(abovePlots$sample)[focus],
                                                            levels(abovePlots$sample)[!levels(abovePlots$sample) %in%
                                                                                        levels(abovePlots$sample)[focus]]))
  volgorde <- match(names(PESCAnlist), levels(abovePlots$sample))
  belowPlots$sample <- factor(belowPlots$sample, levels = belownames[volgorde])

  # higlassCol <- c("white", "#f5a623", "#d0021b", "black")
  # divcol = c('#7f3b08','#f7f7f7','#2d004b')
  # tercol = terrain.colors(100)
  spectCol = c('#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd')


  plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value), interpolate = smooth) +
    ggplot2::facet_grid(. ~ sample) + ggplot2::coord_fixed() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA",  colour = NA)) +
    ggplot2::scale_x_continuous(breaks = c(tickPosDownstream,size.banks + 1, tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(tickPosDownstream,
                                           size.banks + 1, tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "5'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::labs(title = title, x = "", y = "", fill = "O/E") +
    #viridis::scale_fill_viridis( limits = z)
    #ggplot2::scale_fill_gradientn(colours = spectCol, limits = z)
    ggplot2::scale_fill_gradient2(limits = z, midpoint = 1, low = "#2166ac", mid = "white", high = "#b2182b")



  z2 <- NULL
  if (is.null(zBottom)) {
    z2 <- c(quantile(na.exclude(belowPlots$value), 0.001),
            quantile(na.exclude(belowPlots$value), 0.999))

    z2 <- abs(z2[which.max(abs(z2))] )
    z2 <- unname(c(z2*-1, z2))

  } else {
    z2 <- zBottom
  }

  belowPlots$value[belowPlots$value > z2[2]] <- z2[2]
  belowPlots$value[belowPlots$value < z2[1]] <- z2[1]
  plot2 <- ggplot2::ggplot(belowPlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value), interpolate = smooth) +
    ggplot2::facet_grid(. ~ sample) + ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(limits = z2, midpoint = 0, low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA)) +
    ggplot2::scale_x_continuous(breaks = c(tickPosDownstream,size.banks + 1, tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(tickPosDownstream,
                                           size.banks + 1, tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "5'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::labs(x = "", y = "", fill = "Difference")
  grid::grid.newpage()
  grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2), size = "last"))
}
