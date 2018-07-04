#' PESCAn_covert
#'
#' From a ChIP bed file, a HiCpro matrix and a HiCpro bed file calculate a PE-scan like data structure. Run this per chromosome.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minComparables The minimal amount of bed-entries for a given chromosome. If this threshold is not reached, PE-SCAn will skip this chromosome.
#' @param minDist The minimal distance
#' @param maxDist The maximal distance
#' @param rmOutlier to outlier-correction
#' @param outlierCutOff outlierCutOff
#' @param add Add constant value to bed-start and -end.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor).
#' @return A score-matrix.
#' @import data.table
#' @examples
#' # Run PE-SCAn on a bed of super-enhancers, using WT Hi-C data.
#' WT_PE_OUT = PESCAn(exp = WT_40kb, bed = superEnhancers)
#'
#' # Plot using visualise.PESCAn.ggplot
#' visualise.PESCAn.ggplot(PESCAnlist = list(WT = WT_PE_OUT),
#'                         resolution = 40e3,
#'                         smooth = F)
#'
#' # Plot using persp
#' persp(SE_vs_WT/SE_vs_WT_perm, phi = 30, theta = 30, col = 'skyblue')
#'
PESCAn_covert <- function( experiment, bed, minComparables = 10, rmOutlier = F,
                           minDist = 5e6, maxDist = Inf,
                           size = 500e3, add = 0 , outlierCutOff = 0.995,
                           verbose = T){
  #sorting the bed file is essential for the analysis
  bed <- bed[order(bed[,1],bed[,2]),]
  count = 0

  # there could be circumstances where there is only one bed-entry for a specific chromosome!
  chromsomesToLookAt <- names(which(table(bed[,1]) > 1))
  for( chr in chromsomesToLookAt ){
    if(verbose){message("Analyzing ", chr)}

    BED <- bed[bed[,1]==chr,]
    if(nrow(BED) < minComparables){
      next()
      }
    pe.res <- cov2d(experiment = experiment, bed = BED, minDist = minDist,
                    maxDist = maxDist,
                    size = size, add =  add, rmOutlier = rmOutlier,
                    outlierCutOff = outlierCutOff, verbose = verbose)
    if(is.null(pe.res)){
      next()
    }
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
#' From a ChIP-peaks.BED dataframe, calculate the all-vs-all Hi-C contacts.
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minDist The minimal distance
#' @param maxDist The maximal distance
#' @param shift Set to X bp for circular permutation. Set to zero for just getting the signal-matrix.
#' @param size Size in bp of window.
#' @param minComparables Minimum amount of BED-entries per chromosome.
#' @param rmOutlier Perform outlier-correction
#' @param outlierCutOff The severity of outliers. We compute the [outlierCutOff]
#' percentile per pixel and set values bigger than that to this value.
#' @return An list with a O/E score-matrix (if shift != 0), otherwise an observed
#' score-matrix, the underlying signal and background-matrices and the shift used.
#' @import data.table
#' @examples
#' # Run PE-SCAn on a bed of super-enhancers,
#' using WT Hi-C data and a circular permutation of 1Mb
#' WT_PE_OUT = PESCAn(exp = WT_40kb,
#'                    bed = superEnhancers,
#'                    shift = 1e6)
#'
#' # Plot using visualise.PESCAn.ggplot
#' visualise.PESCAn.ggplot(PESCAnlist = list(WT = WT_PE_OUT),
#'                         resolution = 40e3,
#'                         smooth = F)
#'
#' # Plot using persp
#' RES = 40e3 # resolution of the Hi-C
#' persp(list(x = seq(-1*(RES*10),(RES*10), length.out = 21)/1e6, # x-ticks (MB)
#'            y = seq(-1*(RES*10),(RES*10), length.out = 21)/1e6, # y-ticks (MB)
#'            z = WT_PE_OUT)
#' @export
PESCAn = function(exp, bed, shift = 1e6, minComparables = 10, mindist = 5e+06, maxDist = Inf,
                  size = 4e+05, rmOutlier = F, outlierCutOff = 0.995,
                  verbose = T){

  # Get signal
  if(verbose){message("Computing observed of ", exp$NAME)}
  signal = PESCAn_covert(experiment = exp, minComparables = minComparables,
                         bed = bed, minDist = mindist,
                         size = size, rmOutlier = rmOutlier, maxDist = maxDist,
                         outlierCutOff = outlierCutOff, verbose = verbose)

  # Get O/E
  background = NULL
  OE = NULL
  if(shift == 0){
    OE = signal
  } else { # if a shift value is given
    # Get background
    if(verbose){message("Computing shifted of ", exp$NAME)}
    background = PESCAn_covert(experiment = exp, minComparables = minComparables,
                               bed = bed, add = shift,
                               minDist = mindist, size = size,
                               maxDist = maxDist, outlierCutOff = outlierCutOff,
                               verbose = verbose)
    medianBackground = median(background)
    OE = signal/medianBackground
  }

  # Return OE-matrix
  return(list('mat' = OE, 'signal' = OE, 'background' = background, "shift" = shift))

}

#' visualise.PESCAn.ggplot
#'
#' Plot the PE-SCAn-results and differentials.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param PESCAnlist A list of results from PESCAn. Should have the same [shift].
#' @param title Text to plot
#' @param Focus Which sample does need to be the to-compare sample?
#' @param zTop The min and max values for the first (observed or OE) row of plots.
#' @param zBottom The min and max values for the bottom (differential) row of plots.
#' @return A grid object, containing two ggplot-objects.
#' @examples
#' # Run PE-SCAn on a bed of super-enhancers,
#' using WT Hi-C data and a circular permutation of 1Mb
#' WT_PE_OUT = PESCAn(exp = WT_40kb,
#'                    bed = superEnhancers,
#'                    shift = 1e6)
#'
#' # Plot using visualise.PESCAn.ggplot
#' visualise.PESCAn.ggplot(PESCAnlist = list(WT = WT_PE_OUT),
#'                         resolution = 40e3,
#'                         smooth = F)
#' @export
visualise.PESCAn.ggplot = function (PESCAnlist, resolution, title = "PE-SCAn", zTop = NULL, zBottom = NULL, focus = 1, smooth = F, ...) {
  require(ggplot2)

  # check if all have OE or O:
  OE = F
  Olist = c()
  for(i in 1:length(PESCAnlist)){
    Olist = unique(c(Olist, PESCAnlist[[1]]$shift))
  }
  if(length(Olist) == 1){
    if(Olist == 0){
      OE = F
    } else {
      OE = T
    }
  } else {
    stop('All samples should have the same shift-value.')
  }


  size <- dim(as.data.frame(PESCAnlist[[1]]$mat))[1]
  size.banks <- (size - 1)/2
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
    firstMat <- as.data.frame(PESCAnlist[[i]]$mat)

    firstMat <- t(apply(firstMat, 2, rev))
    colnames(firstMat) <- 1:size
    rownames(firstMat) <- 1:size
    secondMat <- as.data.frame(PESCAnlist[[focus]]$mat)

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


  spectCol = c('#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd')


  plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value), interpolate = smooth) +
    ggplot2::facet_grid(. ~ sample) + ggplot2::coord_fixed() +
    GENOVA_THEME() +
    ggplot2::scale_x_continuous(breaks = c(tickPosDownstream,size.banks + 1,
                                           tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "3'",
                                           paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(tickPosDownstream,
                                           size.banks + 1, tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "5'",
                                           paste0(tickLabelDownstream, "kb"))) +
    ggplot2::labs(title = title, x = "", y = "", fill = "O/E")
  if(OE == T){
    plot1 <- plot1 +
      ggplot2::scale_fill_gradient2(trans = 'log2' , limits = z, midpoint = 0,
                                    low = "#2166ac", mid = "white", high = "#b2182b")
  } else {
    plot1 <- plot1 +
      ggplot2::scale_fill_gradientn(colours = spectCol)
  }


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
    ggplot2::scale_fill_gradient2(limits = z2, midpoint = 0, low = "#2166ac",
                                  mid = "white", high = "#b2182b") +
    GENOVA_THEME() +
    ggplot2::scale_x_continuous(breaks = c(tickPosDownstream,size.banks + 1,
                                           tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "3'",
                                           paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(tickPosDownstream,
                                           size.banks + 1, tickPosUpstream),
                                labels = c(paste0(tickLabelUpstream, "kb"), "5'",
                                           paste0(tickLabelDownstream, "kb"))) +
    ggplot2::labs(x = "", y = "", fill = "Difference")

  grid::grid.newpage()
  grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2),
                        size = "last"))
}
