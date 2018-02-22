#' differential.TAD.dotplot
#'
#' Create a dot plot showing the log-ratio of the intra and
#' inter TAD contact frequencies for different TAD distances.
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp1 Output of intra.inter.TAD.contacts for experiment 1.
#' @param exp2 Output of intra.inter.TAD.contacts for experiment 2.
#' @param color.fun A color-function like rainbow().
#' @param yRange Set the y-axis limits, a two-number vector,
#' otherwise use a quantile-based approach to find yRange
#' @param pch Which plot-characters?
#' @param title Add a title to the plot.
#' @param cex How big do you want the points?
#' @return A dotplot
#' @import data.table
#' @examples
#' # get scores for WT and WAPL data
#' TAD_N_WT   <- intra.inter.TAD.contacts(TAD = WT_TADs,
#'                                        max.neighbor = 10,
#'                                        exp = Hap1_WT_10kb)
#' TAD_N_WAPL <- intra.inter.TAD.contacts(TAD = WT_TADs,
#'                                        max.neighbor = 10,
#'                                        exp = Hap1_WAPL_10kb)
#'
#' # plot a TAD+N dotplot
#' differential.TAD.dotplot(exp1 = TAD_N_WT, # denominator
#'                          exp2 = TAD_N_WAPL) # numerator
#' @export
differential.TAD.dotplot <- function( exp1, exp2, color.fun=NULL, yRange =NULL,
                                      pch = '.', title = NULL, cex = 1, ... ){
  comb.exp <- merge(exp1$hic, exp2$hic, by=c("x","y"))
  max.neighbor <- max(comb.exp[,2]-comb.exp[,1])
  #if no user-defined color function is presented take jet.colors
  if(is.null(color.fun)){
    PAL = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
            "yellow", "#FF7F00", "red", "#7F0000")
    color.fun <- colorRampPalette(PAL)
  }
  col.vec <- color.fun(max.neighbor + 1)
  #add one so that we can get the color from the col.vec
  tad.dist <- comb.exp[,2]-comb.exp[,1] + 1
  log.ratio <- log2(comb.exp[,4]/comb.exp[,3])
  yLab = paste0('log2 contact frequency ( ',
                exp2$sampleName,
                " / ",
                exp1$sampleName,
                " )")

  Q = quantile(log.ratio, c(0.0015,0.9985))
  Q = max(abs(Q))
  if(is.null(yRange)){
    plot( tad.dist + runif(nrow(comb.exp), -0.4, 0.4),
          log.ratio,
          pch=pch,
          col=col.vec[tad.dist],
          cex = cex,
          ylab=yLab,
          xlab="TAD distance",
          xaxt='n',
          ylim = c(-1*Q, Q),
          main = title, ...)
    axis(1,
         at=1:(max.neighbor+1),
         lab=paste0("n + ", 0:max.neighbor), las=2)
  } else {
    plot( tad.dist + runif(nrow(comb.exp), -0.4, 0.4),
          log.ratio,
          pch=pch,
          col=col.vec[tad.dist],
          cex = cex,
          ylab=yLab,
          xlab="TAD distance",
          xaxt='n',
          ylim = yRange,
          main = title, ...)
    axis(1, at=1:(max.neighbor+1), lab=paste0("n + ", 0:max.neighbor), las=2)
  }

  invisible(comb.exp)
}
