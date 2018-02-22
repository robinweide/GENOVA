#' differential.TAD.scatterplot
#'
#' Create a scatter plot of the inter and intra TAD contact frequencies.
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp1 Output of intra.inter.TAD.contacts for experiment 1.
#' @param exp2 Output of intra.inter.TAD.contacts for experiment 2.
#' @param line Add a diagonal line?
#' @param allData Set to FALSE to get a zoomed in view.
#' @param color.fun A color-function like rainbow().
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @return A scatterplot.
#' @import data.table
#' @examples
#' # get scores for WT and WAPL data
#' TAD_N_WT   <- intra.inter.TAD.contacts(TAD = WT_TADs,
#'                                        max.neighbor = 10,
#'                                        exp = WT_10kb)
#' TAD_N_WAPL <- intra.inter.TAD.contacts(TAD = WT_TADs,
#'                                        max.neighbor = 10,
#'                                        exp = WAPL_10kb)
#'
#' # plot a TAD+N scatterplot
#' differential.TAD.scatterplot(exp1 = TAD_N_WT, # x
#'                              exp2 = TAD_N_WAPL, # y
#'                              allData = T)
#' @export
#'
differential.TAD.scatterplot <- function( exp1, exp2, line = T, allData = T,
                                          color.fun=NULL, xlab=NULL,
                                          ylab=NULL, ... ){
  #combine the two experimental conditions
  comb.exp <- merge(exp1$hic, exp2$hic, by=c("x","y"))
  max.neighbor <- max(comb.exp[,2]-comb.exp[,1])

  if(!allData){ # remove outliers. Also makes the plot square
    dat = c(comb.exp[,3], comb.exp[,4])
    Q = quantile(dat, c(.001,.999))

    comb.exp[comb.exp[,3] < Q[1],3] = Q[1]
    comb.exp[comb.exp[,3] > Q[2],3] = Q[2]

    comb.exp[comb.exp[,4] < Q[1],4] = Q[1]
    comb.exp[comb.exp[,4] > Q[2],4] = Q[2]
  }

  #if no user-defined color function is presented take jet.colors
  if(is.null(color.fun)){
    PAL = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
            "yellow", "#FF7F00", "red", "#7F0000")
    color.fun <- colorRampPalette(PAL)
  }
  if(is.null(xlab)){
    xlab="Experiment 1"
  }
  if(is.null(ylab)){
    ylab="Experiment 2"
  }
  col.vec <- color.fun(max.neighbor + 1)
  #add one so that we can get the color from the col.vec
  tad.dist <- comb.exp[,2]-comb.exp[,1] + 1
  plot(comb.exp[,3], comb.exp[,4], pch='.', col=col.vec[tad.dist],
       log='xy', xlab=xlab, ylab=ylab, ...)
  if(line){
    abline(a = 0, b = 1, lty = 3)
  }
  invisible(comb.exp)
}
