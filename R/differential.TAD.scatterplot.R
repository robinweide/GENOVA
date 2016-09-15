#' differential.TAD.scatterplot 
#'
#' Create a scatter plot of the inter and intra TAD contact frequencies.
#'
#' @param exp1 Output of intra.inter.TAD.contacts for experiment 1.
#' @param exp2 Output of intra.inter.TAD.contacts for experiment 2.
#' @param color.fun A color-function like rainbow().
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @return A scatterplot.
#' @import data.table
#' @export
#'
differential.TAD.scatterplot <- function( exp1, exp2, color.fun=NULL, xlab=NULL, ylab=NULL, ... ){
  #combine the two experimental conditions
  comb.exp <- merge(exp1, exp2, by=c("x","y"))
  max.neighbor <- max(comb.exp[,2]-comb.exp[,1])
  #if no user-defined color function is presented take jet.colors
  if(is.null(color.fun)){
    color.fun <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  if(is.null(xlab)){
    xlab="Experiment 1"
  }
  if(is.null(ylab)){
    ylab="Experiment 2"
  }       
  col.vec <- color.fun(max.neighbor + 1)
  tad.dist <- comb.exp[,2]-comb.exp[,1] + 1 #add one so that we can get the color from the col.vec
  plot(comb.exp[,3], comb.exp[,4], pch='.', col=col.vec[tad.dist], log='xy', xlab=xlab, ylab=ylab, ...)
  invisible(comb.exp)
}
