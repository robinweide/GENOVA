#' differential.TAD.dotplot 
#'
#' Create a dot plot showing the log-ratio of the intra and inter TAD contact frequencies for different TAD distances.
#'
#' @param exp1 Output of intra.inter.TAD.contacts for experiment 1.
#' @param exp2 Output of intra.inter.TAD.contacts for experiment 2.
#' @param color.fun A color-function like rainbow().
#' @return A dotplot
#' @import data.table
#' @export
#'
differential.TAD.dotplot <- function( exp1, exp2, color.fun=NULL, ... ){
  comb.exp <- merge(exp1, exp2, by=c("x","y"))
  max.neighbor <- max(comb.exp[,2]-comb.exp[,1])
  #if no user-defined color function is presented take jet.colors
  if(is.null(color.fun)){
    color.fun <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }
  col.vec <- color.fun(max.neighbor + 1)
  tad.dist <- comb.exp[,2]-comb.exp[,1] + 1 #add one so that we can get the color from the col.vec
  log.ratio <- log2(comb.exp[,4]/comb.exp[,3])
  plot( tad.dist + runif(nrow(comb.exp), -0.4, 0.4), log.ratio, pch='.', col=col.vec[tad.dist], ylab="Log2-ratio contact frequency", xlab="TAD distance", xaxt='n', ...)
  axis(1, at=1:(max.neighbor+1), lab=paste0("n + ", 0:max.neighbor), las=2)
  invisible(comb.exp)
}
