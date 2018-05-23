#' Plot pyramid-style HiC-regions.
#'
#' Extracts a sub-matrix from given coordinates and plots them pyramid-style (i.e. the Hi-C diagonal becomes horizontal).
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chrom Chromosome.
#' @param start Start position in bp.
#' @param end End position in bp.
#' @param window The window-size (default: 10e3).
#' @param q.val Quantile cutoff
#' @param cut.off Score cutoff
#' @param extraPadding Extra room needed for the ends of the piramid.
#' @param ylim The y-start and -stop: zero is at the diagonal,
#' @return A plot.
#' @import data.table
frequency.from.matrix <- function( experiment, chrom, start, end, window = 10e3, q.val=0.95, cut.off=NULL, extraPadding = 1.75, ylim = c(0, 250), shinyAxis = T, ... ){
  # Check size: if less than 50 bins: makes no sense...
  resolution <- experiment$RES
  mat <- select.subset(experiment, chrom, start- 2e7, end+ 2e7)
  pos <- which(mat$z >0 , arr.ind=T)
  mat.start <- min(mat$x)
  cnt <- as.vector(mat$z[mat$z > 0])
  cnt <- cnt[pos[,2]-pos[,1] > 0]
  pos <- pos[pos[,2]-pos[,1] > 0,]
  if(is.null(cut.off)){
    cut.off = max(quantile(cnt, .995))
    warning("No cut.off was given: using 99.5% percentile: ", round(cut.off), ".")
  }
  cnt[cnt > cut.off ] <- cut.off
  cnt <- cnt/cut.off
  num.rect <- length(unique(pos[,1]))
  cex = (400/num.rect)*2.5
  if(shinyAxis){
    plot( mat.start+window*(pos[,1]+pos[,2])/2, pos[,2]-pos[,1], pch=18, col=rgb(1, 1-cnt, 1-cnt), xlab="", ylab="", ylim = ylim,cex=cex, xlim=c(start,end), axes=F, ...); box(lwd=1); axis(1, at=seq(0,3e9, by=5e5), lab=NA, lwd=1)
  } else {
    plot( mat.start+window*(pos[,1]+pos[,2])/2, pos[,2]-pos[,1], pch=18, col=rgb(1, 1-cnt, 1-cnt), xlab="", ylab="", ylim = ylim,cex=cex, xlim=c(start,end), ...)
  }
}
