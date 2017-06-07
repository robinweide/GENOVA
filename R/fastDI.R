#' Get DI-scores for a specific region
#'
#' Extracts a sub-matrix from given coordinates and calculates the Directionality-Index.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chrom Chromosome.
#' @param start Start position in bp.
#' @param end End position in bp.
#' @param max.i Maximal size of flanking regions to calculate DI.
#' @return A vector of DI-scores per Hi-C bin.
#' @import data.table
fastDI <- function(experiment, chrom, start, end, max.i = 100){
  resolution <- experiment$RES
  mat <- select.subset(experiment, chrom, start, end)$z

  data.vec <- as.vector(mat)
  iv <- as.vector(col(mat))
  dv <- as.vector(row(mat)-col(mat))
  #calculate the stuff for A
  data.A <- data.vec[dv > 0 & dv < max.i]
  i.A <- iv[dv > 0 & dv < max.i]
  d.A <- dv[dv > 0 & dv < max.i]

  #calculate the stuff for B
  data.B <- data.vec[dv < 0 & dv > -max.i]
  i.B <- iv[dv < 0 & dv > -max.i]
  d.B <- dv[dv < 0 & dv > -max.i]

  A <- tapply(data.A, i.A, sum)
  B <- tapply(data.B, i.B, sum)

  #expected
  E <- (A+B)/2
  term1 <- (B-A)/abs(B-A)
  term2 <- ((A-E)**2)/E + ((B-E)**2)/E
  di <- term1*term2
  di <- c(di,0)

  #remove the extremely low coverage scores
  cov.vec <- apply(mat,1, sum)
  cut <- median(cov.vec)-mad(cov.vec)*8
  sel <- unique(c(which(cov.vec < cut), which(cov.vec < cut)-1))
  sel <- sel[sel > 0]
  sel <- sel[sel < length(di)]
  di[sel] <- 0
  -di
}

# v <- fastDI(WT.sub$z, 100); plot(WT.sub$x, v, type='h',lwd = 1.5,ylim=c(-400,400), xlim=xlim, axes=F); axis(1, at=seq(0,3e9, by=5e5), lab=NA, lwd=3)
