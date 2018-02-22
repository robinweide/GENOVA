#' plot.chrom.comparison
#'
#' Plot a genome-wide overview-matrix.
#'
#' @param mat Output of chromosomeMatrix
#' @param color.fun Optional color-function
#' @param z.max Adjust the maximum value of the color-scale
#' @return A plot of the matrix
chrom.comparison.plot <- function(mat, color.fun = NULL, z.max = NULL){

  if(is.null(color.fun)){
    color.fun <- colorRampPalette(c('white', '#f5a623', '#d0021b', 'black'))
  }

  log.mat <- log10(mat$rawCounts / mat$normMat)

  #get z.max from the data if not defined
  if(is.null(z.max)){
    z.max <- max(log.mat[lower.tri(log.mat)])
    warning("No cut.off was given: using 99.5% percentile: ",
            round(z.max),
            ".")
  }

  log.mat[log.mat > z.max] <- z.max
  #plot the matrix and use the chromosome names as row ids
  image(log.mat,
        axes = F,
        col = color.fun(1000),
        lwd = 2)
  axis(1,
       at  = seq(0,1, len = nrow(mat$rawCounts)),
       lab = rownames(mat$normMat),
       las = 2,
       lwd = 2)
  axis(2,
       at  = seq(0,1, len = nrow(mat$rawCounts)),
       lab = rownames(mat$normMat),
       las = 2,
       lwd = 2)

  box(lwd = 2)
}
