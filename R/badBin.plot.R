#' Print RCP_base versus RCP-test
#'
#' Used witin badBin.find
#'
#' @param hicMat The Hi-C Ice-data.table
#' @param res The resolution
#' @param bin The bin to plot
#' @return A list containing a vector with the bad bin-id's (`badBins`). Optionally, it can output a bed-file (`bed`).
badBin.plot <- function(hicMat, res, bin, RCP_base = RCP_base){
  binStop <- (1e6 / res) + bin
  x <- rep(bin, (1e6 / res))
  y <- seq(bin,binStop-1)
  vals <- hicMat[list(x,y)]
  vals[is.na(vals$V3)]$V3   <- 0
  plot(RCP_base , col = 'blue', type='l', xlab = "RCP_test", ylab = 'RCP', main = paste0("Results for bin ",bin))
  lines(vals$V3, col = 'red')
}
