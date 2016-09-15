#' PESCAn
#'
#' From a ChIP bed file, a HiCpro matrix and a HiCpro bed file calculate a PE-scan like data structure. Run this per chromosome.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param bed A bed-dataframe.
#' @param minDist The minimal distance
#' @param add Add constant value to bed-start and -end.
#' @param size The amount of Hi-C bins to take into account (i.e. a score of 21 yield an output with 10 Hi-C bins up- and downstream of the anchor). 
#' @return A score-matrix.
#' @import data.table
#' @export
#' 
PESCAn <- function( experiment, bed, minDist = 5e6, size = 500e3, add = 0 ){
  #sorting the bed file is essential for the analysis
  bed <- bed[order(bed[,1],bed[,2]),]
  count = 0
  for( chr in unique(bed[,1]) ){
    cat("Analyzing ", chr, "\r")
    pe.res <- cov2d(experiment, bed[bed[,1]==chr,], minDist, size, add)
    if(exists("score.mat")){
      score.mat <- score.mat + pe.res$score
      count = count + pe.res$count
    }else{
      score.mat <- pe.res$score
      count = pe.res$count
    }
  }
  score.mat/count
}  
