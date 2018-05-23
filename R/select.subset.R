#' Get matrix from a BED-like entry.
#'
#' Extracts a symmetric matrix around the diagonal from start to stop in chromosome chrom.
#'
#' @param exp The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chrom Chromosome.
#' @param start Start position in bp.
#' @param end End position in bp.
#' @examples
#' # get the TP53-locus in an experiment mapped to hg19
#' ss_out <- select.subset(experiment = WT, chrom = 'chr17', start = 7.5e6, end = 7.6e6)
#'
#' # plot the region
#' image(ss_out$z)
#' @return A list with the X and Y locations and a matrix (Z) containing the contacts.
#' @export
select.subset <- function(exp, chrom, start, end){
  sel <- exp$ABS[exp$ABS[,1]==chrom & exp$ABS[,2] >= start & exp$ABS[,2] <= end,4]
  start.i <- min(sel); end.i <- max(sel)
  #create matrix
  r1.s <- select.sub(exp$ICE, start.i, end.i )
  r2.s <- r1.s
  sub.mat <- matrix(0, ncol=length(start.i:end.i), nrow=length(start.i:end.i))
  #fill the top triangle of the matrix
  index <- cbind(r1.s$V1-start.i+1,r1.s$V2-start.i+1)
  sub.mat[index] <- r1.s$V3
  #fill the bottom triangle of the matrix
  index <- cbind(r2.s$V2-start.i+1,r2.s$V1-start.i+1)
  sub.mat[index] <- r2.s$V3
  #take the middle of the window as the position int the genome
  pos <- (exp$ABS[exp$ABS[,4]%in%(start.i:end.i),2]+exp$ABS[exp$ABS[,4]%in%(start.i:end.i),3])/2
  list(x=pos, y=pos, z=sub.mat)
}
