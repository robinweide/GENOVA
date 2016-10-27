#' Get matrix from a BED-like entry.
#'
#' Extracts a symmetric matrix around the diagonal from *start* to *stop* in chromosome *chrom*.
#'
#' @param r1 A HiC-pro matrix.
#' @param bed A HiC-pro index.
#' @param chrom Chromosome.
#' @param start Start position in bp.
#' @param end End position in bp.
#' @return A data.table with normalised counts.
select.subset <- function(r1, chrom, start, end, bed  ){
  sel <- bed[bed[,1]==chrom & bed[,2] >= start & bed[,2] <= end,4]
  start.i <- min(sel); end.i <- max(sel)
  #create matrix
  r1.s <- select.sub(r1, start.i, end.i )
  r2.s <- r1.s
  sub.mat <- matrix(0, ncol=length(start.i:end.i), nrow=length(start.i:end.i))
  #fill the top triangle of the matrix
  index <- cbind(r1.s$V1-start.i+1,r1.s$V2-start.i+1)
  sub.mat[index] <- r1.s$V3
  #fill the bottom triangle of the matrix
  index <- cbind(r2.s$V2-start.i+1,r2.s$V1-start.i+1)
  sub.mat[index] <- r2.s$V3
  #take the middle of the window as the position int the genome
  pos <- (bed[bed[,4]%in%(start.i:end.i),2]+bed[bed[,4]%in%(start.i:end.i),3])/2
  list(x=pos, y=pos, z=sub.mat)
}
