#' chromosomeMatrix
#'
#' Plot a genome-wide overview-matrix.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param color.fun Optional color-function
#' @param z.max Adjust the maximum value of the color-scale
#' @param remove A vector of (partial) names of chromosomes to leave out.
#' @note Please use only low-resolution matrices with this (e.g. 1Mb resolution).
#' @return A plot of the chromosome matrix and a assignable matrix
#' @export
chromosomeMatrix <- function( exp, color.fun = NULL, z.max = NULL, remove = NULL ){
  exp$ABS[,1] <- as.character(exp$ABS[,1])
  chrom <- c(exp$ABS[1,1],exp$ABS[which(head(exp$ABS[,1],-1) != tail(exp$ABS[,1],-1))+1,1])

  chr1 <- factor(exp$ABS[exp$ICE$V1,1], levels=chrom)
  chr2 <- factor(exp$ABS[exp$ICE$V2,1], levels=chrom)
  #calculate the number of interactions per/between chromosome(s)
  chrom.count <- aggregate(V3 ~ chr1 + chr2, data=exp$ICE, sum)
  chrom.count.inv <- chrom.count[chrom.count[,1]!=chrom.count[,2],c(2,1,3)]
  names(chrom.count.inv) <- names(chrom.count)
  chrom.count <- rbind(chrom.count, chrom.count.inv)

  #make sure all the chromosomes are covered
  all.combinations <- cbind(rep(chrom,each=length(chrom)), rep(chrom, length(chrom)))
  chrom.count <- merge(all.combinations, chrom.count, by=c(1,2), all.x=T, sort=F)
  chrom.count[is.na(chrom.count[,3]),3] <- 0
  #for some reason two merge operation are necessary to get
  #the correct order
  chrom.count <- merge(all.combinations, chrom.count, by=c(1,2), all.x=T, sort=F)

  #create a normalization matrix
  chrom.n <- table(factor(exp$ABS[,1], levels=chrom))
  norm.mat <- outer(chrom.n, chrom.n, "*")
  chrom.mat <- matrix(chrom.count[,3], ncol=length(chrom), byrow=T)

  if(!is.null(remove)){
    remove.vec <- NULL
    for(remove_item in remove){
      hitsToRemove <- grep(remove_item, rownames(norm.mat))
      remove.vec <- c(remove.vec, hitsToRemove)
    }
    chrom.mat <- chrom.mat[-remove.vec,-remove.vec]
    norm.mat <- norm.mat[-remove.vec,-remove.vec]
  }


  mat <- list(rawCounts=chrom.mat , normMat = norm.mat)

  chrom.comparison.plot(mat, color.fun = NULL, z.max = NULL)
  invisible(mat)
}
