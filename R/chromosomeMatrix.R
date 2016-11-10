#' chromosomeMatrix
#'
#' Plot a genome-wide overview-matrix.
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param color.fun Optional color-function
#' @param z.max Adjust the maximum value of the color-scale
#' @note Please use only low-resolution matrices with this.
#' @return A plot of the chromosome matrix and a assignable matrix
#' @export
chromosomeMatrix <- function( exp, color.fun = NULL, z.max = NULL ){
	exp$ABS[,1] <- as.character(exp$ABS[,1])
  #get the chromosome ids (in the order they appear in the bed file)
	chrom <- c(exp$ABS[1,1],exp$ABS[which(head(exp$ABS[,1],-1) != tail(exp$ABS[,1],-1))+1,1])

	chr1 <- factor(exp$ABS[exp$ICE$V1,1], levels=chrom)
	chr2 <- factor(exp$ABS[exp$ICE$V2,1], levels=chrom)
	#calculate the number of interaction/chromosome
	chrom.count <- aggregate(V3 ~ chr1 + chr2, data=exp$ICE, sum)
	chrom.count.inv <- chrom.count[chrom.count[,1]!=chrom.count[,2],c(2,1,3)]
	names(chrom.count.inv) <- names(chrom.count)
	chrom.count <- rbind(chrom.count, chrom.count.inv)
	chrom.count <- chrom.count[order(chrom.count[,1], chrom.count[,2]),]

	#create a normalization matrix
	chrom.n <- table(factor(exp$ABS[,1], levels=chrom))
	norm.mat <- outer(chrom.n, chrom.n, "*")

	chrom.mat <- matrix(chrom.count[,3], ncol=length(chrom), byrow=T)

	mat <- list(rawCounts=chrom.mat , normMat = norm.mat)

	plot.chrom.comparison(mat, color.fun=color.fun, z.max=z.max)
	invisible(mat)
}
