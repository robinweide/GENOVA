plot.chrom.comparison <- function( mat, color.fun = NULL, z.max = NULL){
	
	if(is.null(color.fun)){
		color.fun <- colorRampPalette(c("blue","red", "white"))
	}

	log.mat <- log10(mat$rawCounts/mat$normMat)
	#get z.max from the data if not defined
	if(is.null(z.max)){
		z.max <- max(log.mat[lower.tri(log.mat)])
	}	
	log.mat[log.mat > z.max] <- z.max
  #plot the matrix and use the chromosome names as row ids
	image(log.mat, axes=F, col=color.fun(1000), lwd=2)
	axis(1, at=seq(0,1, len=nrow(mat$rawCounts)), lab=rownames(mat$normMat), las=2, lwd=2)
	axis(2, at=seq(0,1, len=nrow(mat$rawCounts)), lab=rownames(mat$normMat), las=2, lwd=2)
	box(lwd=2)
}	



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

	plot.chrom.comparison(mat)
	invisible(mat)
}
