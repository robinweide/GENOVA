#' plot.insulation.single
#'
#' Plot a single insulation
#'
#' @param exp1 The Hi-C experiment object: produced by construct.experiment().
#' @param chrom Chromosome
#' @param start Start position of the region of interest
#' @param end End position of the region of interest
#' @param cut.off The cut.off for Hi-C scores
#' @param window.size The sliding square size
#' @param local Per chromosome?
#' @return A plot
#' @export
plot.insulation.single <- function( exp1, chrom, start, end, cut.off=0, window.size = 21, local = T ){
	#create a plotting layout
	w = 6
	lay <- matrix(4, nrow=w, ncol=w)
	lay[2:w,2:w] <- 1; lay[1,] <- 2; lay[,1] <- 3; lay[1,1] <- 4
	layout(lay)
	par(mar=rep(1,4), xaxs="i", yaxs="i")
	#layout

	#get a matrix from the experiment
	mat1 <- select.subset( exp1$ICE, chrom, start, end, exp1$ABS)
	mat1$z[mat1$z > cut.off] <- cut.off
	wr <- colorRampPalette(c("white","red"))
	image( mat1, col=wr(256), axes=F, ylim=rev(range(mat1$x)) )
	box(lwd=2)
	size.region <- diff(range(mat1$x))
	if( size.region > 2e6){
		axis(2, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6, lwd=2, cex.axis=1.6)
	}else{
		lab <- seq(0,3e9, by=500e3)/1e6
		lab <- sprintf("%.1f",lab)
		axis(2, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=1.6)
	}

	extend <- exp1$RES*window.size
	ins.score <- insulation.score( exp1, window.size, chrom, start-extend, end+extend, local=local )
	#remove NA elements, because the can generate problems in the plotting
	#of lines etc.
	ins.score <- ins.score[is.finite(ins.score[,2]),]

	plot(ins.score[,1], ins.score[,2], xlim=range(mat1$x), type='l', axes=F)
	insulation.polygon(ins.score, rotate=F)
	axis(2)

	plot(ins.score[,2], ins.score[,1], ylim=rev(range(mat1$x)), xlim=rev(range(ins.score[,2])), type='l', axes=F)
	insulation.polygon(ins.score, rotate=T)
	axis(3)
}

#' plot.insulation.dual
#'
#' Plot insulation-scores of two experiments
#'
#' @param exp1 A Hi-C experiment object: produced by construct.experiment().
#' @param exp2 A Hi-C experiment object: produced by construct.experiment().
#' @param chrom Chromosome
#' @param start Start position of the region of interest
#' @param end End position of the region of interest
#' @param cut.off The cut.off for Hi-C scores
#' @param window.size The sliding square size
#' @param local Per chromosome?
#' @param delta plot the differential insulation scores
#' @return A plot
#' @export
plot.insulation.dual <- function( exp1, exp2, chrom, start, end, cut.off=0, window.size = 21, local = T, delta = F ){
	#make sure the resolutions are the same
	if(exp1$RES != exp2$RES){
		stop("The Hi-C matrices should have the same resolution")
	}

	if(!all(exp1$ABS[,4] == exp2$ABS[,4])){
		stop("Not all ICE indexes are the same. Are you these experiments were mapped to the same genome (build)?")
	}

	#create a plotting layout
	w = 6
	lay <- matrix(4, nrow=w, ncol=w)
	lay[2:w,2:w] <- 1; lay[1,] <- 2; lay[,1] <- 3; lay[1,1] <- 4
	layout(lay)
	par(mar=rep(1,4), xaxs="i", yaxs="i")
	#layout

	#get a matrix from the experiment
	mat1 <- select.subset( exp1$ICE, chrom, start, end, exp1$ABS)
	mat2 <- select.subset( exp2$ICE, chrom, start, end, exp2$ABS)

	mat1$z[lower.tri(mat1$z)] <- mat2$z[lower.tri(mat2$z)]

	mat1$z[mat1$z > cut.off] <- cut.off
	wr <- colorRampPalette(c("white","red"))
	image( mat1, col=wr(256), axes=F, ylim=rev(range(mat1$x)) )
	box(lwd=2)
	size.region <- diff(range(mat1$x))
	if( size.region > 2e6){
		axis(2, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6, lwd=2, cex.axis=1.6)
	}else{
		lab <- seq(0,3e9, by=500e3)/1e6
		lab <- sprintf("%.1f",lab)
		axis(2, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=1.6)
	}

	#plot the insulation score for the matrix
	extend <- exp1$RES*window.size
	#remove NA elements, because the can generate problems in the plotting
	#of lines etc.

	if(delta){ #plot the differential insulation scores
		ins.score1 <- insulation.score( exp1, window.size, chrom, start-extend, end+extend, local=local )
		ins.score2 <- insulation.score( exp2, window.size, chrom, start-extend, end+extend, local=local )
		delta.score <- ins.score1
		delta.score[,2] <- delta.score[,2] - ins.score2[,2]
		plot(delta.score[,1], delta.score[,2], xlim=range(mat1$x), type='l', axes=F)
		insulation.polygon(delta.score, rotate=F)
		axis(2)

		plot(delta.score[,2], delta.score[,1], ylim=rev(range(mat1$x)), xlim=rev(range(delta.score[,2])), type='l', axes=F)
		insulation.polygon(delta.score, rotate=T)
		axis(3)

	}else{ #plot the seperate insulation score value on the top and the side

		#plot the top side insulation
		ins.score <- insulation.score( exp2, window.size, chrom, start-extend, end+extend, local=local )
		ins.score <- ins.score[is.finite(ins.score[,2]),]
		plot(ins.score[,1], ins.score[,2], xlim=range(mat1$x), type='l', axes=F)
		insulation.polygon(ins.score, rotate=F)
		axis(2)

		#plot the left side insulation
		ins.score <- insulation.score( exp1, window.size, chrom, start-extend, end+extend, local=local )
		plot(ins.score[,2], ins.score[,1], ylim=rev(range(mat1$x)), xlim=rev(range(ins.score[,2])), type='l', axes=F)
		insulation.polygon(ins.score, rotate=T)
		axis(3)
	}
}


#' insulation.domainogram
#'
#' create a domainogram for the insulation scores
#'
#' @param exp1 A Hi-C experiment object: produced by construct.experiment().
#' @param chrom Chromosome
#' @param start Start position of the region of interest
#' @param end End position of the region of interest
#' @param cut.off The cut.off for Hi-C scores
#' @param window.size1 The minimal sliding square size
#' @param window.size2 The maximal sliding square size
#' @param step Thesliding square step-size
#' @param local Per chromosome?
#' @param zlim Zlims of the plot
#' @return A plot
#' @export
insulation.domainogram <- function( exp1, chrom, start, end, window.size1 = 5, window.size2 = 101, step = 2, local = T, zlim=c(-1, 0.5)){
	if(step %% 2){
		stop("Step size has to be even")
	}
	if(!window.size1 %% 2){
		stop("Window sizes has to be odd")
	}
	color.fun <- colorRampPalette(c("#f03b20", "#ffeda0", "white",  "#31a354") )
	color.vec <- color.fun(1000)
	plot( c(start,end), c(window.size1, window.size2), type='n', xaxs="i", yaxs='i', axes=F, ylab="window size", xlab="chromosomal position")
	window.vec <- seq(window.size1, window.size2, by=step)
	for( window.size in window.vec){
		extend <- exp1$RES*window.size
		ins.score <- insulation.score( exp1, window.size, chrom, start-extend, end+extend, local=local )
		#set the limits on the values of the insulation score
		ins.score[ins.score[,2] < zlim[1],2] <- zlim[1]
		ins.score[ins.score[,2] > zlim[2],2] <- zlim[2]
		x1 <- c(ins.score[1,1], (head(ins.score[,1],-1)+tail(ins.score[,1],-1))/2)
		x2 <- c((head(ins.score[,1],-1)+tail(ins.score[,1],-1))/2, tail(ins.score[,1],1))
		y1 <- window.size-step/2
		y2 <- window.size+step/2
		color.index <- (1000-1)*(ins.score[,2] - zlim[1])/(zlim[2]-zlim[1]) + 1
		rect(x1, y1, x2, y2, col=color.vec[color.index], border=NA)
	}
	at <- seq(1e6*floor(start/1e6), 1e6*floor(end/1e6), by=1e6)
	axis(1, at=at, lab=at/1e6, lwd=2)
	axis(2, lwd=2)
	box(lwd=2)
}


#create a domainogram for the insulation scores
delta.insulation.domainogram <- function( exp1, exp2, chrom, start, end, window.size1 = 5, window.size2 = 101, step = 2, local = T, zlim=c(-0.4, 0.4)){
	if(step %% 2){
		stop("Step size has to be even")
	}
	if(!window.size1 %% 2){
		stop("Window sizes has to be odd")
	}
	#set the differential color palette
	color.fun <- colorRampPalette(c("blue", "white",  "red") )
	color.res <- 1000
	color.vec <- color.fun(color.res)
	#create an empty plot for the domainogram
	plot( c(start,end), c(window.size1, window.size2), type='n', xaxs="i", yaxs='i', axes=F, ylab="window size", xlab="chromosomal position")
	window.vec <- seq(window.size1, window.size2, by=step)
	for( window.size in window.vec){
		extend <- exp1$RES*window.size
		delta.score <- insulation.score( exp1, window.size, chrom, start-extend, end+extend, local=local )
		score2      <- insulation.score( exp2, window.size, chrom, start-extend, end+extend, local=local )
		delta.score[,2] <- delta.score[,2]-score2[,2]
		#set the limits on the values of the insulation score
		delta.score[delta.score[,2] < zlim[1],2] <- zlim[1]
		delta.score[delta.score[,2] > zlim[2],2] <- zlim[2]
		x1 <- c(delta.score[1,1], (head(delta.score[,1],-1)+tail(delta.score[,1],-1))/2)
		x2 <- c((head(delta.score[,1],-1)+tail(delta.score[,1],-1))/2, tail(delta.score[,1],1))
		y1 <- window.size-step/2
		y2 <- window.size+step/2
		color.index <- (color.res-1)*(delta.score[,2] - zlim[1])/(zlim[2]-zlim[1]) + 1
		rect(x1, y1, x2, y2, col=color.vec[color.index], border=NA)
	}
	at <- seq(1e6*floor(start/1e6), 1e6*floor(end/1e6), by=1e6)
	axis(1, at=at, lab=at/1e6, lwd=2)
	axis(2, lwd=2)
	box(lwd=2)
}


#' insulation.domainogram
#'
#' create a domainogram for the insulation scores
#'
#' @param ins.score Scores.
#' @return happiness
insulation.polygon <- function( ins.score, rotate=F ){

	x <- c(ins.score[1,1],ins.score[,1],tail(ins.score[,1],1))
	y.up   <- c(0,ifelse(ins.score[,2] < 0, 0, ins.score[,2]),0)
	y.down <- c(0,ifelse(ins.score[,2] > 0, 0, ins.score[,2]),0)

	if(rotate){
		polygon(y.up, x, col=rgb(1,0,0,0.3), border=NA)
		polygon(y.down, x, col=rgb(0,0,1,0.3), border=NA)

	}else{
		polygon(x, y.up, col=rgb(1,0,0,0.3), border=NA)
		polygon(x, y.down, col=rgb(0,0,1,0.3), border=NA)
	}
}

#' insulation.score
#'
#' Get insulation scores
#'
#' @param exp1 A Hi-C experiment object: produced by construct.experiment().
#' @param chrom Chromosome
#' @param start Start position of the region of interest
#' @param end End position of the region of interest
#' @param cut.off The cut.off for Hi-C scores
#' @param window.size The sliding square size
#' @param step Thesliding square step-size
#' @param local Per chromosome?
#' @param diag.add Add values to diaginal
#' @return A plot
#' @export
insulation.score <- function( hic, window.size, chrom, start, end, diag.add = 0, local = T ){
	if((end - start)/hic$RES > 1000 + 2*window.size){
		stop("Please use a region that is smaller than 1000 times the resolution (+2 times the window size)")
	}
	if(window.size %% 2 == 0){
		stop("Please use an odd window size")
	}
	mat <- select.subset( hic$ICE, chrom, start, end, hic$ABS)
	ins.score <- matrix.insulation( mat, window.size )
	if(local){
		ins.score[,2] <- log2(ins.score[,2]/mean(ins.score[,2], na.rm=T))
	}else{
		chrom.data <- chromosome.wide.insulation( hic, window.size, chrom )
		ins.score[,2] <- log2(ins.score[,2]/median(chrom.data[,2], na.rm=T))
	}
	ins.score
}

#' chromosome.wide.insulation
#'
#' get insulation for matrix
#'
#' @param mat exp
#' @param window.size The sliding square size
#' @return happiness
matrix.insulation <- function( mat, window.size ){
	m <- matrix(1, ncol=window.size, nrow=window.size)
	id <- which(m > 0, arr.ind=T)
	x <- id[,1]; y <- id[,2]
	num.windows <- nrow(mat$z) - 2*window.size
	x.add <- rep((window.size):(num.windows+window.size), each = window.size**2)
	y.add <- x.add - window.size
	x <- rep(x, num.windows+1)
	y <- rep(y, num.windows+1)

	pos.id <- window.size:(num.windows+window.size)
	#use x.add as an index for the window
	score.vec <- mat$z[cbind(x+x.add, y+y.add)]
	insulation <- tapply(score.vec, x.add, mean)
	#insulation is the interaction score divided by the average interaction score of the entire region
	data.frame(mat$x[pos.id],insulation)
}


#' chromosome.wide.insulation
#'
#' loop over the entire chromosome and calcute the insulation score for every position, use this to correct if requeste
#'
#' @param hic exp
#' @param window.size The sliding square size
#' @return happiness
chromosome.wide.insulation <- function( hic, window.size, chrom ){
	max.pos <- max(hic$ABS[hic$ABS[,1]==chrom,3])
	window = hic$RES*1000
	chrom.ins.vec <- NULL
	for( start in seq(0, max.pos, by = window ) ){
		#nothing to be gained here new matrix will extend beyond the bound box
		#the 2-fold is a quick fix, we should probably think of a better method
		if(start + 2*window.size*hic$RES > max.pos){
			break
		}
		end <- start+window+window.size
		mat <- select.subset( hic$ICE, chrom, start, end, hic$ABS)
		ins.vec <- matrix.insulation( mat, window.size )
		chrom.ins.vec <- rbind(chrom.ins.vec, ins.vec )
	}
	chrom.ins.vec
}

#' genome.wide.insulation
#'
#' Get genome-wide insulation
#'
#' @param hic exp
#' @param window.size The sliding square size
#' @return happiness
genome.wide.insulation <- function( hic, window.size, normalize.genome = F ){
	chrom.vec <- c(hic$ABS[1,1],hic$ABS[which(head(hic$ABS[,1],-1) != tail(hic$ABS[,1],-1))+1,1])
	chrom.save <- c(); pos.save <- c(); ins.vec <- c()
	for( chrom in chrom.vec ){
		cat("Currently analyzing ", chrom, "\r")
		#do not analyse chromosomes that are twice the window size
		if(sum(hic$ABS[,1]==chrom) <= 2*window.size ){
			next
		}
		ins.data <- chromosome.wide.insulation( hic, window.size, chrom )
		#normalize per chromosome
		if( ! normalize.genome ){
			ins.data[,2] <- log2(ins.data[,2]/median(ins.data[,2]))
		}
		chrom.save <- c(chrom.save,rep(chrom, nrow(ins.data)))
		pos.save <- c(pos.save, ins.data[,1])
		ins.vec <- c(ins.vec, ins.data[,2])
	}
	if(normalize.genome){
		ins.vec <- log2(ins.vec/median(ins.vec))
	}
	#note substract/add half of the resolution because we create a bed-like structure
	data.frame(chrom=chrom.save, start=pos.save-hic$RES/2, end=pos.save+hic$RES/2, insulation=ins.vec)
}
