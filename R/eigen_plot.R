#create an observed over expected matrix
obs.exp.matrix <- function( mat, correct=F, outlier.correct = 0.995, lowess = F ){
	pos <- which(mat$z > -1, arr.ind=T)
	val <- mat$z[mat$z > -1]
	d <- abs(pos[,1]-pos[,2])
	if(correct){
		norm.factor <- tapply(val, d, function(x) mean( x[x < quantile(x, outlier.correct)] ) )
	}else{
		norm.factor <- tapply(val, d, mean)
	}
	if(lowess){
		x <- as.numeric(names(norm.factor)); y <- norm.factor
		norm.lowess <- lowess(x, y, f = 0.01, iter=20)
		norm.factor <- norm.lowess$y
		names(norm.factor) <- norm.lowess$x
	}

	val <- val/norm.factor[as.character(d)]
	obs.exp <- matrix(val, ncol = length(val)**0.5)
	list(x=mat$x, y=mat$y, z=obs.exp)
}


#input is an observed over expected matrix
get.eigen <- function( mat, with.eigen.value = T, outlier.correct = 0.995 ){
	oe <- obs.exp.matrix( mat )
	#set non-finite values to 1
	oe$z[!is.finite(oe$z)] <- 1
	#remove outliers
	th <- quantile(oe$z, outlier.correct)
	oe$z[oe$z > th] <- th

	ev <- eigen(oe$z - 1)
	ev
}

eigen.sort <- function( mat, outlier.correct = 0.995 ){
	oe <- obs.exp.matrix ( mat )
	#set non-finite values to 1
	oe$z[!is.finite(oe$z)] <- 1
	#remove outliers
	th <- quantile(oe$z, outlier.correct)
	oe$z[oe$z > th] <- th

	ev <- eigen(oe$z - 1)

	oe.sort <- oe
	mat.order <- order(ev$vector[,1])
	oe.sort$z <- oe.sort$z[mat.order,mat.order]
	oe.sort$x <- sort(ev$vector[,1])
	oe.sort$y <- sort(ev$vector[,1])

	oe.sort
}

eigen.struct <- function( mat, outlier.correct = 0.995 ){
	oe <- obs.exp.matrix ( mat )
	#set non-finite values to 1
	oe$z[!is.finite(oe$z)] <- 1
	#remove outliers
	th <- quantile(oe$z, outlier.correct)
	oe$z[oe$z > th] <- th

	ev <- eigen(oe$z - 1)

	oe$ev1 <- ev$vector[,1]*(ev$value[1]**0.5)
	oe$ev2 <- ev$vector[,2]*(ev$value[2]**0.5)

	oe
}

#select a matrix of interactions for between two chromosomes
selectData <- function (exp, chrom1, chrom2){
	bed <- exp$ABS
	data <- exp$ICE
	X <- bed[bed[,1]==chrom1,4]
	Y <- bed[bed[,1]==chrom2,4]
	#the order of the chromosomes matters for the analysis
	#make sure that X is smaller than Y, otherwise switch
	#them around
	if(X[1] > Y[1]){
		temp <- X
		X <- Y
		Y <- temp
		temp <- chrom1; chrom1 <- chrom2; chrom2 <- temp; #switch the chromosomes around as well
	}
	#create x and y vectors that contain the positions of the
	#entries in the matrix that we are creating
	x <- rep(X[1]:X[length(X)],        tail(Y, n=1) - Y[1] + 1)
	y <- rep(Y[1]:Y[length(Y)], each = tail(X, n=1) - X[1] + 1)
	data.sub <- data[base::list(x, y)]
	data.sub <- data.sub[!is.na(data.sub$V3)]
	#create an empty matrix, that has as many rows as the 'X' chromosome has
	#windows and as many columns as the 'Y' chromosome has windows
	mat <- matrix(0, ncol=tail(Y, n=1) - Y[1] + 1, nrow=tail(X, n=1) - X[1] + 1)
	mat[cbind(data.sub$V1-min(X)+1, data.sub$V2-min(Y)+1)] <- data.sub$V3
	x.pos <- bed[bed[,1]==chrom1,2]
	y.pos <- bed[bed[,1]==chrom2,2]
	#create a list that is compatible with the image function
	mat <- list(x=x.pos, y=y.pos, z=mat)
	mat
}

select.cis.arm <- function( mat, cp, arm = NULL ){
	if(is.null(arm)){ stop( "No arm selected" ) }

	if(arm == "p"){
		sel.i <- which(mat$x < cp[1,2]); sel.j <- which(mat$y < cp[1,2]);
	}else if(arm == "q"){
		sel.i <- which(mat$x > cp[1,3]); sel.j <- which(mat$y > cp[1,3]);
	}

	mat.new <- list(x=mat$x[sel.i],y=mat$y[sel.j], z=mat$z[sel.i,sel.j])
	mat.new$z[lower.tri(mat.new$z)] <- t(mat.new$z)[lower.tri(mat.new$z)]
	mat.new
}

#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch <- function( x ){
	temp <- cumsum(c(1,diff(x) - 1))
	temp2 <- rle(temp)
	x[which(temp == with(temp2, values[which.max(lengths)]))]
}

#remove the outlier rows, i.e. the rows that a lower score than the quantile
#divided by some value
remove.outliers <- function( arm, q=0.05, th.mult = 1 ){
	row.scores <- apply(arm$z, 1, mean)
	cut.off <- quantile(row.scores, 0.05)
	th <- cut.off/th.mult
	print(th)

	sel <- which(row.scores < th)

	if(length(sel) > 0){
		arm$x <- arm$x[-sel]
		arm$y <- arm$y[-sel]
		arm$z <- arm$z[-sel,-sel]
	}
	arm
}

#switch the eigen vector based on a chip track of for instance
#active histone marks
switch.EV <- function( ev.data, chip, chrom ){
        start <- min(ev.data$x); end <- max(ev.data$x)
        sub.chip <- chip[chip[,1]==chrom & chip[,2] > start & chip[,2] < end,]
        window.cnt <- findInterval(sub.chip[,2], ev.data$x)
        window.cnt <- table(factor(window.cnt,1:length(ev.data$x)))
        #return true or false depending on whether the median value of
        #ChIP in the up or down compartment is highests (i.e. true if
        #down scores have the highest number of ChIP peaks)
        #up comp.                             down comp.
        median(window.cnt[ev.data$ev1 > 0]) < median(window.cnt[ev.data$ev1 < 0])
}






#' cis.compartment.plot
#'
#' Draw intrachromosomal interaction heatmap for a chromosome (arm) with corresponding compartment scores
#'
#' @param exp1 A GENOVA experiment-object
#' @param exp2 A GENOVA experiment-object. If given, all plots show this exp in the lower-left corner.
#' @param chrom The chromosome that should be drawn
#' @param arm Which chromosome arm: 'p' or 'q' (for acrocentric this should be set to 'q')
#' @param cut.off maximum value for the heatmap
#' @param obs.exp whether an observed over expected matrix should be drawn (default is false)
#' @param invert whether the compartment score should be inverted
#' @param color.scheme color scheme that should be used, defaults to fall, other values result in white-red gradient
#' @param cs.lim y-axis limit for the compart score, if unset (NULL) will default to the maximum absolute value
#' @param chip A data.frame, containg ChIP-seq peaks of active histone marks to correctly orient A/B compartments
#' @param smoothNA Set to TRUE to perform a Nadaraya/Watson normalization. This will try to eliminate white stripes: this is only cosmetic and has no effect on the compartment-scores.
#' Requires fields and is only needed if large white stripes are bothering you.
#' @note
#' # Plot a cis-compartment plot of the q-arm of chromosome 14.
#' cis.compartment.plot(exp = Hap1_WT_40kb, chrom = 'chr14', arm = 'q', cs.lim = 1.75, cut.off = 15, chip = H3K27ac_peaks)
#'
#' # Plot a observed/expected cis-compartment plot of the q-arm of chromosome 14.
#' cis.compartment.plot(exp = Hap1_WT_40kb, chrom = 'chr14', arm = 'q', cs.lim = 1.75, cut.off = 15, chip = H3K27ac_peaks, obs.exp = T)
#' @export
cis.compartment.plot <- function( exp1, exp2 = NULL, chrom, arm="p", cut.off=NULL, obs.exp = F, invert=F, color.scheme="fall", cs.lim=NULL, chip = NULL, smoothNA = F){
	#exp1 = exp

  mat1 <- selectData( exp1, chrom, chrom)
  if(!is.null(exp2)){
    mat2 <- selectData( exp2, chrom, chrom)
  }

	cent <- which(apply(mat1$z,1,sum)==0)
	cent <- largest.stretch(cent)


	centromere.pos <- data.frame(chrom=c(chrom,chrom), start=c(min(cent)*exp1$RES, min(cent)*exp1$RES), end=c(max(cent)*exp1$RES, min(cent)*exp1$RES))
	if( centromere.pos[1,2] < 1e6 && arm == 'p'){
		print("No p arm. Enter q as arm")
		return(NULL)
	}





	#plot the p and q arms

	arm1 <- select.cis.arm(mat1, centromere.pos, arm)
	oe1 <- eigen.struct( arm1 )

	if(!is.null(exp2)){
	  arm2 <- select.cis.arm(mat2, centromere.pos, arm)
	  oe2 <- eigen.struct( arm2 )
	}


	#orient compartment score
	if(!is.null(chip)){
		if(switch.EV( oe1, chip, chrom) ){
			oe1$ev1 = -oe1$ev1
		}
		#note that if you add a chip track the invert option
		#is overridden
		invert=F
	}

	if(!is.null(exp2)){
	  if(!is.null(chip)){
	    if(switch.EV( oe2, chip, chrom) ){
	      oe2$ev1 = -oe2$ev1
	    }
	    #note that if you add a chip track the invert option
	    #is overridden
	    invert=F
	  }
	}

	if(!is.null(exp2)){
	  arm1$z[upper.tri(arm1$z)] <- arm2$z[upper.tri(arm2$z)]
	  oe1$z[upper.tri(oe1$z)] <- oe2$z[upper.tri(oe2$z)]
	}


	#create a plotting layout
	w = 6
	lay <- matrix(4, nrow=w, ncol=w)
	lay[2:w,2:w] <- 1; lay[1,] <- 2; lay[,1] <- 3; lay[1,1] <- 4
	layout(lay)
	par(mar=rep(1,4), xaxs="i", yaxs="i")

	if(smoothNA){
	  require(fields)
	  oe1.5 = oe1
	  oe1.5$z[oe1.5$z == 0] <- NA
	  oe1.5$z = fields::image.smooth(oe1.5$z, theta = 0.25)$z
	  oe1 = oe1.5

	  arm1.5 = arm1
	  arm1.5$z[arm1.5$z == 0] <- NA
	  arm1.5$z = fields::image.smooth(arm1.5$z, theta = 0.25)$z
	  arm1 = arm1.5
	}

	if(obs.exp){
		#color scale used: blue, white, red
		bwr <- colorRampPalette(c("blue","white","red"))

		log.mat <- log2(oe1$z)

		if(is.null(cut.off)){
		  cut.off = max(quantile(log.mat, .99))
		  warning("No cut.off was given: using 99 percentile: ", round(cut.off), ".")
		}

		log.mat[log.mat >  cut.off] <- cut.off
		log.mat[log.mat < -cut.off] <- -cut.off
		oe1$z <- log.mat
		image(oe1, col=bwr(300), zlim=c(-cut.off,cut.off), axes=F, ylim=rev(range(oe1$y)))


	}else{
		if(color.scheme=='fall'){
			col.fun <- colorRampPalette(c("white", "orange", "darkred", "black"))
		}else{
			col.fun <- colorRampPalette(c("white", "red"))
		}

	  if(is.null(cut.off)){
	    cut.off = max(quantile(arm1$z, .99))
	    warning("No cut.off was given: using 99% percentile: ", round(cut.off), ".")
	  }

		arm1$z[arm1$z > cut.off] <- cut.off
		image(arm1, col=col.fun(300), zlim=c(0,cut.off), axes=F, ylim=rev(range(oe1$y)))
	}
	box(lwd=2)

	label <- seq( 10*floor(min(arm1$x)/10e6), 10*floor(max(arm1$x)/10e6), by=10 )
	at <- label*1e6
	axis(3, at, label)
	axis(2, at, label)

	#plot the first compartment
	compartment.score.plot( oe1, invert=invert, cs.lim=cs.lim)

	#plot the second compartment
	if(!is.null(exp2)){
	  compartment.score.plot( oe2, invert=invert, cs.lim=cs.lim, rotate=T)
	} else {
	  compartment.score.plot( oe1, invert=invert, cs.lim=cs.lim, rotate=T)
	}

	oe1$raw <- arm1$z
	invisible(oe1)
}


#select the arm combination of the trans interactions
select.trans.arm <- function( mat, arm1 = NULL, arm2 = NULL, cp ){
	if(is.null(arm1) || is.null(arm2) ){ stop("No arm defined") }

	#select columns values for the first chromosome
	if(arm1 == "p"){
		sel.i <- which(mat$x < cp[1,2]);
	}else if(arm1 == "q"){
		sel.i <- which(mat$x > cp[1,3]);
	}

	#select row values for the second chromosome
	if(arm2 == "p"){
		sel.j <- which(mat$y < cp[2,2]);
	}else if(arm2 == "q"){
		sel.j <- which(mat$y > cp[2,3]);
	}
	mat.new <- list(x=mat$x[sel.i],y=mat$y[sel.j], z=mat$z[sel.i,sel.j])
	mat.new
}


cis.ev <- function( data, chrom, arm="p"){
	mat <- selectData( data, chrom, chrom)

	cent <- which(apply(mat$z,1,sum)==0)
	cent <- largest.stretch(cent)


	centromere.pos <- data.frame(chrom=c(chrom,chrom), start=c(min(cent)*data$RES, min(cent)*data$RES), end=c(max(cent)*data$RES, min(cent)*data$RES))
	if( centromere.pos[1,2] < 1e6 && arm == 'p'){
		print("No p arm. Enter q as arm")
		return(NULL)
	}

	#plot the p and q arms

	arm <- select.cis.arm(mat, centromere.pos, arm)
	oe <- eigen.struct( arm )
	oe$cent <- centromere.pos
	oe
}

#function to determine whether chromosomes need to be switched
switch.chromosomes <- function( data, chrom1, chrom2 ){
	bed <- data$ABS
	X <- bed[bed[,1]==chrom1,4]
	Y <- bed[bed[,1]==chrom2,4]
	return(X[1] > Y[1])
}


#' trans.compartment.plot
#'
#' Draw interchromosomal interaction heatmap for a chromosome (arm) with corresponding (cis) compartment scores
#'
#' @param exp A GENOVA experiment-object
#' @param chrom1,chrom2 Which chromosome-combination should be drawn
#' @param arm1,arm2 which chromosome arm: 'p' or 'q' (for acrocentric this should be set to 'q')
#' @param cut.off maximum value for the heatmap
#' @param invert whether the compartment score should be inverted, this should be a vector with a logical value for each chromosome
#' @param color.scheme color scheme that should be used, defaults to fall, other values result in white-red gradient
#' @param cs.lim y-axis limit for the compartment-score, if unset (NULL) will default to the maximum absolute value
#' @param chip A data.frame, containg ChIP-seq peaks of active histone marks to correctly orient A/B compartments
#' @note
#' # Plot a trans-compartment plot of the q-arms of chromosome 9 and 22
#' trans.compartment.plot(exp = Hap1_WT_40kb, chrom1 = 'chr9', arm1 = 'q', chrom2 = 'chr22', arm2 = 'q', cut.off = 10,chip = H3K27ac_peaks)
#' @export
trans.compartment.plot <- function( exp, chrom1, arm1="p", chrom2, arm2 = "p", cut.off=NULL, invert=c(F,F), color.scheme="fall", cs.lim=NULL, chip=NULL){
  data = exp
	#error handling
	if( length( invert ) != 2 ){
		stop("invert option should be a vector of length 2: for the first chromosome and for the second chromosome")
	}
	#get the compartment scores for the chromosomes
	oe1 <- cis.ev( data, chrom1, arm1)
	oe2 <- cis.ev( data, chrom2, arm2)

	if( is.null(oe1) || is.null(oe2) ){
		return(NULL)
	}

	#orient compartment score
	if(!is.null(chip)){
		if(switch.EV( oe1, chip, chrom1) ){
			oe1$ev1 = -oe1$ev1
		}
		if(switch.EV( oe2, chip, chrom2) ){
			oe2$ev1 = -oe2$ev1
		}
		#note that if you add a chip track the invert option
		#is overridden
		invert=c(F,F)
	}

	centromere.pos <- rbind(oe1$cent[1,], oe2$cent[1,])

	mat <- selectData( data, chrom1=chrom1, chrom2=chrom2 )
	#if chromosomes were added in the "wrong" orientation, change the orientation
	#of the matrix and and the x and y positions
	#chroms <- switch.chromosomes(data, chrom1, chrom2 ) #determine the order of the chromsomes
	if(switch.chromosomes(data, chrom1, chrom2 )){
		temp <- mat$x; mat$x <- mat$y; mat$y <- temp #alternative then the switching of chromosomes is not needed
		mat$z <- t(mat$z)
	}
	arm <- select.trans.arm(mat, arm1=arm1, arm2=arm2, cp = centromere.pos )


	#create a plotting layout
	w = 6
	lay <- matrix(4, nrow=w, ncol=w)
	lay[2:w,2:w] <- 1; lay[1,] <- 2; lay[,1] <- 3; lay[1,1] <- 4
	layout(lay)
	par(mar=rep(1,4), xaxs="i", yaxs="i")

	if(color.scheme=='fall'){
		col.fun <- colorRampPalette(c("white", "orange", "darkred", "black"))
	}else{
		col.fun <- colorRampPalette(c("white", "red"))
	}

	if(is.null(cut.off)){
	  cut.off = max(quantile(arm$z, .99))
	  warning("No cut.off was given: using 99% percentile: ", round(cut.off), ".")
	}

	arm$z[arm$z > cut.off] <- cut.off
	image(arm, col=col.fun(300), zlim=c(0,cut.off), axes=F, ylim=rev(range(arm$y)))
	box(lwd=2)

	label <- seq( 10*floor(min(arm$x)/10e6), 10*floor(max(arm$x)/10e6), by=10 )
	at <- label*1e6
	axis(3, at, label)
	label <- seq( 10*floor(min(arm$y)/10e6), 10*floor(max(arm$y)/10e6), by=10 )
	at <- label*1e6
	axis(2, at, label)

	#plot the first compartment
	compartment.score.plot( oe1, invert=invert[1], cs.lim=cs.lim)

	#plot the second compartment
	compartment.score.plot( oe2, invert=invert[2], cs.lim=cs.lim, rotate=T)

}

compartment.score.plot <- function( oe, invert=F, cs.lim=NULL, rotate=F){
	x.pos <- oe$x; y.pos <- oe$ev1
	#invert the compartment score if necessary
	if(invert){
		y.pos <- -y.pos
	}

	if(is.null(cs.lim)){
		cs.lim <- max(abs(y.pos))
	}
	if(!rotate){
		plot(x.pos, y.pos, type='n', xaxt='n', axes=F, ylim=c(-cs.lim,cs.lim))
		ab.polygon(x.pos, y.pos)
		axis(2)
	}else{
		plot(y.pos, x.pos, type='n', axes=F, xlim=c(cs.lim,-cs.lim), ylim=rev(range(x.pos)))
		ab.polygon(x.pos, y.pos, rotate=T)
		axis(3)
	}
}

#draw a polygon for the compartment scores
#up is a red polygon, down is a blue polygon
ab.polygon <- function( x.pos, y.pos, rotate=F){
	x <- c(x.pos[1], x.pos, tail(x.pos,1))
	y.up   <- c(0, ifelse(y.pos < 0, 0, y.pos), 0 )
	y.down <- c(0, ifelse(y.pos > 0, 0, y.pos), 0 )

	if(rotate){
		polygon(y.up, x, col=rgb(1,0,0,0.8), border=NA)
		polygon(y.down, x, col=rgb(0,0,1,0.8), border=NA)
	}else{
		polygon(x, y.up, col=rgb(1,0,0,0.8), border=NA)
		polygon(x, y.down, col=rgb(0,0,1,0.8), border=NA)
	}
}
