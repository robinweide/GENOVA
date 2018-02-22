
resize.mat <- function(mat, ndim=dim(mat)){
  # input object
  odim <- dim(mat)
  obj <- list(x= 1:odim[1], y=1:odim[2], z= mat)
  # output object
  ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])
  ndim <- dim(ans)
  # rescaling
  ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
  loc <- ncord
  loc[,1] = rescale(ncord[,1], c(1,odim[1]))
  loc[,2] = rescale(ncord[,2], c(1,odim[2]))
  # interpolation
  ans[ncord] <- fields.interp.surface(obj, loc)
	ans[is.na(ans)] <- 0
  ans
}

rescale <- function(x, newrange=range(x)){
  xrange <- range(x)
  mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  newrange[1]+(x-xrange[1])*mfac
}

fields.interp.surface <- function(obj, loc) {

    # obj is a surface or image  object like the list for contour, persp or image.
    # loc a matrix of 2 d locations -- new points to evaluate the surface.
    x <- obj$x
    y <- obj$y
    z <- obj$z
    nx <- length(x)
    ny <- length(y)
    # this clever idea for finding the intermediate coordinates at the new points
    # is from J-O Irisson
    lx <- approx(x, 1:nx, loc[, 1])$y
    ly <- approx(y, 1:ny, loc[, 2])$y
    lx1 <- floor(lx)
    ly1 <- floor(ly)
    # x and y distances between each new point and the closest grid point in the lower left hand corner.
    ex <- lx - lx1
    ey <- ly - ly1
    # fix up weights to handle the case when loc are equal to
    # last grid point.  These have been set to NA above.
    ex[lx1 == nx] <- 1
    ey[ly1 == ny] <- 1
    lx1[lx1 == nx] <- nx - 1
    ly1[ly1 == ny] <- ny - 1
    # bilinear interpolation finds simple weights based on the
    # the four corners of the grid box containing the new
    # points.
    return(z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + z[cbind(lx1 +
        1, ly1)] * ex * (1 - ey) + z[cbind(lx1, ly1 + 1)] * (1 -
        ex) * ey + z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}

#select a matrix of interactions for between two chromosomes
selectTransData <- function (exp, chrom1, chrom2){
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


#select an arm combination of the chromosome, with two chromosomes
#that both have two arms, there are 4 possible arm combinations
#the resulting matrix is reoriented so that the centromere-centromere
#interaction is at (1,1)(bottomleft) and the telomere-telomere interaction is at
#(nrow,ncol)(topright)
select.arm <- function( mat, cp, quadrant = NULL ){
	if(is.null(quadrant)){ stop( "No quadrant selected" ) }

	if(quadrant==1){
		sel.i <- which(mat$x < cp[1,2]); sel.j <- which(mat$y < cp[2,2]);
		sel.i <- rev(sel.i); sel.j <- rev(sel.j); #reorientation
	}else if(quadrant==2){
		sel.i <- which(mat$x > cp[1,3]); sel.j <- which(mat$y < cp[2,2]);
		sel.j <- rev(sel.j) #reorientation
	}else if(quadrant==3){
		sel.i <- which(mat$x < cp[1,2]); sel.j <- which(mat$y > cp[2,3]);
		sel.i <- rev(sel.i); #reorientation
	}else if(quadrant==4){
		sel.i <- which(mat$x > cp[1,3]); sel.j <- which(mat$y > cp[2,3]);
	}

	mat$z[sel.i,sel.j]
}

#find the largest stretch of 0s, which is
#most likely the centromere
largest.stretch <- function( x ){
	temp <- cumsum(c(1,diff(x) - 1))
	temp2 <- rle(temp)
	x[which(temp == with(temp2, values[which.max(lengths)]))]
}

#draw a cartoon chromosome
draw.chromosome <- function(){

	x <- c(0,1)
	y <- c(1.05,1.05)

	#telomere
	phi <- seq(0,pi,len=100)
	x <- c(x, 0.02*sin(phi)+1)
	y <- c(y, 0.015*cos(phi)+1.035)

	x <- c(x, 1, 0)
	y <- c(y, 1.02, 1.02)

	#centromere
	phi <- seq(0,2*pi, len=300)
	phi <- rev(phi[c(37:1,300:114)])
	x <- c(x, 0.02*sin(phi)-0.7*0.02)
	y <- c(y, 0.02062593*cos(phi)+1.035)

	range.val <- abs(diff(range(x)))
	x<-(x - min(x))/range.val
	y <- 1.035+(y-1.035)/range.val
	#horizontal chromosome
	polygon(x,y, col="grey", border="black", lwd=2)
	#vertical chromosome
	polygon(y-1.07,1-x, col="grey", border="black", lwd=2)
}

#' draw.centromere.telomere
#'
#' Plot the log2-ratio of the relative interaction frequencies.
#' The ratio is calculated over the median of the matrix.
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param m matrix result from centromere-telomere
#' @param cut.off log2-ratio cut.off
#' @examples
#' # Get a scaled matric of the interchomosomal interactions between 15 and 19
#' out1519 = centromere.telomere.analysis(WT_40kb,
#'                                        chrom.vec = c('chr15', 'chr19'))
#'
#' # Plot the results
#' draw.centromere.telomere(out1519)
#' @export
draw.centromere.telomere <- function(m, cut.off = 2){
	bwr <- colorRampPalette(c("blue", "white", "red"))

	#normalize the contact frequencies by dividing by the median
	lm <- log2(m/median(m))

	#maximize the maximum and minimum threshold
	lm[lm > cut.off] <- cut.off
	lm[lm < -cut.off] <- -cut.off
	image(lm[,nrow(lm):1],
	      axes=F,
	      xlim=c(-0.06, 1.01),
	      ylim=c(-0.01, 1.06),
	      col=bwr(1000),
	      zlim=c(-cut.off, cut.off))

	#add cartoon chromosomes
	draw.chromosome()
}


#' centromere.telomere.analysis
#'
#' Calculate a scaled matrix of the interchromosomal interactions aligned from
#' centromere to telomere
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp GENOVA data structure for your experiment of interest
#' @param chrom.vec a vector containing the chromosomes that should be
#' considered in the pair-wise interchromosomal interactions
#' @param nrow dimensions of the matrix (number of columns will be the same)
#' @param leave.out two-column matrix or data.frame containing the specific
#' chromosome combinations that need to left out
#' @param q.top top quantile of scores that should be left out of the analysis
#' (because they are outliers)
#' @param verbose Produces a progress-indication.
#' @examples
#' # Get a scaled matric of the interchomosomal interactions between 15 and 19
#' out1519 = centromere.telomere.analysis(WT_40kb,
#'                                        chrom.vec = c('chr15', 'chr19'))
#'
#' # Plot the results
#' draw.centromere.telomere(out1519)
#' @export
centromere.telomere.analysis <- function(exp, chrom.vec, nrow = 100,
                                         leave.out = NULL, q.top = 1e-5,
                                         verbose = F){
	i.vec <- 1:(length(chrom.vec)-1)
	# empty matrix for holding the contact frequencies
	m.total <- matrix(0, nrow = nrow, ncol = nrow)
	# loop over the chromosome combinations
	for(i in i.vec){
		for(j in (i+1):(max(i.vec)+1)){
			chrom1 <- chrom.vec[i]
			chrom2 <- chrom.vec[j]

			# if the leave.out dataset is defined check if the chromosome combination
			# is found
			if(!is.null(leave.out)){
				sum.val <- sum(leave.out[,1] == chrom1 & leave.out[,2] == chrom2) +
				           sum(leave.out[,1] == chrom2 & leave.out[,2] == chrom1)
				if(sum.val > 0){
					next
				}
			}
			if(verbose){
			  message(chrom1, "\t", chrom2, "\r")
			}
			# select the interaction matrix between the trans chromosomes
			trans.mat <- selectTransData(exp, chrom1, chrom2)

			# set the q.top highest values to this value
			th <- quantile(trans.mat$z, 1-q.top)
			trans.mat$z[trans.mat$z > th] <- th

			# get the centromere positions emprically
			cent1 <- which(apply(trans.mat$z,1,sum) == 0)
			cent2 <- which(apply(trans.mat$z,2,sum) == 0)

			cent1 <- largest.stretch(cent1)
			cent2 <- largest.stretch(cent2)

			res <- exp$RES
			centromere.pos <- data.frame(chrom = c(chrom1,chrom2),
			                             start = c(min(cent1)*res, min(cent2)*res),
			                             end   = c(max(cent1)*res, max(cent2)*res))

			# Q1,2,3,4 are the different quadrants of the chromosome-chromosome
			# combinations
			# +-----------+
			# |  1  |  2  |
			# +-----C-----+
			# |  3  |  4  |
			# +-----------+
			for(quadrant in 1:4){
				m.sub <- select.arm(trans.mat, centromere.pos, quadrant)
				# acrocentric chromosomes happen
				if(is.null(dim(m.sub))){
				  next
				}
				m.sub <- resize.mat(m.sub, c(nrow, nrow))
				# calculate the average of the transposed matrix and the regular matrix
				# to get rid of chromosome specific enrichments
				m.total <- m.total + (m.sub + t(m.sub))/2
			}

		}
	}
	m.total
}

inter.arm.data <- function(exp, chrom.vec, nrow = 100,
                           leave.out = NULL, q.top = 1e-5){

	m.total <- matrix(0, nrow = nrow, ncol = nrow)
	for( chrom in chrom.vec ){
		cis.mat <- selectTransData(exp, chrom, chrom )

		#set the q.top highest values to this value
		th <- quantile(cis.mat$z, 1-q.top)
		cis.mat$z[cis.mat$z > th] <- th

		#get the centromere positions emprically
		cent <- which(apply(cis.mat$z,1,sum)==0)
		cent <- largest.stretch(cent)

		if(min(cent) > 1 ){
			rows <- 1:min(cent)
			cols <- max(cent):nrow(cis.mat$z)
			inter.arm <- cis.mat$z[rows,cols]
			inter.arm <- resize.mat(inter.arm, c(nrow, nrow))
			#m.total <- m.total + (inter.arm[nrow:1,] + t(inter.arm[nrow:1,]))
			m.total <- m.total + inter.arm[nrow:1,]
		}
	}
	m.total
}

draw.interarm <- function( exp, chrom.vec, q.top = 1e-5, q.z=0.99 ){
	for( chrom in chrom.vec ){
		cis.mat <- selectTransData( exp, chrom, chrom )

		#set the q.top highest values to this value
		th <- quantile(cis.mat$z, 1-q.top)
		cis.mat$z[cis.mat$z > th] <- th

		#get the centromere positions emprically
		cent <- which(apply(cis.mat$z,1,sum)==0)
		cent <- largest.stretch(cent)

		if(min(cent) > 1 ){
			rows <- 1:min(cent)
			cols <- max(cent):nrow(cis.mat$z)
			inter.arm <- cis.mat$z[rows,cols]

			zmax <- quantile(inter.arm, q.z)
			inter.arm[inter.arm > zmax] <- zmax

			fall <- colorRampPalette(c("white", "orange", "darkred", "black"))
			image(list(x=cis.mat$x[rows], y=cis.mat$y[cols], z=inter.arm), zlim=c(0,zmax), main=chrom, col=fall(1000))
		}
	}
}


