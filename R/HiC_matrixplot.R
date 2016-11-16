draw_exon <- function( genes, chrom, y.pos, width, rotate=F ){
  genes.chr <- genes[genes$chrom==chrom,]

  es <- as.numeric(unlist(strsplit(as.character(genes.chr$exonStart), ",")))
  ee <- as.numeric(unlist(strsplit(as.character(genes.chr$exonEnd), ",")))

  strands <- rep(genes.chr$strand, unlist(lapply(strsplit(as.character(genes.chr$exonStart), ","), length)))

  add <- ifelse(strands=="+", 0.025*width+width, -0.025*width-width)
  add2 <- ifelse(genes.chr$strand=="+", 0.025*width+width, -0.025*width-width)
	if(rotate){
  	rect( (y.pos-width/2)+add, es, (y.pos+width/2)+add, es, col='black')
  	segments(y.pos+add2, genes.chr$txStart, y.pos+add2, genes.chr$txEnd)
	}else{
  	rect(es, (y.pos-width/2)+add, ee, (y.pos+width/2)+add, col='black')
  	segments(genes.chr$txStart, y.pos+add2, genes.chr$txEnd, y.pos+add2)
	}
}

plot.triangle <- function( bed, chrom, y1, y2, start, end, rotate=F ){
	y.scale <- abs(y1-y2)/0.15
  x.wid <- (end-start)*0.012*y.scale

	#first plot positive
  sel.bed <- bed[bed[,1]==chrom & bed[,2] < end & bed[,3] > start & bed[,6]=='+',]
	if(nrow(sel.bed)>0){
		col = "red"
	  add <- ifelse(sel.bed[,6]=='+', x.wid, -x.wid)
		add.list <- unique(add)
		if(rotate){
			triangle <- lapply( sel.bed[,2], function(x) list(yy=c(x,x,x+add.list), xx=c(y1,y2,(y1+y2)/2) ) )
			lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
		}else{
			#with polygon
			triangle <- lapply( sel.bed[,2], function(x) list(xx=c(x,x,x+add.list), yy=c(y1,y2,(y1+y2)/2) ) )
			lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
		}
	}

	#then repeat for negative
	sel.bed <- bed[bed[,1]==chrom & bed[,2] < end & bed[,3] > start & bed[,6]=='-',]
	if(nrow(sel.bed)>0){
		col="blue"
	  add <- ifelse(sel.bed[,6]=='+', x.wid, -x.wid)
		add.list <- unique(add)
		if(rotate){
			triangle <- lapply( sel.bed[,2], function(x) list(yy=c(x,x,x+add.list), xx=c(y1,y2,(y1+y2)/2) ) )
			lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
		}else{
			#with polygon
			triangle <- lapply( sel.bed[,2], function(x) list(xx=c(x,x,x+add.list), yy=c(y1,y2,(y1+y2)/2) ) )
			lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
		}
	}

}

plot.rectangle <- function( bed, chrom, y1, y2, start, end, col, rotate=F ){
	sel.bed <- bed[bed[,1]==chrom & bed[,2] < end & bed[,3] > start,]
	if(nrow(sel.bed)>0){
		if(rotate){
			rect(y1, sel.bed[,2], y2, sel.bed[,3], col=col, border=col)
		}else{
			rect(sel.bed[,2], y1, sel.bed[,3], y2, col=col, border=col)
		}
	}
}

plot.genes <- function( genes, chrom, start, end, y.pos, horiz=T){
	#plot the horizontal top info
	if(horiz){
		draw_exon( genes, chrom=chrom, y.pos = y.pos, width =0.05)
	}else{
	#plot the vertical top info
		draw_exon( genes, chrom=chrom, y.pos = y.pos, width =0.05, rotate=T)
	}
}

features.bed <- function(mat1, chrom, genes=NULL, chip1=NULL, chip2=NULL, y.values, type1="triangle", type2="triangle", col1="red", col2="blue",   rotate=F ){

	#calculate the appropriate y position values
	if(max(mat1$x)-min(mat1$x) > 2e6){
		gene.pos = -0.5
		chip1.y1 <- -0.85; chip1.y2 <- -0.80;
		chip2.y1 <- -0.75; chip2.y2 <- -0.70;
	}else{
		gene.pos = -0.2
		chip1.y1 <- -0.85; chip1.y2 <- -0.70;
		chip2.y1 <- -0.65; chip2.y2 <- -0.50;
	}
	y.values <- list(gene.pos=gene.pos, chip1.y1=chip1.y1, chip1.y2=chip1.y2, chip2.y1=chip2.y1, chip2.y2=chip2.y2)

	if(rotate){
		plot(0, type='n', ylim=rev(range(mat1$x)), xlim=c(0,-1), axes=F, xlab="", ylab="" )
	}else{
		plot(0, type='n', xlim=range(mat1$x), ylim=c(-1,0), axes=F, xlab="", ylab="" )
	}

	#plot the genes
	if(!is.null(genes)){
		plot.genes(genes, chrom, start, end, y.pos=y.values$gene.pos, rotate=rotate)
	}

	#plot the first chip dataset
	if(!is.null(chip1)){
		if(type1=="triangle"){
			plot.triangle( chip1, chrom=chrom, y1=y.values$chip1.y1, y2=y.values$chip1.y2, start=min(mat1$x), end=max(mat1$x), rotate=rotate)
		}else if(type1=="rectangle"){
			plot.rectangle(chip1, chrom=chrom, y1=y.values$chip1.y1, y2=y.values$chip1.y2, start=min(mat1$x), end=max(mat1$x), col=col1, rotate=rotate)
		}
	}

	#plot the second chip dataset
	if(!is.null(chip2)){
		if(type2=="triangle"){
			plot.triangle( chip2, chrom=chrom, y1=y.values$chip2.y1, y2=y.values$chip2.y2, start=min(mat1$x), end=max(mat1$x), rotate=rotate)
		}else if(type2=="rectangle"){
			plot.rectangle(chip2, chrom=chrom, y1=y.values$chip2.y1, y2=y.values$chip2.y2, start=min(mat1$x), end=max(mat1$x), col=col2, rotate=rotate)
		}
	}
	if(is.null(NULL)){
		#prevents printing random nulls to to stdout
	}
}

features <- function( mat1, chrom, genes=NULL, chip1=NULL, chip2=NULL, bed.col=c("red","blue"), bw.col=c("navy","darkred"), rotate=F, type=c("triangle","triangle") ){

	if(!is.null(chip1) && !is.null(chip2) && typeof(chip1) != typeof(chip2) ){
		stop("Feature tracks should both be bed structure or both be file path to bigWig file. Not mixed")
	}

	#if chip1 and/or chip2 are file paths do the bigwig analysis
	if(typeof(chip1)=="character" || typeof(chip2)=="character"){
		features.bw(mat1, chrom, chip1=chip1, chip2=chip2, col1=bw.col[1], col2=bw.col[2],  rotate=rotate)
	}


	if(typeof(chip1)=="list" || typeof(chip2)=="list"){
		features.bed(mat1, chrom, genes=genes, chip1=chip1, chip2=chip2, col1=bed.col[1], col2=bed.col[2], type1=type[1], type2=type[2], rotate=rotate)
	}
}

features.bw <- function( mat1, chrom, chip1=NULL, chip2=NULL, col1="navy", col2="darkred", rotate=F){
	start <- min(mat1$x); end <- max(mat1$x)
	if(rotate){
		plot(0, type='n', ylim=rev(range(mat1$x)), xlim=c(0,-1), axes=F, xlab="", ylab="" )
	}else{
		plot(0, type='n', xlim=range(mat1$x), ylim=c(-1,0), axes=F, xlab="", ylab="" )
	}

	if(!is.null(chip1)){
		plot.bw(chip1, chrom, start, end, -0.45, 0, col=col1, rotate=rotate)
	}
	if(!is.null(chip2)){
		plot.bw(chip2, chrom, start, end, -0.92, -0.47, col=col2, rotate=rotate)
	}
}

plot.bw <- function( file, chrom, start, end, y1, y2, col, rotate=F){
  if(require("xtable") == F){stop("Please install github.com/jayhesselberth/bigwrig\n")  }

	#leave here in case something changes
	#d <- read_bigwig( file, chrom=chrom, start=start, end=end) #official call but don't need it
	#d <- as.data.frame(d)

	d <- bigwrig::read_bigwig_impl( file, chrom=chrom, start=start, end=end)

	#to speed up plotting we are going to collapse the same consecutive
	#values
	i <- which(diff(d[,4])!=0) #select points where the consecutive values differ
	sel <- floor(c(0,head(i,-1))+i)/2 #select the middle of a set of consecutive values
	d <- d[sel,]

	max.val <- max(d[,4])
	y.range <- abs(y1-y2)
	y.val <- y.range*d[,4]/max.val
	if(rotate){
		segments(y1, d[,2], y1+y.val, d[,2], col=col)
	}else{
		segments(d[,2], y1, d[,2], y1+y.val, col=col)
	}
}

#' hic.matrixplot
#'
#' Plot a matrix (or two) for a region of interest with annotations
#'
#' @param exp1 The control Hi-C experiment object: produced by construct.experiment(). (bottom)
#' @param exp2 Optional: the treatment Hi-C experiment object: produced by construct.experiment().
#' @param chrom Chromosome
#' @param start Start position of the region of interest
#' @param end End position of the region of interest
#' @param cut.off The cut.off for the hic-matrix plot, in the diff option the negative of this is the lower bound
#' @param chip A list of feature tracks, can be bed structure (i.e. data frames) or a path to bigwig file (i.e. character variable), maximum length of 4, first two and last two have to be the same (i.e. 1:bed,2:bed,3:path,4:path, but 1:bed,2:NULL,3:path,4:NULL is also allowed)
#' @param bed.col Color of the bed track (max.len is 4)
#' @param bw.col Same as bed col, but for the bigwig track
#' @param type Should a rectangle or a triangle be drawn? Not that for a triangle a 6th strand column should be included
#' @param coplot When drawing together two experiments, dual is bottom triangle exp1, top triangle exp2; diff plots a substraction of exp2-exp1
#' @param genes Structure with gene information, will only be combined with bed structure
#' @return A matrix-plot
#' @export
hic.matrixplot <- function( exp1, exp2=NULL, chrom, start, end, cut.off=0, chip=list(NULL,NULL,NULL,NULL), bed.col=rep(c("red","blue"),2), bw.col=rep(c("navy","darkred"),2), type=rep("triangle",4), coplot="dual", genes=NULL){

	#some error handling
	if(!is.null(exp2)){
		#make sure the resolutions are the same
		if(exp1$RES != exp2$RES){
			stop("The Hi-C matrices should have the same resolution")
		}

		if(!all(exp1$ABS[,4] == exp2$ABS[,4])){
			stop("Not all ICE indexes are the same. Are you these experiments were mapped to the same genome (build)?")
		}
	}

	#if only one color is given use it for all feature tracks
	if(length(bw.col)==1){
		bw.col <- rep(bw.col,4)
	}

	#if only one color is given use it for all feature tracks
	if(length(bed.col) == 1){
		bed.col <- rep(bed.col,4)
	}

	#if only one color is given use it for all feature tracks
	if(length(type) == 1){
		type <- rep(type,4)
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
	if(!is.null(exp2)){
		mat2 <- select.subset( exp2$ICE, chrom, start, end, exp2$ABS)
	}

	if(is.null(exp2)){
		mat1$z[mat1$z > cut.off] <- cut.off
		wr <- colorRampPalette(c("white","red"))
		image( mat1, col=wr(256), axes=F, ylim=rev(range(mat1$x)) )
	}else{
		if(coplot == "dual"){
			mat1$z[lower.tri(mat1$z)] <- mat2$z[lower.tri(mat2$z)]
			mat1$z[mat1$z > cut.off] <- cut.off
			wr <- colorRampPalette(c("white","red"))
			image( mat1, col=wr(256), axes=F, ylim=rev(range(mat1$x)) )
		}else{
			mat1$z <- mat2$z - mat1$z
			mat1$z[mat1$z > cut.off] <- cut.off
			mat1$z[mat1$z < -cut.off] <- -cut.off
			bwr <- colorRampPalette(c("blue","white","red"))
			image( mat1, col=bwr(500), axes=F, ylim=rev(range(mat1$x)) )
		}

	}

	#draw pretty axes and boxes
	box(lwd=2)
	size.region <- diff(range(mat1$x))
	if( size.region > 40e6){
		axis(2, at=seq(0,3e9, by=10e6), lab=seq(0,3e9, by=10e6)/1e6, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=10e6), lab=seq(0,3e9, by=10e6)/1e6, lwd=2, cex.axis=1.6)
	}else if( size.region > 2e6){
		axis(2, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6, lwd=2, cex.axis=1.6)
	}else{
		lab <- seq(0,3e9, by=500e3)/1e6
		lab <- sprintf("%.1f",lab)
		axis(2, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=1.6)
		axis(3, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=1.6)
	}

	#fill up empty elements
	if(length(chip) < 4){
		for(i in (length(chip)+1):4){
			chip[i] <- list(NULL)
		}
	}
	#plot the features horizontal
	features( mat1, chrom, genes, chip[[1]], chip[[2]], bed.col[1:2], bw.col[1:2], type=type[1:2] )
	#if the feature entries 3 and 4 are empty, clone 1 and two
	if(length(chip) < 3 || (is.null(chip[[3]]) && is.null(chip[[4]])) ){
		if(!is.null(chip[[1]])){
			chip[[3]] <- chip[[1]]
		}
		if(!is.null(chip[[2]])){
			chip[[4]] <- chip[[2]]
		}
	}
	features( mat1, chrom, genes, chip[[3]], chip[[4]], bed.col[3:4], bw.col[3:4], type=type[3:4], rotate=T )
}

