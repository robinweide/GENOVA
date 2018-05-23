#' visualise.APA.ggplot
#'
#' Plot the APA-results and the differential results.
#'
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param APAlist A list of results from the APA-function.
#' @param title Title-text to plot.
#' @param focus Which sample will be the to-compare sample in the differential row?
#' @param color.fun1 Optional color function for the first row
#' @param color.fun2 Optional color function for the second row
#' @param abs.max Optional Z-max of the first row
#' @param rel.max Optional Z-max of the second row
#' @return A plot
#' @details
#' By substracting the values from each sample with the values from the \code{focus}-sample, we generate the differentials.
#' A positive value in the differential plots thus means an enrichment in that sample versus the \code{focus}-sample.
#' @export
visualise.APA.base <- function( APAlist, title, focus = 1, color.fun1=NULL, color.fun2=NULL, abs.max=NULL, rel.max=NULL, ...){
	#set the margins and store the old settings
	opar <- par(mar=rep(3,4))

	#create a proper plotting device
	numSamples <- length(APAlist)
	layout.matrix <- matrix(1:(numSamples*2), nrow=2, byrow=2)
	layout(layout.matrix)

	#set the standard color functions if none is given by the user
	if(is.null(color.fun1)){
		col.vec <- c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2")
		color.fun1 <- colorRampPalette(rev(col.vec))
	}
	if(is.null(color.fun2)){
		col.vec <- c('#ef8a62','#f7f7f7','#67a9cf')
		color.fun2 <- colorRampPalette(rev(col.vec))
	}

	#set the maximum absolute value if this is not already set
	if(is.null(abs.max)){
		abs.max = 0
		for(i in 1:numSamples){
			if(max(APAlist[[i]]$APAxy$z) > abs.max){
				abs.max = max(APAlist[[i]]$APAxy$z)
			}
		}
	}

	#set the maximum relative value if this is not already set
	if(is.null(rel.max)){
		rel.max = 0
		for(i in 1:numSamples){
			if(max(abs(APAlist[[i]]$APAxy$z-APAlist[[focus]]$APAxy$z)) > rel.max){
				rel.max = max(abs(APAlist[[i]]$APAxy$z-APAlist[[focus]]$APAxy$z))
			}
		}
	}

	#loop over the samples and plot the coverage scores.
	for(i in 1:numSamples){
		max.len <- max(APAlist[[i]]$APAxy$x)
		if(max.len < 1e6){
			APAlist[[i]]$APAxy$x <- APAlist[[i]]$APAxy$x/1e3
			APAlist[[i]]$APAxy$y <- APAlist[[i]]$APAxy$y/1e3
			xlab="distance (kb)"
		}else{
			APAlist[[i]]$APAxy$x <- APAlist[[i]]$APAxy$x/1e6
			APAlist[[i]]$APAxy$y <- APAlist[[i]]$APAxy$y/1e6
			xlab="distance (Mb)"
		}
		abs.mat <- APAlist[[i]]$APAxy
		abs.mat$z[abs.mat$z > abs.max] <- abs.max
		image(abs.mat, col=color.fun1(200), xlab=xlab, mgp=c(2,1,0), zlim=c(0,abs.max))
	}


	#differential scores
	for( i in 1:numSamples){
		diff.mat <- APAlist[[i]]$APAxy
		diff.mat$z <- diff.mat$z - APAlist[[focus]]$APAxy$z
		diff.mat$z[diff.mat$z < -rel.max] <- -rel.max
		diff.mat$z[diff.mat$z >  rel.max] <-  rel.max
		if(i == focus){
			image(diff.mat, axes=F, col='white')
		}else{
			image(diff.mat, col=color.fun2(200), xlab=xlab, zlim=c(-rel.max, rel.max))
		}
	}
	par(opar)
}
