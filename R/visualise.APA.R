visualise.APA <- function( APAlist, title, zCoef = 1, focus = 1, color.fun1=NULL, color.fun2=NULL, ...){
  #set the margins and store the old settings
  opar <- par(mar=rep(3,4))

  #create a proper plotting device
  numSamples <- length(APAlist)
  layout.matrix <- matrix(1:(numSamples*2), nrow=2, byrow=2)
  layout(layout.matrix)


  if(is.null(color.fun1)){
    col.vec <- c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2")
    color.fun1 <- colorRampPalette(rev(col.vec))
  }
  if(is.null(color.fun2)){
    col.vec <- c('#ef8a62','#f7f7f7','#67a9cf')
    color.fun2 <- colorRampPalette(rev(col.vec))
  }


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
  }



  #absolute scores
  for( i in 1:numSamples){
    image(APAlist[[i]]$APAxy, col=color.fun1(200), xlab=xlab, mgp=c(2,1,0))
  }



  #differential scores
  for( i in 1:numSamples){
    diff.mat <- APAlist[[i]]$APAxy
    diff.mat$z <- diff.mat$z - APAlist[[focus]]$APAxy$z
    if(i == focus){
      image(diff.mat, axes=F, col='white')
    }else{
      image(diff.mat, col=color.fun2(200), xlab=xlab)
    } 
  }
  par(opar)
} 
