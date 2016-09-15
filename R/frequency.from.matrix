#' Plot pyramid-style HiC-regions.
#'
#' Extracts a sub-matrix from given coordinates and plots them pyramid-style (i.e. the Hi-C diagonal becomes horizontal).
#'
#' @param experiment The Hi-C experiment object of a sample: produced by construct.experiment().
#' @param chrom Chromosome.
#' @param start Start position in bp.
#' @param end End position in bp.
#' @param window 
#' @param q.val
#' @param cut.off
#' @param extraPadding
#' @param ylim The y-start and -stop: zero is at the diagonal,
#' @return A plot.
#' @import data.table
frequency.from.matrix <- function( experiment, chrom, start, end, window = 10e3, q.val=0.95, cut.off=0, extraPadding = 1.75, ylim = c(0, 250), shinyAxis = T, ... ){
  # Check size: if less than 50 bins: makes no sense...
  resolution <- experiment$RES
  mat <- select.subset(experiment$ICE, chrom, start- 2e7, end+ 2e7, experiment$ABS  )
  pos <- which(mat$z >0 , arr.ind=T)
  mat.start <- min(mat$x)
  cnt <- as.vector(mat$z[mat$z > 0])
  cnt <- cnt[pos[,2]-pos[,1] > 0]
  pos <- pos[pos[,2]-pos[,1] > 0,]
  if(cut.off == 0){
    cut.off <- quantile(cnt, q.val)
  }       
  cnt[cnt > cut.off ] <- cut.off
  cnt <- cnt/cut.off
  num.rect <- length(unique(pos[,1]))
  cex = (400/num.rect)*2.5
  if(shinyAxis){
    plot( mat.start+window*(pos[,1]+pos[,2])/2, pos[,2]-pos[,1], pch=18, col=rgb(1, 1-cnt, 1-cnt), xlab="", ylab="", ylim = ylim,cex=cex, xlim=c(start,end), axes=F, ...); box(lwd=1); axis(1, at=seq(0,3e9, by=5e5), lab=NA, lwd=1)
  } else {
    plot( mat.start+window*(pos[,1]+pos[,2])/2, pos[,2]-pos[,1], pch=18, col=rgb(1, 1-cnt, 1-cnt), xlab="", ylab="", ylim = ylim,cex=cex, xlim=c(start,end), ...)
  }
  #gc()
} 

# 
# wt <- construct.experiment(ICEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/WT_10000_iced.matrix", BEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/WaplKO_combined_10000_abs.bed", SAMPLENAME = "WT", COLOR = 'black')
# wapl <- construct.experiment(ICEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/WaplKO_combined_10000_iced.matrix", BEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/WaplKO_combined_10000_abs.bed", SAMPLENAME = "WT", COLOR = 'red')
# ssc<- construct.experiment(ICEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/SCC4KO_10000_iced.matrix", BEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/WaplKO_combined_10000_abs.bed", SAMPLENAME = "WT", COLOR = 'green')
# dko <- construct.experiment(ICEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/DKO_10000_iced.matrix", BEDpath = "~/dataDisk/Projects/WAPL/finalDataSet/mats/WaplKO_combined_10000_abs.bed", SAMPLENAME = "WT", COLOR = 'blue')
# 
# chrom <- "chr1"
# start <- 237180000
# end <- 241750000
# 
# par(mfrow=c(1,1), xaxs="i", yaxs="i", mar=rep(1,4))
# frequency.from.matrix(experiment = wt, chrom = chrom,start = start, end = end, ylim = c(0,300),cut.off = 7)
# frequency.from.matrix(experiment = wapl, chrom = chrom,start = start, end = end, ylim = c(0,300),cut.off = 7)
# frequency.from.matrix(experiment = ssc, chrom = chrom,start = start, end = end, ylim = c(0,300),cut.off = 7)
# frequency.from.matrix(experiment = dko, chrom = chrom,start = start, end = end, ylim = c(0,300),cut.off = 7)
