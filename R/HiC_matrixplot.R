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
  sel.bed <- bed[all(bed[,1] == chrom,
                     bed[,2] < end,
                     bed[,3] > start,
                     bed[,6] == '+'),]
  if(nrow(sel.bed)>0){
    col = "red"
    add <- ifelse(sel.bed[,6]=='+', x.wid, -x.wid)
    add.list <- unique(add)
    if(rotate){
      triangle <- lapply( sel.bed[,2], function(x) list(yy=c(x,x,x+add.list),
                                                        xx=c(y1,y2,(y1+y2)/2)))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
    }else{
      #with polygon
      triangle <- lapply( sel.bed[,2], function(x) list(xx=c(x,x,x+add.list),
                                                        yy=c(y1,y2,(y1+y2)/2)))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
    }
  }

  #then repeat for negative
  sel.bed <- bed[all(bed[,1] == chrom,
                     bed[,2] < end,
                     bed[,3] > start,
                     bed[,6] == '-'),]
  if(nrow(sel.bed)>0){
    col="blue"
    add <- ifelse(sel.bed[,6]=='+', x.wid, -x.wid)
    add.list <- unique(add)
    if(rotate){
      triangle <- lapply( sel.bed[,2], function(x) list(yy=c(x,x,x+add.list),
                                                        xx=c(y1,y2,(y1+y2)/2) ))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
    }else{
      #with polygon
      triangle <- lapply( sel.bed[,2], function(x) list(xx=c(x,x,x+add.list),
                                                        yy=c(y1,y2,(y1+y2)/2) ))
      lapply(triangle, function(x) polygon(x$xx, x$yy, col=col, border=col) )
    }
  }

  # check if there are any without orientation and plot these as rectangles
  sel.bed <- bed[all(bed[,1] == chrom,
                     bed[,2] < end,
                     bed[,3] > start,
                     bed[,6] != '-',
                     bed[,6] != '+' ),]
  if(nrow(sel.bed)>0){
    col="black"
    add <- ifelse(sel.bed[,6] != '+', 0, 0)
    add.list <- unique(add)
    if(rotate){
      rect(y1, sel.bed[,2], y2, sel.bed[,3], col=col, border=col)
    }else{
      rect(sel.bed[,2], y1, sel.bed[,3], y2, col=col, border=col)
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

plot.genes <- function( genes, chrom, start, end, y.pos, rotate=F){

  draw_exon( genes, chrom=chrom, y.pos = y.pos, width =0.05, rotate=rotate)
}



features <- function(mat1, chrom, yMax = NULL, genes=NULL, chip1=NULL,
                     chip2=NULL,  autoCHIP = T, rotate=F,
                     type=c("triangle","triangle"), bw.col=NULL, bed.col=NULL){
  # - set bounds for three tracks
  if(is.null(bed.col)){
    bed.col <- 'black'
  }
  if(!is.null(genes)){
    # BW in 1,  bed/null in 2
    if(typeof(chip1) == 'character' && typeof(chip2) != typeof(chip1)){
      gene.pos = 1.05
      chip1.y1 <- 0; chip1.y2 <- 0.5
      chip2.y1 <- 0.6; chip2.y2 <- 0.8
      # BW in 2,  bed/null in 1
    } else if(typeof(chip2) == 'character' && typeof(chip2) != typeof(chip1)){
      gene.pos = 1.05
      chip1.y1 <- 0; chip1.y2 <- 0.2
      chip2.y1 <- 0.3; chip2.y2 <- 0.8
    }  else {
      gene.pos = 1.05
      chip1.y1 <- 0; chip1.y2 <- 0.375
      chip2.y1 <- 0.425; chip2.y2 <- 0.8
    }

  } else {
    # BW in 1,  bed/null in 2
    if(typeof(chip1) == 'character' && typeof(chip2) != typeof(chip1)){
      gene.pos = 0
      chip1.y1 <- 0; chip1.y2 <- 0.675
      chip2.y1 <- 0.725; chip2.y2 <- 1.1
      # BW in 2,  bed/null in 1
    } else if(typeof(chip2) == 'character' && typeof(chip2) != typeof(chip1)){
      gene.pos = 0
      chip1.y1 <- 0; chip1.y2 <- 0.375
      chip2.y1 <- 0.425; chip2.y2 <- 1.1
    }  else {
      gene.pos = 0
      chip1.y1 <- 0; chip1.y2 <- 0.5
      chip2.y1 <- 0.7; chip2.y2 <- 1.1
    }


  }
  y.values <- list(gene.pos=gene.pos, chip1.y1=chip1.y1, chip1.y2=chip1.y2,
                   chip2.y1=chip2.y1, chip2.y2=chip2.y2)


  # make chip smaller

  if(max(mat1$x)-min(mat1$x) > 2e6){

    if(typeof(chip1) == 'list'){
      newSize = diff(unlist(y.values[2:3])) * 0.75
      D = diff(unlist(y.values[2:3])) - newSize
      #y.values[[2]] = y.values[[2]] + D/2
      y.values[[3]] = y.values[[3]] - D
    }

    if(typeof(chip2) == 'list'){
      newSize = diff(unlist(y.values[4:5])) * 0.75
      D = diff(unlist(y.values[4:5])) - newSize
      y.values[[4]] = y.values[[4]] + D
      #y.values[[5]] = y.values[[5]] - D/2
    }
  }



  # - plot an empty canvas c(0,1.3)

  if(rotate){
    plot(0, type='n', ylim=rev(range(mat1$x)), xlim=c(1.3,0), axes=F,
         xlab="", ylab="" )
  }else{
    plot(0, type='n', xlim=range(mat1$x), ylim=c(0,1.3), axes=F,
         xlab="", ylab="" )
  }


  # - plot inner track c(0.1,0.4) -> chipInner

  if(typeof(chip1) == 'list'){ # BED!

    if(ncol(chip1) < 6){
      type[1] = "rectangle"
      # there is no orientation in column 6
    } else if(autoCHIP & !any(c("-","+") %in% chip1[,6])){
      type[1] = "rectangle"
    }

    if(type[1]=="triangle"){
      # +
      tmp_step = (y.values$chip1.y2 - y.values$chip1.y1)/2
      tmp_step = tmp_step * 0.9
      blabla = plot.triangle( chip1[chip1[,6] == '+',  ], chrom=chrom,
                              y1=y.values$chip1.y1,
                              y2=y.values$chip1.y1 + tmp_step,
                              start=min(mat1$x), end=max(mat1$x), rotate=rotate)
      # -
      blabla = plot.triangle( chip1[chip1[,6] == '-',  ], chrom=chrom,
                              y1=y.values$chip1.y2 - tmp_step,
                              y2=y.values$chip1.y2 , start=min(mat1$x),
                              end=max(mat1$x), rotate=rotate)
      # # other
      blabla = plot.rectangle(chip1[chip1[,6] != '-' & chip1[,6] != '+',],
                              chrom=chrom, y1=y.values$chip1.y1,
                              y2=y.values$chip1.y2, start=min(mat1$x),
                              end=max(mat1$x), col=bed.col[1], rotate=rotate)

    }else if(type[1]=="rectangle"){
      blabla = plot.rectangle(chip1, chrom=chrom, y1=y.values$chip1.y1,
                              y2=y.values$chip1.y2, start=min(mat1$x),
                              end=max(mat1$x), col=bed.col[1], rotate=rotate)
    }

  } else if(typeof(chip1) == 'character'){ # BW!

    blabla =  suppressWarnings(plot.bw(chip1, chrom,  min(mat1$x),
                                       max(mat1$x), y.values$chip1.y1,
                                       y.values$chip1.y2, col=bw.col[1],
                                       rotate=rotate, yMax = yMax))

  }

  # - plot middle track c(0.5,0.8) -> chipOuter

  if(typeof(chip2) == 'list'){ # BED!

    if(ncol(chip2) < 6){
      type[2] = "rectangle"
      # there is no orientatio in column 6
    } else if(autoCHIP & !any(c("-","+") %in% chip2[,6])){
      type[2] = "rectangle"
    }

    if(type[2]=="triangle"){
      # +
      tmp_step = (y.values$chip2.y2 - y.values$chip2.y1)/2
      tmp_step = tmp_step * 0.9
      blabla = plot.triangle( chip2[chip2[,6] == '+',  ], chrom=chrom,
                              y1=y.values$chip2.y1,
                              y2=y.values$chip2.y1 + tmp_step,
                              start=min(mat1$x), end=max(mat1$x), rotate=rotate)
      # -
      blabla = plot.triangle( chip2[chip2[,6] == '-',  ], chrom=chrom,
                              y1=y.values$chip2.y2 - tmp_step,
                              y2=y.values$chip2.y2 , start=min(mat1$x),
                              end=max(mat1$x), rotate=rotate)
      # # other
      blabla = plot.rectangle(chip2[chip2[,6] != '-' & chip2[,6] != '+',],
                              chrom=chrom, y1=y.values$chip2.y1,
                              y2=y.values$chip2.y2, start=min(mat1$x),
                              end=max(mat1$x), col='black', rotate=rotate)

    }else if(type[2]=="rectangle"){
      blabla =  plot.rectangle(chip2, chrom=chrom, y1=y.values$chip2.y1,
                               y2=y.values$chip2.y2, start=min(mat1$x),
                               end=max(mat1$x), col=bed.col[2], rotate=rotate)
    }

  } else if(typeof(chip2) == 'character'){ # BW!

    blabla = suppressWarnings(plot.bw(chip2, chrom, min(mat1$x), max(mat1$x),
                                      y.values$chip2.y1, y.values$chip2.y2,
                                      col=bw.col[2], rotate=rotate,
                                      yMax = yMax))

  }


  # - plot outer track c(0.9,1.2) -> genes


  if(!is.null(genes)){
    blabla = plot.genes(genes, chrom, min(mat1$x), max(mat1$x),
                        y.pos=y.values$gene.pos, rotate=rotate)
  }


  # - if rotate, plot an empty canvas for quadrant 4

  if(rotate){
    plot(0, type='n', xlim=range(mat1$x), ylim=c(0,1.3), axes=F,
         xlab="", ylab="" )
  }

}





plot.bw <- function( file, chrom, start, end, y1, y2, col,
                     yMax = NULL,rotate=F){
  if(require("bigwrig") == F){
    stop("Please install github.com/jayhesselberth/bigwrig\n")
    }

  #leave here in case something changes
  #official call but don't need it
  #d <- read_bigwig( file, chrom=chrom, start=start, end=end)
  #d <- as.data.frame(d)

  d <- bigwrig::read_bigwig_impl( file, chrom=chrom, start=start, end=end)

  #to speed up plotting we are going to collapse the same consecutive
  #values
  i <- which(diff(d[,4])!=0) #select points where the consecutive values differ
  #select the middle of a set of consecutive values
  sel <- floor(c(0,head(i,-1))+i)/2
  d <- d[sel,]


  max.val <- max(d[,4])
  d[,4] = as.numeric(d[,4])
  if(!is.null(yMax)){
    d[d[,4] > yMax,] = yMax
    max.val = yMax
  }

  y.range <- abs(y1-y2)
  y.val <- y.range*d[,4]/max.val
  if(rotate){
    segments(y1, d[,2], y1+y.val, d[,2], col=col)
  }else{
    segments(d[,2], y1, d[,2], y1+y.val, col=col)
  }
}

#overlay TAD positions with the hi-c data
draw.tads <- function( tads, chrom, tads.type="lower", tads.color ="#006837",
                       lwd=2){

  # check if tads is a list of dataframes or a data.frame
  if(inherits(tads, "list") ){ #this seems to be a tad!
    if(all(unlist(lapply(tads, inherits, "data.frame")))){ # all DFs in list!
      tadList = T
    } else {
      stop("tads is not a (list of) data.frame!")
    }
  } else if(inherits(tads, "data.frame")) {
    tmpList = list() # cool. there is one df of tads. Put in a list
    tmpList[[1]] = tads
    tads = tmpList
    rm(tmpList)
  } else {
    stop("tads is not a (list of) data.frame!")
  }
  # if you survived this, you will now have a list of dfs!

  # get the tads color, resize and type
  if(length(tads) != length(tads.type)){
    tads.type = rep(tads.type[1], length(tads))
  }

  if(length(tads) != length(tads.color)){
    tads.color = rep(tads.color[1], length(tads))
  }


  for(listIDX in 1:length(tads)){

    tad = tads[[listIDX]]

    # tads has min. 6 cols : 1 and 4 are the chrom. 2 is outer-edge 5' anchor.
    # 6 is outer edge 3' anchor
    tad <- tad[tad[,1]==chrom,]

    if(tads.type[listIDX] == "both"){
      rect(tad[,2], tad[,2], tad[,3], tad[,3], border=tads.color[listIDX],
           lwd=lwd)
    }else if(tads.type[listIDX] == "lower"){
      segments(tad[,2], tad[,2], tad[,2], tad[,3], col=tads.color[listIDX],
               lwd=lwd)
      segments(tad[,2], tad[,3], tad[,3], tad[,3], col=tads.color[listIDX],
               lwd=lwd)
    }else if(tads.type[listIDX] == "upper"){
      segments(tad[,2], tad[,2], tad[,3], tad[,2], col=tads.color[listIDX],
               lwd=lwd)
      segments(tad[,3], tad[,2], tad[,3], tad[,3], col=tads.color[listIDX],
               lwd=lwd)
    }else{
      stop("Wrong option for TAD plot type: upper, lower and both are allowed")
    }

  }


}

#overlay loop positions with the hi-c data
draw.loops <- function( loops, chrom, loops.type="both", loops.color = "green",
                        lwd=2 , loops.resize ){

  # check if loops is a list of dataframes or a data.frame
  if(inherits(loops, "list") ){ #this seems to be a loop!
    if(all(unlist(lapply(loops, inherits, "data.frame")))){ # all DFs in list!
      loopList = T
    } else {
      stop("Loops is not a (list of) data.frame!")
    }
  } else if(inherits(loops, "data.frame")) {
    tmpList = list() # cool. there is one df of loops. Put in a list
    tmpList[[1]] = loops
    loops = tmpList
    rm(tmpList)
  } else {
    stop("Loops is not a (list of) data.frame!")
  }
  # if you survived this, you will now have a list of dfs!

  # get the loops color, resize and type
  if(length(loops) != length(loops.type)){
    loops.type = rep(loops.type[1], length(loops))
  }

  if(length(loops) != length(loops.color)){
    loops.color = rep(loops.color[1], length(loops))
  }

  if(length(loops) != length(loops.resize)){
    loops.resize = rep(loops.resize[1], length(loops))
  }

  # loops over de list and plot
  for(listIDX in 1:length(loops)){

    loop = loops[[listIDX]]

    # loops has min. 6 cols : 1 and 4 are the chrom. 2 is outer-edge 5' anchor.
    # 6 is outer edge 3' anchor
    loop <- loop[loop[,1]==chrom,]

    # resize loops
    loop[,2] <- loop[,2] - loops.resize[listIDX]
    loop[,3] <- loop[,3] + loops.resize[listIDX]
    loop[,5] <- loop[,5] - loops.resize[listIDX]
    loop[,6] <- loop[,6] + loops.resize[listIDX]

    if(loops.type[listIDX] == "both"){
      rect(xleft = loop[,5], ybottom = loop[,3], xright = loop[,6],
           ytop = loop[,2], border=loops.color[listIDX], lwd=lwd)
      rect(xleft = loop[,2], ybottom = loop[,6], xright = loop[,3],
           ytop = loop[,5], border=loops.color[listIDX], lwd=lwd)
    }else if(loops.type[listIDX] == "lower"){

      rect(xleft = loop[,2], ybottom = loop[,6], xright = loop[,3],
           ytop = loop[,5], border=loops.color[listIDX], lwd=lwd)

    }else if(loops.type[listIDX] == "upper"){
      rect(xleft = loop[,5], ybottom = loop[,3], xright = loop[,6],
           ytop = loop[,2], border=loops.color[listIDX], lwd=lwd)
    }else{
      stop("Wrong option for loop plot type: upper, lower and both are allowed")
    }

  }

}




#' hic.matrixplot
#'
#' Plot a matrix (or two) for a region of interest with annotations
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @author Elzo de Wit, \email{e.d.wit@nki.nl}
#' @param exp1 The Hi-C experiment object: produced by construct.experiment().
#' @param exp2 Optional: a Hi-C experiment object.
#' @param chrom Chromosome of the region of interest
#' @param start Start position of the region of interest
#' @param end End position of the region of interest
#' @param cut.off The cut.off for the hic-matrix plot, in the diff option the
#' negative of this is the lower bound
#' @param chip A list of feature tracks, can be bed structure
#' (i.e. data frames) or a path to bigwig file (i.e. character variable),
#' maximum length of 4. Placement is inner-top, outer-top, outer-left,
#' inner-left.
#' @param inferno White/Red/black or White/Red coloscale?
#' @param cexTicks Change size of numbers on the axis
#' @param bed.col Color of the bed track (max.len is 4)
#' @param bw.col Same as bed col, but for the bigwig track
#' @param yMax Set the maximum height of all biwigs.
#' @param type Should a rectangle or a triangle be drawn?
#' Note that for a triangle a 6th strand column should be included
#' @param guessType If an element in the chip is a dataframe,
#' infer the type from column 6? If true, column six should hold "+" and "-".
#' If a row has other characters, we view this entry as no-oriention and
#' thus plot a rectangle.
#' @param coplot When drawing together two experiments: [dual] is bottom
#' triangle exp2, top triangle exp1; [diff] plots a substraction of exp2-exp1
#' @param genes Structure with gene information,
#' will only be combined with bed structure
#' @param tads BED-like dataframe or a list of these data.frames
#' @param tads.type How to show TADS: upper, lower and or both
#' @param loops BED-like dataframe or a list of these data.frames
#' @param loops.type How to show loops: upper, lower and or both
#' @param loops.resize Make the loop-boxes bigger by X bp. Can help visibility.
#' @param loops.color Which color do you want the loops to have?
#' @param skipAnn Do not plot outside-annotation. Can be used to plot other
#' things next to the matrix.
#' @param symmAnn Put features 1&2 also on verical (ignore chip-entries 3&4)
#' @param check.genome Check if reference genomes in exp1 and exp2 are the same
#' @param antoni Logical: plot an explorer of the microscopic world
#' @note
#' To plot genes, a gene-model data.frame must be made. This can be done via a
#' multitude of ways (e.g. biomart, UCSC table browser). The resulting
#' data.frame must have one exon per row with the following columns:
#' geneID | transcriptID | chrom | txStart | txEnd | exonStart | exonEnd |
#' strand.
#' Alternatively, if you download the knowngene table from UCSC, you can
#' directly use this table (with exons combined per row), by renaming
#' exonStarts and exonEnds to exonStart and exonEnd.
#' @examples
#' # plot two matrices of Hap-1 cells, including their respective loop-calls
#' hic.matrixplot(exp1 = Hap1_Haarhuis2017_10kb,
#'                exp2 = Hap1_Sanborn2015_10kb,
#'                chrom = 'chr7',
#'                start = 25e6,
#'                end=30e6,
#'                loops = list(Haarhuis2017_Loops, sanborn2015_Loops),
#'                loops.color = c('blue','green'),
#'                loops.type = c('upper','lower'),
#'                loops.resize = c(20e3,20e3), # expand for visibility
#'                cut.off = 25) # upper limit of contacts
#' @return A matrix-plot
#' @export
hic.matrixplot <- function( exp1, exp2=NULL, chrom, start, end, cut.off=NULL,
                            chip=list(NULL, NULL, NULL, NULL), inferno = T,
                            cexTicks = 1, bed.col=rep("black",4),
                            bw.col="black", yMax = NULL,
                            type=rep("triangle",4), guessType = T,
                            coplot="dual", genes=NULL,
                            tads=NULL, tads.type="lower", loops=NULL,
                            loops.type="lower", tads.color = "#3288bd",
                            loops.resize = 0, loops.color = "#3288bd",
                            skipAnn = F, symmAnn = F,
                            check.genome = T, antoni = F){

  if(length(chip) < 3){
    symmAnn = T
  }
  #some error handling
  if(!is.null(exp2)){
    if(any(exp2$RMCHROM, exp1$RMCHROM)){
      check.genome = F
    }
    #make sure the resolutions are the same
    if(exp1$RES != exp2$RES){
      stop("The Hi-C matrices should have the same resolution")
    }

    if(check.genome){
      if(!all(exp1$ABS[,4] == exp2$ABS[,4])){
        msg = paste0("Not all ICE indexes are the same.\nAre you these ",
                     "experiments were mapped to the same genome (-build)?")
        stop(msg)
      }
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
  mat1 <- select.subset( exp1, chrom, start, end)
  if(!is.null(exp2)){
    mat2 <- select.subset( exp2, chrom, start, end)
  }

  if(antoni){
    load(url("https://github.com/robinweide/GENOVA/raw/dev/R/antoni.Rdat"))
    CR = colorRampPalette(c('#ffffff','#f0f0f0','#d9d9d9',
                            '#bdbdbd','#969696','#737373','#525252',
                            '#252525','#000000'))

    image(  antoni[,nrow(antoni):1 ], col = rev(CR(10)),
            axes=F,ylim=rev(range(antoni)))
  } else if(is.null(exp2)){
    if(is.null(cut.off)){
      cut.off = max(quantile(mat1$z, .99))
      warning("No cut.off was given: using 99% percentile: ",
              round(cut.off), ".")
    }
    mat1$z[mat1$z > cut.off] <- cut.off
    if(inferno){
      higlassCol <- c('white', '#f5a623', '#d0021b', 'black')
      wr <- colorRampPalette(higlassCol)

      image( mat1, col=wr(1e4), axes=F, ylim=rev(range(mat1$x)) )

    } else {
      wr <- colorRampPalette(c("white","red"))
      image( mat1, col=wr(1e4), axes=F, ylim=rev(range(mat1$x)) )
    }

  } else{
    if(coplot == "dual"){

      mat1$z[upper.tri(mat1$z)] <- mat2$z[upper.tri(mat2$z)]

      if(is.null(cut.off)){
        cut.off = max(quantile(mat1$z, .995))
        warning("No cut.off was given: using 99.5% percentile: ",
                round(cut.off), ".")
      }
      mat1$z[mat1$z > cut.off] <- cut.off


      if(inferno){
        higlassCol <- c('white', '#f5a623', '#d0021b', 'black')
        wr <- colorRampPalette(higlassCol)
        image( mat1, col=wr(1e4), axes=F, ylim=rev(range(mat1$x)) )
      } else {
        wr <- colorRampPalette(c("white","red"))
        image( mat1, col=wr(1e4), axes=F, ylim=rev(range(mat1$x)) )
      }

    }else{
      mat1$z <- mat2$z - mat1$z
      if(is.null(cut.off)){
        cut.off = max(quantile(mat1$z, .995))
        warning("No cut.off was given: using 99.5% percentile: ",
                round(cut.off), ".")
      }
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
    axis(2, at=seq(0,3e9, by=10e6), lab=seq(0,3e9, by=10e6)/1e6,
         lwd=2, cex.axis=cexTicks)
    axis(3, at=seq(0,3e9, by=10e6), lab=seq(0,3e9, by=10e6)/1e6,
         lwd=2, cex.axis=cexTicks)
  }else if( size.region > 2e6){
    axis(2, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6,
         lwd=2, cex.axis=cexTicks)
    axis(3, at=seq(0,3e9, by=1e6), lab=seq(0,3e9, by=1e6)/1e6,
         lwd=2, cex.axis=cexTicks)
  }else{
    lab <- seq(0,3e9, by=500e3)/1e6
    lab <- sprintf("%.1f",lab)
    axis(2, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=cexTicks)
    axis(3, at=seq(0,3e9, by=500e3), lab=lab, lwd=2, cex.axis=cexTicks)
  }

  #draw tads on the image plot
  if(!is.null(tads)){
    draw.tads( tads, chrom, tads.type=tads.type, tads.color = tads.color)
  }

  #draw loops on the image plot
  if(!is.null(loops)){
    draw.loops( loops, chrom, loops.type=loops.type,
                loops.resize = loops.resize, loops.color = loops.color)
  }

  #fill up empty elements
  if(length(chip) < 4){
    for(i in (length(chip)+1):4){
      chip[i] <- list(NULL)
    }
  }
  if(!skipAnn){

    #plot the features horizontal
    features( mat1, chrom, genes, chip1 = chip[[1]], chip2 = chip[[2]],
              autoCHIP = guessType ,yMax = yMax,bed.col = bed.col[1:2],
              bw.col = bw.col[1:2], type=type[1:2] )
    # clone 1 and two to feature entries 3 and 4
    if(symmAnn){
      if(!is.null(chip[[1]])){
        chip[[3]] <- chip[[1]]
      }
      if(!is.null(chip[[2]])){
        chip[[4]] <- chip[[2]]
      }
    }

    features( mat1, chrom, genes, chip[[3]], chip[[4]], bed.col = bed.col[3:4],
              autoCHIP = guessType ,yMax = yMax,bw.col = bw.col[3:4],
              type=type[3:4], rotate=T )
  }
}
