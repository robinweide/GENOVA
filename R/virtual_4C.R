VP = HICCUPS[2498,1:3]
xlim      <- VP[,2:3]
xlim[1]   <- xlim[1]-5e6
xlim[2]   <- xlim[2]+5e6

virtual_4C(exp = WT_new, viewpoint = VP, xlim = xlim)


virtual_4C <- function(exp, viewpoint, xlim = NULL){
  
  vp_idx <- median(bed2idx(exp$IDX, viewpoint))
  
  signal <- NULL
  if( is.null(xlim) ){
    # run genome-wide ==========================================================
    signal <- exp$MAT[V1 == vp_idx | V2 == vp_idx]
    
    upstream_signal   <- signal[V2 %in% vp_idx][,c(1,3)]
    downstream_signal <- signal[V1 %in% vp_idx][,2:3]
  } else {
    # run for a region =========================================================
    flank_downstream <- floor((viewpoint[,2] - xlim[1])/attr(exp, 'resolution'))
    flank_upstream <- floor((xlim[2] - viewpoint[,3])/attr(exp, 'resolution'))
    
    range_idx <- unlist(vp_idx - flank_upstream):unlist(vp_idx + flank_upstream)
    
    upstream_signal <- exp$MAT[list(range_idx,vp_idx), nomatch = 0][,c(1,3)]
    downstream_signal <- exp$MAT[list(vp_idx, range_idx), nomatch = 0][,2:3]
  }
  
  signal <- rbind(upstream_signal, downstream_signal, use.names = FALSE)
  colnames(signal) = c('V4','signal')
  setkey(signal, 'V4')
  
  # set basepairs --------------------------------------------------------------
  bed_ranges <- exp$IDX[match(signal$V4, exp$IDX$V4)]
  setkey(bed_ranges, 'V4')
  signal <- signal[bed_ranges]
  signal$mid = rowMeans(signal[,4:5])
  
  if( !is.null(xlim) ){
    signal <- signal[V1 == viewpoint[1,1]]
  }
 
  
  # output ---------------------------------------------------------------------
  signal <- signal[,c(3,6,2)]
  colnames(signal) <- c('chromosome','mid','signal')
  
  signal <- structure(signal, 
                      'viewpoint' = viewpoint, 
                      'xlim' = xlim,
                      'sample' = attr(exp, 'sample'),
                      class = "virtual4C_discovery",
                      package = "GENOVA")
  
  signal
}



visualise.virtual4C_discovery(){
  
  
}


library(ggplot2)
bed = HICCUPS[,1:3]
bins = 200
blackout_region <- 2.5e5
data <- signal
draw_blackout <- T

if( is.null(bins) ) {
  bins = nrow(data)
}

if( !is.null(bed)) {
  bed = bed[bed[,1] == attr(data, 'viewpoint_chromosome'),2:3]
}
  
blackout_up   <- attr(data, 'viewpoint_start') - (blackout_region/2)
blackout_down <- attr(data, 'viewpoint_end') + (blackout_region/2)

data <- data[!(data$mid > blackout_up & data$mid < blackout_down)]

breaks <- seq(min(data$mid), max(data$mid), length.out = bins)
smooth <- data[, mean(signal),by = findInterval(data$mid, breaks)]
smooth$mid = breaks[unlist(smooth[,1])] + unique(diff(breaks)/2)
smooth[,1] = NULL
colnames(smooth) = c("signal","mid")

p = ggplot(data, aes(x= mid/1e6, y = signal)) +
  annotate('rect', 
           fill = "black",
           xmin = bed[,1]/1e6, 
           xmax = bed[,2]/1e6,
           ymin = -ceiling(max(smooth$signal))/100,
           ymax = 0)  +
  geom_col(data = smooth, fill = 'black') +
  coord_cartesian(xlim = unlist(attr(signal, 'xlim'))/1e6 ,expand = F) +
  theme_classic() 

if( draw_blackout ){
  p = p + annotate('rect',
                   fill =  "#D8D8D8",
                   xmin = blackout_up/1e6,
                   xmax = blackout_down/1e6,
                   ymin = 0,
                   ymax = ceiling(max(smooth$signal))) 
}
 
p








