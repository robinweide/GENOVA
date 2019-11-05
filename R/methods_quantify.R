
# Documentation -----------------------------------------------------------

#' @title Quantification of results
#' @name quantify
#' @description Here be a placholder description
#' @param discovery A \code{discovery} object as returned by GENOVA analysis functions.
#' @param signal_size The width/height of the signal (e.g. a value of 3 will 
#' @param ... further arguments passed to or from other methods.
#' take the middle 3x3 matrix of the APA).

# Functions ---------------------------------------------------------------

# If you fail it is best to fail graciously

#' @export
#' @rdname quantify
#' @usage NULL
quantify.default <- function(discovery, ...) {
  stop("No quantify method for class '", class(discovery),
       "' has been implemented.", call. = FALSE)
}


#' @rdname quantify
#' @export
quantify.APA_discovery <- function(discovery, signal_size = 3, ...) {

  mid = unique(floor(dim(discovery$signal[,,1])/2)+1)
  mid_range = (mid - ((signal_size - 1)/2)):(mid + ((signal_size - 1)/2))
  
  median_donuthole_lfc = apply(discovery$signal[mid_range,mid_range,] ,3, median)
  median_donut_lfc = apply(discovery$signal[-mid_range,-mid_range,] ,3, median)
  
  within_sample_lfc = log2(median_donuthole_lfc/median_donut_lfc)
  within_sample_lfc = data.table(as.data.frame(within_sample_lfc),keep.rownames = T)

  eg = t(combn(1:nrow(within_sample_lfc), 2))
  out = cbind(within_sample_lfc[eg[,1],], 
              within_sample_lfc[eg[,2],])
  
  colnames(out) = c('exp1', 'exp1_lfc', 'exp2', 'exp2_lfc')
  out$contrast_lfc = log2(out$exp2_lfc/out$exp1_lfc)

  
  return(out[,c(1,3,2,4,5)])
  
}

#' @rdname visualise
#' @export
quantify.saddle_discovery <- function(x, ...){
  dat <- x$saddle
  
  # get bins
  MAXbin <- max(dat$q1)
  binsTOse = floor(MAXbin * .2)
  binsTOse = max(1, binsTOse)
  
  # transform from log2-space
  dat$unLog = 2 ** dat$mean
  
  # assign cat
  dat$CC <- 'XXX'
  dat[dat$q1 <= binsTOse & dat$q2 <= binsTOse,"CC"] = "BB"
  dat[dat$q1 >= MAXbin-binsTOse+1 & dat$q2 >= MAXbin-binsTOse+1,"CC"] = "AA"
  dat[dat$q1 <= binsTOse & dat$q2 >= MAXbin-binsTOse+1,"CC"] = "AB"
  dat <- dat[! dat$CC == 'XXX',]
  
  # summarise
  dat[ , score := mean(unLog), by = c("exp", "chr", "CC")]
  dat <- unique(dat[,-c(3:6)])
  
  # compute strength
  SPLT <- split(dat, list(dat$exp, dat$chr))
  SPLT <- SPLT[lapply(SPLT, nrow) == 3]
  SPLT <- lapply(SPLT, function(y){
    # check if all exist
    
    data.frame(exp = y[1,1], 
               chr = y[1,2], 
               strength = log2( 
                 (y[y$CC == 'AA',score] * y[y$CC == 'BB',score]) / 
                   (y[y$CC == 'AB',score]**2) 
               )) 
  })
  strength <- data.table::rbindlist(SPLT)
  dat <- merge(dat, strength, by = c('exp', 'chr'))
  dat$exp <- factor(dat$exp, levels = unique(x$saddle$exp))
  dat$CC  <- factor(dat$CC, levels = c('AA', 'BB', 'AB'))
  # maybe add class?
  return(dat)
}