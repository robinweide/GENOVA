
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
  out_summarised = cbind(within_sample_lfc[eg[,1],], 
              within_sample_lfc[eg[,2],])
  
  colnames(out_summarised) = c('exp1', 'exp1_lfc', 'exp2', 'exp2_lfc')
  out_summarised$contrast_lfc = log2(out_summarised$exp2_lfc/out_summarised$exp1_lfc)

  out_summarised = out_summarised[,c(1,3,2,4,5)]

  raw_data_melt <- lapply(discovery$signal_raw, function(dat) {
    idx <- seq_along(dim(dat))
    idx <- setNames(idx, paste0("Var", idx))
    idx <- lapply(idx, function(i) {
      dimnames(dat)[[i]][as.vector(slice.index(dat, i))]
    })
    idx[c(2,3)] <- lapply(idx[c(2,3)], as.numeric)
    as.data.table(idx)[, value := as.vector(dat)]
  })
  raw_data_melt <- data.table::rbindlist(raw_data_melt, idcol = 'sample_name')
  
  bp_bins <- unique(raw_data_melt$Var2)
  pixel_range <- which(bp_bins == 0)
  pixel_range <- c(pixel_range - 1, pixel_range, pixel_range+1)
  bp_bins <- bp_bins[pixel_range]
  
  median_donuthole_lfc <- raw_data_melt[Var2 %in% bp_bins & Var3 %in% bp_bins, 
                                        median(value, na.rm = T), 
                                        by = 'sample_name,Var1']
  median_donut_lfc <- raw_data_melt[(!Var2 %in% bp_bins) & (!Var3 %in% bp_bins), 
                                    median(value, na.rm = T), 
                                    by = 'sample_name,Var1']
  median_lfc <- merge.data.table(median_donuthole_lfc,
                                 median_donut_lfc, 
                                 by = c('sample_name','Var1'))
  colnames(median_lfc)[2:4] <- c('loop_ID', 'pixel', 'background')
  
  out <- median_lfc[ , list('difference' = pixel/background,
                            'lfc' = log2(pixel/background)), 
                     by = 'sample_name,loop_ID']
  out$sample_name <- factor(out$sample_name, levels = names(discovery$signal_raw))
  
  bed <- as.data.table(tstrsplit(out$loop_ID, "\\:|\\-|;"))
  bed <- bed[, list("chr_up" = V1, "start_up" = as.numeric(V2), 
                    "end_up" = as.numeric(V3),
                    "chr_down" = V4, "start_down" = as.numeric(V5), 
                    "end_down" = as.numeric(V6))]
  bed$sample_name <- out$sample_name
  bed$difference <- out$difference
  bed$lfc <- out$lfc
  
  return(list('per_sample' = out_summarised, 'per_loop' = bed)) 
}

#' @rdname quantify
#' @export
quantify.saddle_discovery <- function(discovery, ...){
  dat <- discovery$saddle
  
  # get bins
  MAXbin <- max(dat$q1, na.rm = T)
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
  dat$exp <- factor(dat$exp, levels = unique(discovery$saddle$exp))
  dat$CC  <- factor(dat$CC, levels = c('AA', 'BB', 'AB'))
  # maybe add class?
  return(dat)
}

#' @rdname quantify
#' @export
quantify.IIT_discovery <- function(discovery, ...) {
  data <- discovery$results
  data <- melt.data.table(data, id.vars = c("x", "y"))
  data$distance <- data$y - data$x
  
  summ <- data[, as.list(summary(value)), by = c("distance", "variable")]
  setnames(summ, 2, "sample")
  as.data.frame(summ)
}


