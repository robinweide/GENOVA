
# Documentation -----------------------------------------------------------

#' @title Quantification of results
#' @name quantify
#' 
#' @description A good amount of \code{discovery} objects can be quantified.
#' What exactly is to be quantified differs per \code{discovery} type.
#' \describe{
#'   \item{saddle disoveries}{can be used to compute compartment strengths.}
#'   \item{ARMLA discoveries}{such as APA, PE-SCAn, ATA and ARA compare
#'   different regions of their outputs.}
#'   \item{IIT discoveries}{summarise their values by neighbours.}
#' }
#' 
#' 
#' @param discovery A \code{discovery} object as returned by GENOVA analysis functions.
#' @param ... further arguments passed to or from other methods.
#' take the middle 3x3 matrix of the APA).
#' 
#' @section Shapes:
#' 
#' The quantification of ARMLA discoveries require a shape to distinguish 
#' regions to quantify.
#' 
#' \subsection{ARA, PESCAn and CSCAn}{
#' APA, PESCAn and CSCAn require one of the following:
#' \itemize{
#'   \item{\code{"center_vs_quadrants"}}
#'   \item{\code{"center_vs_rest"}}
#'   \item{\code{"circle"}}
#' }
#' The \code{size} parameter determines the number of bins of the central 
#' foreground.
#' 
#' In the illustrations below, red is considered foreground and blue is 
#' considered background.
#' 
#' The \code{"center_vs_quadrants"} option does not include region directly
#' horizontal or vertical of the centre as background.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-center-vs-quandrants.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-center-vs-quandrants.png}\out{\end{center}}
#'   }
#'
#' The \code{"center_vs_rest"} option sees everything but the centre as background.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-center-vs-rest.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-center-vs-rest.png}\out{\end{center}}
#'   }
#' 
#' The \code{"circle"} option is like \code{"center_vs_rest"} but rounds corners 
#' of the central foreground. Note that for \code{size <= 3} these two options 
#' are equivalent.
#' 
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-circle.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-circle.png}\out{\end{center}}
#'   }
#' }
#' 
#' \subsection{ATA}{
#' ATA requires one of the following:
#' \itemize{
#'   \item{\code{"insulation"}}
#'   \item{\code{"cornerpeak"}}
#'   \item{\code{"checker"}}
#' }
#' 
#' In the illustrations below, red is considered foreground and blue is 
#' considered background. The line indicates the diagonal.
#' 
#' The \code{"insulation"} option compares within-TAD contacts to between-TAD
#' contacts.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-insulation.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-insulation.png}\out{\end{center}}
#'   }
#' 
#' The \code{"cornerpeak"} option compares the intersection of boundaries versus
#' within-TAD contacts.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-cornerpeak.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-cornerpeak.png}\out{\end{center}}
#'   }
#'   
#' The \code{"checker"} option compares the within-TAD contacts to the 
#' contacts between it's immediate neighbours.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-checker.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-checker.png}\out{\end{center}}
#'   }
#' }
#' 
#' \subsection{ARA}{
#' ARA requires one of the following:
#' \itemize{
#'   \item{\code{"ARA"}}
#'   \item{\code{"stripes"}}
#' }
#' The \code{size} argument controls the width of the stripes.
#' 
#' The \code{"ARA"} option reports about both 3' and 5' stripes and regions as
#' well as the bins that span the locus. They are indicated in different colours 
#' below.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-ARA.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-ARA.png}\out{\end{center}}
#'   }
#' 
#' The \code{"stripes"} options reports the values and distances of the stripes.
#' The 5' distances are encoded as negative, whereas 3' distances are positive.
#' 
#' \if{html}{
#'   \out{<div style="text-align: center">}\figure{quant-stripes.png}{options: style="width:62px;max-width:100\%;"}\out{</div>}
#'   }
#' \if{latex}{
#'   \out{\begin{center}}\figure{quant-stripes.png}\out{\end{center}}
#'   }
#' 
#' }
#'
#' @export
#' @examples 
#' NULL
quantify <- function(discovery, ...) {
  UseMethod("quantify", discovery)
}

# Functions ---------------------------------------------------------------

# If you fail it is best to fail graciously

#' @export
#' @rdname quantify
#' @usage NULL
quantify.default <- function(discovery, ...) {
  stop("No quantify method for class '", class(discovery),
       "' has been implemented.", call. = FALSE)
}


#' Standard quantification of ARMLA output.
#' 
#' This is a common interface for standard ARMLA quantifications. It is for
#' internal use only. Specific methods should extract the relevant fields and
#' feed it into this function. This is just documented so people who want to
#' read to code have some context.
#'
#' @param aggregate The summarised output in a discovery object.
#' @param raw The raw output in a discovery object.
#' @param expnames A character with experiment names
#' @param shape A function that when given the dimensions of 'aggregate' should
#'   produce a list with 'foreground' and 'background' matrices matching those 
#'   dimensions. See the 'shape_*' functions down this document.
#' @param IDX The IDX of the original experiment. Just used when the rownames
#'   of the raw data indicate bins instead of loci. This happens in for example
#'   the extended loops APA.
#' @param fun A summarising function that takes a numeric vector as input and
#'   outputs a length 1 summary. Typically "mean" or "median".
#' @param ... Just here for S3 method consistency and catch bad arguments. 
#'   Doesn't do anything.
#'
#' @return A list with a per sample quantification and a per feature 
#' quantification.
#' 
#' @noRd
#' @examples \dontrun{
#' # See quantify methods for APA/PESCAn/ATA for usage.
#' }
quantify.ARMLA <- function(
  aggregate, 
  raw = NULL, 
  expnames = NULL,
  shape = shape_center_vs_quadrants(3),
  IDX = NULL,
  fun = median,
  ...
) {
  if (!(is.array(aggregate) | is.matrix(aggregate))) {
    stop("The aggregate should be an array.", call. = FALSE)
  }
  dim <- dim(aggregate)
  if (length(expnames) == 0) {
    expnames <- paste0("exp", seq_len(tail(dim, 1)))
  }
  
  if (!is.function(shape)) {
    stop("The shape argument should be a function.", call. = FALSE)
  }
  
  shape <- shape(dim)
  foreground <- shape$foreground
  background <- shape$background
  
  # Split signal into samples
  samples <- split(aggregate, slice.index(aggregate, 3))
  # Calculate medians per sample
  metrics  <- vapply(samples, function(x) {
    c(fun(x[foreground]), fun(x[background]))
  }, numeric(2))
  # Format global sample stats
  global <- data.frame(
    sample = dimnames(aggregate)[[3]],
    foreground = metrics[1, ],
    background = metrics[2, ],
    foldchange = metrics[1, ] / metrics[2, ],
    difference = metrics[1, ] - metrics[2, ]
  )
  
  if (is.null(raw)) {
    warning("No raw data found. Returning global summaries.")
    return(list(per_sample = global))
  }
  
  # Loop over experiments
  local <- lapply(seq_along(raw), function(i) {
    arr <- raw[[i]]
    arr[is.na(arr)] <- 0
    # Split sample into individual loops
    loops <- split(arr, slice.index(arr, 1))
    # Calculate medians per loops
    metrics <- vapply(loops, function(x) {
      c(fun(x[foreground], na.rm = TRUE), fun(x[background], na.rm = TRUE))
    }, numeric(2))
    # Format loop stats
    loc <- dimnames(arr)[[1]]
    loc <- if (is.null(loc)) seq_len(nrow(arr)) else loc
    data.table(
      sample = i,
      loc = loc,
      foreground = metrics[1, ],
      background = metrics[2, ],
      foldchange = metrics[1, ] / metrics[2, ],
      difference = metrics[1, ] - metrics[2, ]
    )
  })
  local <- rbindlist(local)
  local$sample <- expnames[local$sample]
  
  # Format location
  loc <- interpret_location_string(local$loc, IDX)
  local$loc <- NULL
  local <- cbind(local, loc)
  local <- as.data.frame(local)
  
  return(list(per_sample = global, per_entry = local))
}

#' @rdname quantify
#' @param size An \code{integer} of length one to determine the size of features
#'   of interest in bins.
#' @param metric Either \code{"median"} or \code{"mean"} to summarise features.
#' @param shape A character of length 1 specifying what shape to use. See the 
#'   section shapes.
#' @param IDX The \code{IDX} part of a contacts object. Used only in converting 
#'   features expressed in bins back to genomic space. This is rarely needed,
#'   but is useful for APAs ran with extended loops where features aren't 1:1
#'   traceable to the input.
#' @export
quantify.APA_discovery <- function(
  discovery, size = 3, 
  metric = "median",
  shape = "center_vs_quadrants", IDX = NULL,
  ...
) {
  
  shape <- parse_shape_arg(
    shape, size, 
    c("center_vs_quadrants", "center_vs_rest", "circle")
  )
  
  metric <- match.arg(metric, c("mean", "median"))
  metric <- switch(
    metric,
    mean = mean.default,
    median = median.default
  )
  
  out <- quantify.ARMLA(
    aggregate = discovery$signal, 
    raw = discovery$signal_raw, 
    expnames = expnames(discovery),
    shape = shape,
    IDX = IDX,
    fun = metric,
    ...
  )
  
  if (length(out) == 2) {
    names(out) <- c(names(out)[1], "per_loop")
  }
  
  return(out)
}

#' @rdname quantify
#' @export
quantify.PESCAn_discovery <- function(
  discovery, size = 5, 
  metric = "median",
  shape = "circle", IDX = NULL,
  ...
) {
  
  metric <- match.arg(metric, c("mean", "median"))
  metric <- switch(
    metric,
    mean = mean.default,
    median = median.default
  )
  
  shape <- parse_shape_arg(
    shape, size,
    c("center_vs_quadrants", "center_vs_rest", "circle")
  )
  
  # Normalize row values for grand median background
  raw <- lapply(seq_along(discovery$signal_raw), function(i) {
    bg <- median(discovery$shifted[, , i])
    raw <- discovery$signal_raw[[i]]
    raw[is.na(raw)] <- 0
    raw[] <- raw / bg
    raw
  })
  
  out <- quantify.ARMLA(
    aggregate = discovery$obsexp, 
    raw = raw, 
    expnames = expnames(discovery),
    shape = shape,
    IDX = IDX,
    fun = metric,
    ...
  )
  
  if (length(out) == 2) {
    names(out) <- c(names(out)[1], "per_interaction")
  }
  out <- lapply(out, function(x) {
    x$foldchange <- NULL # Fold change from obs/exp over obs/exp is meaningless
    x
  })
  
  return(out)
}

#' @rdname quantify
#' @export
quantify.CSCAn_discovery <- function(
  discovery, size = 5,
  metric = "median",
  shape = "circle", IDX = NULL,
  ...
) {
  metric <- match.arg(metric, c("mean", "median"))
  metric <- switch(
    metric,
    mean = mean.default,
    median = median.default
  )
  
  shape <- parse_shape_arg(
    shape, size,
    c("center_vs_quadrants", "center_vs_rest", "circle")
  )
  
  groups <- lapply(discovery$signal_raw, attr, "group")
  
  # Normalize row values for grand median background
  raw <- lapply(seq_along(discovery$signal_raw), function(i) {
    bg <- median(discovery$shifted[, , , i])
    raw <- discovery$signal_raw[[i]]
    raw[is.na(raw)] <- 0
    raw[] <- raw / bg
    raw
  })
  
  agg <- discovery$obsexp
  dnames <- dimnames(agg)
  collapse <- CJ(dnames[[3]], dnames[[4]])
  dnames[[3]] <- do.call(paste, 
                         c(unname(as.list(collapse)),
                           sep = "~&~"))
  dim <- dim(agg)
  dim(agg) <- c(dim[1], dim[2], prod(dim[3:4]))
  dimnames(agg) <- dnames[1:3]
  
  out <- quantify.ARMLA(
    aggregate = agg,
    raw = raw,
    expnames = expnames(discovery),
    shape = shape,
    IDX = IDX,
    fun = metric,
    ...
  )
  
  if (length(out) == 2) {
    names(out) <- c(names(out)[1], "per_interaction")
  }
  out <- lapply(out, function(x) {
    x$foldchange <- NULL # Fold change for obsexp is pointless
    x
  })
  sample <- out$per_sample$sample
  sample <- tstrsplit(sample, "~&~")
  out$per_sample$sample <- sample[[2]]
  out$per_sample$group <- sample[[1]]
  
  groups <- attr(raw[[1]], "group")
  out$per_interaction$group <- groups[out$per_interaction$feature_id]
  
  return(out)
}

#' @rdname quantify
#' @export
quantify.ATA_discovery <- function(
  discovery, size = 3, 
  metric = "median",
  shape = "insulation", 
  IDX = NULL,
  ...
) {
  
  metric <- match.arg(metric, c("mean", "median"))
  metric <- switch(
    metric,
    mean = mean.default,
    median = median.default
  )
  
  padding <- attr(discovery, "padding")
  padding <- if (is.null(padding)) 1 else padding
  shape <- parse_shape_arg(
    shape, padding,
    c("insulation", "cornerpeak", "checker")
  )
  
  out <- quantify.ARMLA(
    aggregate = discovery$signal, 
    raw = discovery$signal_raw, 
    expnames = expnames(discovery),
    shape = shape,
    IDX = IDX,
    fun = metric,
    ...
  )
  
  if (length(out) == 2) {
    names(out) <- c(names(out)[1], "per_TAD")
  }
  
  return(out)
}

#' @rdname quantify
#' @export
quantify.ARA_discovery <- function(discovery, size = 3, shape = "ARA", ...) {
  aggregate <- discovery$obsexp
  
  dim <- dim(aggregate)
  
  expnames <- expnames(discovery)
  if (length(expnames) == 0) {
    expnames <- paste0("exp", seq_len(tail(dim, 1)))
  }
  
  shape <- match.arg(shape, c("ARA", "stripes"))
  shape <- switch(shape,
                  ARA = shape_ARA(size),
                  stripes = shape_stripes(size, resolution(discovery))
  )
  shape <- shape(dim)
  
  ncat <- unique(as.vector(shape))
  ncat <- length(ncat[!is.na(ncat)])
  
  samples <- split(aggregate, slice.index(aggregate, 3))
  
  metrics <- vapply(samples, function(x) {
    x <- vapply(split(x, shape), mean, numeric(1))
  }, numeric(ncat))
  
  global <- data.frame(
    sample = expnames[as.vector(col(metrics))],
    feature = as(rownames(metrics)[as.vector(row(metrics))], 
                 typeof(shape)),
    value = as.vector(metrics)
  )
  global
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


# Retired
quantify_old.APA_discovery <- function(discovery, signal_size = 3, ...) {
  
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


# Shapes ------------------------------------------------------------------

parse_shape_arg <- function(shape, size, valid) {
  if (is.function(shape)) {
    return(shape)
  }
  shape <- match.arg(shape, valid)
  shape <- switch(
    shape,
    center_vs_quadrants = shape_center_vs_quadrants(size),
    center_vs_rest = shape_center_vs_rest(size),
    circle = shape_circle(size),
    insulation = shape_TAD_insulation(size),
    cornerpeak = shape_TAD_cornerpeak(size),
    checker = shape_TAD_checker(size),
    stop("Invalid shape argument.", call. = FALSE)
  )
  return(shape)
}



shape_center_vs_quadrants <- function(size) {
  force(size)
  function(dim) {
    # Select central rows / cols
    mid <- floor(dim[1] / 2) + 1 + c(-1, 1) * (size - 1)/2
    mid <- seq.int(mid[1], mid[2], by = 1L)
    
    # Set foreground/background
    foreground <- background <- matrix(FALSE, dim[1], dim[2])
    foreground[mid, mid] <- TRUE
    background[-mid, -mid] <- TRUE
    
    return(list(foreground = foreground, background = background))
  }
}

shape_center_vs_rest <- function(size) {
  force(size)
  function(dim) {
    # Select central rows / cols
    mid <- floor(dim[1] / 2) + 1 + c(-1, 1) * (size - 1)/2
    mid <- seq.int(mid[1], mid[2], by = 1L)
    
    # Set foreground/background
    foreground <- background <- matrix(FALSE, dim[1], dim[2])
    foreground[mid, mid] <- TRUE
    background <- !foreground
    
    return(list(foreground = foreground, background = background))
  }
}

shape_circle <- function(size) {
  force(size)
  function(dim) {
    size <- size / 2
    mat <- matrix(FALSE, dim[1], dim[2])
    mid <- (dim + 1)/2
    dist <- sqrt((row(mat) - mid[1])^2 + (col(mat) - mid[2])^2)
    
    foreground <- dist < size
    
    return(list(foreground = foreground, background = !foreground))
  }
}

shape_TAD_insulation <- function(padding) {
  force(padding)
  function(dim) {
    borders <- cumsum(c(padding - 0.5, 1)) / (2 * padding) * dim[1]
    borders <- round(borders)
    
    foreground <- background <- matrix(FALSE, dim[1], dim[2])
    fg <- seq(borders[1] + 1, borders[2] - 1, by = 1)
    foreground[fg, fg] <- TRUE
    
    bg <- c(seq(1, borders[1] - 1, by = 1), seq(borders[2] + 1, dim[1], by = 1))
    background[bg, fg] <- TRUE
    background[fg, bg] <- TRUE
    
    return(list(foreground = foreground, background = background))
  }
}

shape_TAD_cornerpeak <- function(padding) {
  force(padding)
  function(dim) {
    borders <- cumsum(c(padding - 0.5, 1)) / (2 * padding) * dim[1]
    borders <- round(borders)
    
    corner <- expand.grid(borders[1] + c(-2:2), borders[2] + c(-2:2))
    corner <- as.matrix(corner)
    
    foreground <- background <- matrix(FALSE, dim[1], dim[2])
    bg <- seq(borders[1] + 2, borders[2] - 2, by = 1)
    background[bg, bg] <- TRUE
    
    foreground[corner] <- TRUE
    foreground[cbind(corner[, 2], corner[, 1])] <- TRUE
    
    return(list(foreground = foreground, background = background))
  }
}

shape_TAD_checker <- function(padding) {
  force(padding)
  function(dim) {
    borders <- cumsum(c(padding - 0.5, 1)) / (2 * padding) * dim[1]
    borders <- round(borders)
    
    foreground <- background <- matrix(FALSE, dim[1], dim[2])
    fg <- seq(borders[1] + 1, borders[2] - 1, by = 1)
    foreground[fg, fg] <- TRUE
    
    begin <- seq(1, borders[1] - 1, by = 1)
    end <- seq(borders[2] + 1, dim[1], by = 1)
    
    background[end, begin] <- TRUE
    background[begin, end] <- TRUE
    
    return(list(foreground = foreground, background = background))
  }
}

# Only works in ARA
shape_ARA <- function(size) {
  force(size)
  function(dim) {
    template <- matrix(NA_character_, dim[1], dim[2])
    mid <- floor(dim[1] / 2) + 1 + c(-1, 1) * (size - 1)/2
    mid <- seq.int(mid[1], mid[2], by = 1L)
    
    up <- seq(1, mid[1] - 1, by = 1)
    down <- seq(tail(mid, 1) + 1, dim[1], by = 1)
    
    template[up, up] <- "5' Region"
    template[up, mid] <- "5' Stripe"
    template[up, down] <- "Span"
    template[mid, down] <- "3' Stripe"
    template[down, down] <- "3' Region"
    template[row(template) >= col(template)] <- NA
    return(template)
  }
}

# Only works in ARA
shape_stripes <- function(size, resolution) {
  force(size)
  force(resolution)
  function(dim) {
    template <- matrix(NA_integer_, dim[1], dim[2])
    midpoint <- floor(dim[1] / 2) + 1
    mid <- midpoint + c(-1, 1) * (size - 1)/2
    mid <- seq.int(mid[1], mid[2], by = 1L)
    
    up <- seq(1, mid[1] - 1, by = 1)
    down <- seq(tail(mid, 1) + 1, dim[1], by = 1)
    
    template[down, mid] <- down - midpoint
    template <- t(template)
    template[up, mid] <- up - midpoint
    return(template * resolution) 
  }
}