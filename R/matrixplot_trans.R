#' Trans matrix plot
#' 
#' This function plots a matrix of contacts data for the intersection of two
#' regions. Using two regions on different chromosomes plots trans data.
#'
#' @param exp1 A GENOVA \code{contacts} object.
#' @param exp2 A GENOVA \code{contacts} object compatible with \code{exp1} 
#'   (optional) or \code{NULL}. If not \code{NULL} the \code{exp2} argument will 
#'   be subtracted from the \code{exp1} data.
#' @param chrom_up,chrom_down One of the following: \itemize{
#'   \item{A \code{character} of length 1 giving the chromosome name of the 
#'   region of interest}
#'   \item{A 3-column, 1-row \code{data.frame} in BED-format}
#'   \item{A single \code{character} describing a locus in UCSC notation, e.g.
#'   \code{"chr1:30,000,000-40,000,000"}.}
#' } The latter two options automatically provide the \code{start_*} and 
#'   \code{end_*} arguments too.
#' @param start_up,start_down An \code{integer} with the start position of 
#'   the regions of interest.
#' @param end_up,end_down An \code{integer} with the end position of the regions
#'   of interest.
#' @param colour_lim A \code{numeric} of length two giving the minimum and
#'   maximum values in contacts. Out of bounds values will be set to the nearest
#'   limits.
#' @param rasterise Set to \code{TRUE} to use a bitmap raster instead of 
#'   polygons.
#'   
#' @section Resolution recommendation: A resolution in the ballpark of 
#'   \code{(end - start) / 500}.
#' @param colour_bar A \code{logical} of length 1, indicating whether a
#'   colour-bar legend should be drawn at the right.
#'
#' @return No values are returned but a plot is outputted to the graphics 
#'   device.
#' @export
#' 
#' @family matrix plots
#'
#' @examples
#' \dontrun{
#' trans_matrixplot(exp, 
#'                  chrom_up = "chr1:30,000,000-40,000,000",
#'                  chrom_down = "chr2:50,000,000-60,000,000")
#' }
trans_matrixplot <- function(
  exp1, exp2 = NULL,
  chrom_up, start_up, end_up,
  chrom_down, start_down, end_down,
  colour_lim = NULL, rasterise = FALSE,
  colour_bar = FALSE
) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  dat <- select_arbitrary(
    exp1,
    chrom_up, start_up, end_up,
    chrom_down, start_down, end_down
  )
  div_scale <- attr(exp1, "znorm")
  if (!is.null(exp2)) {
    if (attr(exp2, "znorm") != div_scale) {
      stop("Cannot combine z-normalised experiments with regular one.",
           call. = FALSE)
    }
    exp2 <- select_arbitrary(
      exp2,
      chrom_up, start_up, end_up,
      chrom_down, start_down, end_down
    )
    checks1 <- list(dim(dat$z), dat$x, dat$y, dat$chrom)
    checks2 <- list(dim(exp2$z), exp2$x, exp2$y, exp2$chrom)
    if (!identical(checks1, checks2)) {
      stop("`exp1` and `exp2` should have identical resolution and indices.",
           call. = FALSE)
    }
    dat$z <- dat$z - exp2$z
    div_scale <- TRUE
    if (is.null(colour_lim)) {
      colour_lim <- max(abs(dat$z)) * c(-1, 1)
    }
  }
  
  
  
  if (!is.null(colour_lim)) {
    if (length(colour_lim) == 2 && is.numeric(colour_lim)) {
      dat$z[] <- pmin(dat$z, colour_lim[2])
      dat$z[] <- pmax(dat$z, colour_lim[1])
    } else {
      stop("The `colour_lim` argument should be numeric and of length 2.")
    }
  } else {
    colour_lim <- range(dat$z)
  }
  
  chroms <- dat$chrom
  dat$chrom <- NULL
  
  colours <- if (div_scale) {
    bezier_corrected_divergent
  } else {
    .choose_palette()
  }
  colours <- colorRampPalette(colours)
  if (colour_bar) {
    lay <- matrix(c(1,2), 1, 2)
    layout(lay, widths = c(1, 0.2))
  }
  
  
  image(dat, col = colours(255), axes = FALSE, ylim = rev(range(dat$y)),
        useRaster = literalTRUE(rasterise), zlim = colour_lim,
        xlab = chroms[[1]], ylab = chroms[[2]])
  
  box(lwd = 2)
  
  region <- range(dat$x)
  if (diff(region) > 40e6) {
    region <- region / {mult <- 10e6}
  } else if (diff(region) > 2e6) {
    region <- region / {mult <- 1e6}
  } else {
    region <- region / {mult <- 1e5}
  }
  region <- c(floor(region[1]), ceiling(region[2]))
  region <- seq(region[1], region[2])
  axis(
    1, at = region * mult, 
    labels = region * ifelse(mult == 10e6, 10, 1), 
    lwd = 2, cex.axis = 1
  )
  
  region <- range(dat$y)
  if (diff(region) > 40e6) {
    region <- region / {mult <- 10e6}
  } else if (diff(region) > 2e6) {
    region <- region / {mult <- 1e6}
  } else {
    region <- region / {mult <- 1e5}
  }
  region <- c(floor(region[1]), ceiling(region[2]))
  region <- seq(region[1], region[2])
  axis(
    2, at = region * mult, 
    labels = region * ifelse(mult == 10e6, 10, 1), 
    lwd = 2, cex.axis = 1
  )
  
  if (colour_bar) {
    par(plt = c(0, 0.25, 0.15, 0.85))
    m <- t(as.matrix(seq(colour_lim[1], colour_lim[2], length.out = 255)))
    image(1, m[1,], m, col = colours(255), xaxt = "n", yaxt = "n",
          xlab = "", ylab = "")
    axis(side = 4, lwd = 0, lwd.ticks = 1, lend = 1)
    title <- "Contacts"
    if (div_scale) {
      if (is.null(exp2)) {
        title <- "Z-score"
      } else {
        title <- "Difference"
      }
    }
    mtext(title, side = 4, line = 2.5)
  }
}



select_arbitrary <- function(
  exp,
  chrom_up, start_up, end_up,
  chrom_down, start_down, end_down
) {
  loc <- rbind(
    standardise_location(chrom_up, start_up, end_up, singular = TRUE),
    standardise_location(chrom_down, start_down, end_down, singular = TRUE)
  )
  bins <- mapply(
    seq.int, by = 1,
    from = bed2idx(exp$IDX, loc, "start"),
    to = bed2idx(exp$IDX, loc, "end"),
    SIMPLIFY = FALSE
  )
  idx <- do.call(CJ, bins)
  idx2 <- idx[, list(V1 = pmin(V1, V2), V2 = pmax(V1, V2))]
  i <- idx[, cbind(V1 - min(V1) + 1, V2 - min(V2) + 1)]
  
  dat <- exp$MAT[idx2,]
  dat[is.na(V3), V3 := 0]
  
  mat <- matrix(0, nrow = length(bins[[1]]), ncol = length(bins[[2]]))
  mat[i] <- dat$V3
  list(
    x = exp$IDX[bins[1], (V2 + V3) / 2, on = "V4"],
    y = exp$IDX[bins[2], (V2 + V3) / 2, on = "V4"],
    z = mat,
    chrom = loc$chrom
  )
}



