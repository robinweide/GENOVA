#' Insulation matrix plot.
#' 
#' Display a Hi-C matrix with insulation scores annotated in the margins.
#'
#' @param exp1 A GENOVA \code{contacts} object.
#' @param exp2 A GENOVA \code{contacts} object compatible with \code{exp1} 
#'   (optional) or \code{NULL}. If not \code{NULL}, the \code{exp2} argument 
#'   will be displayed alongside the \code{exp1} data.
#' @param IS_discovery A \code{IS_discovery} object, as created by the 
#'   \code{\link[GENOVA]{insulation_score}} function. When \code{NULL}, the 
#'   insulation score will be calculated on the fly for this region only.
#' @param chrom One of the following: \itemize{
#'     \item{A \code{character} of length 1 giving the chromosome name of the 
#'     region of interest.}
#'     \item{A 3-column, 1-row \code{data.frame} in BED-format.}
#'     \item{A single \code{character} describing a locus in UCSC notation, e.g.
#'     \code{"chr1:30,000,000-40,000,000"}.}
#'   } The latter two options automatically provide the \code{start} and 
#'   \code{end} arguments too.
#' @param start,end An \code{integer} with the start- and end-positions of the
#'   region of interest in basepairs. 
#' @param colour_lim A \code{numeric} of length two giving the minimum and 
#'   maximum values in contact values. Out of bounds values will be set to the
#'   nearest limits.
#' @param rasterise Set to \code{TRUE} to use a bitmap raster instead of 
#'   polygons. 
#' @param colour_bar A \code{logical} of length 1, indicating whether a
#'   colour-bar legend should be drawn at the right.
#'   
#' @details The 'on-the-fly' calculation of the insulation score uses the 
#'  default settings. Therefore, if for example the window size is to be 
#'  changed, please use the \code{insulation_score()} function and provide the
#'  result as the \code{IS_discovery} argument in this function.
#' 
#'  Because the 'on-the-fly' calculation of the insulation score only
#'  computes the insulation score for this locus, there might a small offset
#'  to the scores due to the limited scope for normalisation.
#'  
#' @section Recommended resolution: 10kb-40kb
#'
#' @return No values are returned but a plot is outputted on the graphics 
#'   device.
#' @export
#' 
#' @family matrix plots
#' @seealso The \code{\link[GENOVA]{insulation_score}} function for computation
#'   of insulation scores.
#'
#' @examples
#' \dontrun{
#' insulation_matrixplot(
#'   WT_20kb, KO_20kb,
#'   chrom ="chr1:40,000,000-60,000,000"
#' )
#' }
insulation_matrixplot <- function(
  exp1, exp2 = NULL, 
  IS_discovery = NULL,
  chrom, start, end,
  colour_lim = NULL, rasterise = FALSE,
  colour_bar = FALSE
) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  
  loc <- standardise_location(chrom, start, end, singular = TRUE)
  
  if (is_contacts(exp2)) {
    explist <- list(exp1, exp2)
  } else {
    explist <- exp1
  }
  res <- resolution(exp1)
  explist <- check_compat_exp(explist)
  expnames <- expnames(explist)
  expdat <- lapply(explist, select_subset, 
                   chrom = loc$chrom, start = loc$start, end = loc$end)
  if (length(expdat) == 1) {
    expdat <- expdat[[1]]
  } else {
    mat <- expdat[[1]]$z
    expdat <- expdat[[2]]
    mat[upper.tri(mat)] <- expdat$z[upper.tri(expdat$z)]
    expdat$z <- mat
  }
  
  if (is.null(IS_discovery)) {
    is_data <- lapply(explist, subset, chrom = loc$chrom,
                      start = loc$start - 31 * res,
                      end = loc$end + 31 * res)
    IS_discovery <- insulation_score(is_data)
  }
  
  if (!all(expnames %in% expnames(IS_discovery))) {
    stop("Sample names in the `IS_discovery` do not match contacts object(s).",
         call. = FALSE)
  }
  ins_dat <- as.data.table(IS_discovery$insula_score)
  ins_dat <- ins_dat[chrom == loc$chrom & start >= loc$start & end <= loc$end]
  
  loclim <- c(loc$start, loc$end)
  locbreaks <- scales::breaks_pretty()(loclim)
  inslim <- lapply(expnames, function(i){ins_dat[, eval(as.symbol(i))]})
  inslim <- range(do.call(c, inslim))
  inslim <- scales::expand_range(inslim, 0.1)
  insbreaks <- scales::breaks_width(0.5)(inslim)
  # Setup layout
  if (length(explist) == 2) {
    lay <- matrix(1:4, 2)
    if (colour_bar) {
      layout(cbind(lay, c(5,6)), widths = c(0.2, 1, 0.2), 
             heights = c(0.2, 1), respect = TRUE)
    } else {
      layout(lay, widths = c(0.2, 1), heights = c(0.2, 1), respect = TRUE)
    }
    par(mar = c(0,0,0,0))
    plot.new()
    par(mar = c(3, rep(0.2, 3)))
    plot(
      x = ins_dat[, eval(as.symbol(expnames[2]))],
      y = ins_dat$end,
      col = attr(exp2, "colour"),
      type = 'l',
      ylim = rev(loclim), xlim = rev(inslim), 
      xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n"
    )
    axis(1, insbreaks, insbreaks, cex.axis = 1)
  } else {
    lay <- matrix(1:2, 2)
    if (colour_bar) {
      layout(cbind(lay, 3:4), 
             widths = c(1, 0.2), heights = c(0.2, 1), respect = TRUE)
    } else {
      layout(lay, widths = 1, heights = c(0.2, 1), respect = TRUE)
    }
  }
  par(mar = c(rep(0.2, 3), 3))
  plot(
    x = ins_dat$end,
    y = ins_dat[, eval(as.symbol(expnames[1]))],
    col = attr(exp1, "colour"),
    type = 'l',
    xlim = loclim, ylim = inslim,
    xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n"
  )
  axis(2, insbreaks, insbreaks, cex.axis = 1, las = 1)
  par(mar = c(3, 0.2, 0.2, 3))
  
  colours <- .choose_palette()
  colours <- colorRampPalette(colours)
  
  if (is.null(colour_lim)) {
    colour_lim <- c(0, quantile(expdat$z[!as.logical(diag(ncol(expdat$z)))], 
                                0.99))
    message(paste0("Limits set to [", colour_lim[1], ",", 
                   format(colour_lim[2], digits = 2), "]"))
  }
  expdat$z[] <- pmax(expdat$z, colour_lim[1])
  expdat$z[] <- pmin(expdat$z, colour_lim[2])
  
  image.default(expdat, ylim = rev(loclim), col = colours(255),
                xaxt = "n", yaxt = "n", axes = FALSE,
                useRaster = literalTRUE(rasterise))
  box(lwd = 1)
  axis(1, locbreaks, locbreaks / 1e6, col = "transparent", col.ticks = "black")
  axis(4, locbreaks, locbreaks / 1e6, col = "transparent", col.ticks = "black",
       las = 1)
  
  if (length(explist) == 2) {
    txt <- seq(loclim[1], loclim[2], length.out = 20)
    txt <- txt[c(2, 19)]
    text(x = rev(txt)[1], y = txt[1], labels = expnames[1], adj = c(1, 1))
    text(x = rev(txt)[2], y = txt[2], labels = expnames[2], adj = c(0, 0))
  }
  
  if (colour_bar) {
    plot.new()
    par(plt = c(0, 0.25, 0.1, 0.95))
    m <- t(as.matrix(seq(colour_lim[1], colour_lim[2], length.out = 255)))
    image(1, m[1,], m, col = colours(255), xaxt = "n", yaxt = "n",
          xlab = "", ylab = "")
    axis(side = 4, lwd = 0, lwd.ticks = 1, lend = 1)
    mtext("Contacts", side = 4, line = 2.5)
  }
}




# polygonise <- function(x, y, t = F, n = 81, lim = NULL) {
#   
#   if (!t) {
#     l <- list(x = c(x[1], x, tail(x, 1)),
#               y = c(0, y, 0))
#   } else {
# 
#     l <- list(x = c(0, x, 0), 
#               y = c(y[1], y, tail(y, 1)))
#   }
#   
#   if (is.null(lim)) {
#     seq <- max(abs(l$y)) * c(-1, 1)
#   } else {
#     seq <- max(abs(lim)) * c(-1, 1)
#   }
#   seq <- seq(seq[1], seq[2], length.out = n)
#   cols <- bezier_corrected_divergent
#   cols <- rev(colorRampPalette(cols)(n-1))
#   
#   if (t) {
#     for (i in head(seq_along(seq), -1)) {
#       if (all(l$x < seq[i]) | all(l$x > seq[i + 1])) {
#         next()
#       }
#       xx <- .clip_poly(l$x, l$y, xrange = c(seq[i], seq[i + 1]))
#       polygon(xx$x, xx$y, col = cols[i], border = NA)
#     }
#   } else {
#     for (i in head(seq_along(seq), -1)) {
#       if (all(l$y < seq[i]) | all(l$y > seq[i + 1])) {
#         next()
#       }
#       xx <- .clip_poly(l$x, l$y, yrange = c(seq[i], seq[i + 1]))
#       polygon(xx$x, xx$y, col = cols[i], border = NA)
#     }
#   }
#   lines(x, y)
# }
