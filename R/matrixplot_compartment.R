#' Compartment score matrix plot
#' 
#' Display a Hi-C matrix with compartment scores annotated in the margins.
#'
#' @inheritParams insulation_matrixplot
#' @param CS_discovery A \code{CS_discovery} object, as created by the 
#'   \code{\link[GENOVA]{compartment_score}} function. When \code{NULL}, the
#'   compartment score will be calculated on the fly for this arm only.
#' @param chrom A \code{character} of length 1 with a chromosome name.
#' @param arm Either \code{"p"} or \code{"q"} denoting the chromosome arm of 
#'   which to display the data from.
#' @param metric One of the follow: \describe{
#'   \item{\code{"contacts"}}{Displays contacts.}
#'   \item{\code{"obsexp"}}{Displays observed over expected by distance.}
#'   \item{\code{"correlation"}}{Displays Pearson correlation for the observed
#'     over expected by distance.}
#' }
#' @param ... Arguments passed to the \code{\link[GENOVA]{compartment_score}}
#'   function.
#' 
#' @details The 'on-the-fly' calculation of the compartment score uses the 
#'   default settings. Therefore, if for example a different eigenvector is to
#'   be retrieved, please use the \code{compartment_score} function and provide
#'   the result as the \code{CS_discovery} argument in this function.
#'   
#' @section Recommended resolution: 100kb-150kb.
#'
#' @return No values are returned but a plot is outputted on the graphics 
#'   device.
#' @export
#' 
#' @family matrix plots
#' @seealso The \code{\link[GENOVA]{compartment_score}} function for computation
#'   of compartment scores.
#'
#' @examples
#' \dontrun{ 
#' compartment_matrixplot(
#'   WT_150kb, KO_150kb,
#'   chrom = "chr18", arm = "q",
#'   metric = "correlation"
#' )
#' }
compartment_matrixplot <- function(
  exp1,
  exp2 = NULL,
  CS_discovery = NULL,
  chrom, arm = "p",
  colour_lim = NULL,
  rasterise = FALSE,
  metric = c("contacts", "obsexp", "correlation"),
  ...
) {
  metric <- match.arg(metric, c("contacts", "obsexp", "correlation"))
  if (!(chrom %in% exp1$CENTROMERES$chrom)) {
    stop("No centromere information found for this chromosome.",
         call. = FALSE)
  }
  
  # Finding arm definitions
  centro <- exp1$CENTROMERES$chrom == chrom
  centro <- exp1$CENTROMERES[centro]
  chr <- exp1$IDX[V1 == chrom]
  if (arm == "p") {
    chr <- chr[V4 < centro$start]
  } else {
    chr <- chr[V4 > centro$end]
  }
  if ((tail(chr$V3, 1) - head(chr$V2, 1) < 1e6) | nrow(chr) < 10) {
    stop("No `", arm, "` found. Try the other")
  }
  chr <- chr[, list(V1 = V1[1], V2 = V2[2], V3 = tail(V3, 1))]
  
  if (is_contacts(exp2)) {
    explist <- list(exp1, exp2)
  } else {
    explist <- exp1
  }
  explist <- check_compat_exp(explist)
  expnames <- expnames(explist)
  expdat <- lapply(explist, select_subset, 
                   chrom = chr$V1, start = chr$V2, end = chr$V3)
  
  if (is.null(CS_discovery)) {
    cs_data <- lapply(explist, subset, chrom = chr$V1,
                      start = chr$V2, end = chr$V3)
    cs_data[[1]]$CENTROMERES <- centro
    CS_discovery <- compartment_score(cs_data, ...)
  }
  
  if (!all(expnames %in% expnames(CS_discovery))) {
    stop("Sample names in the `IS_discovery` do not match contacts object(s).",
         call. = FALSE)
  }
  comp_dat <- as.data.table(CS_discovery$compart_scores)
  comp_dat <- comp_dat[chrom == chr$V1 & start >= chr$V2 & end <= chr$V3]
  
  # Change metrics
  div_scale <- FALSE
  if (metric != "contacts") {
    expdat <- lapply(expdat, function(xp) {
      mat <- xp$z
      id <- abs(row(mat) - col(mat)) + 1
      keep <- mat != 0
      means <- vapply(split(mat[keep], id[keep]), mean, numeric(1))
      mat[keep] <- mat[keep] / means[id[keep]] 
      mat[!keep] <- 1
      xp$z <- mat
      return(xp)
    })
    div_scale <- TRUE
  }
  if (metric == "correlation") {
    expdat <- lapply(expdat, function(xp) {
      xp$z <- suppressWarnings(cor(xp$z))
      return(xp)
    })
  }

  if (length(expdat) == 1) {
    expdat <- expdat[[1]]
  } else {
    mat <- expdat[[1]]$z
    expdat <- expdat[[2]]
    mat[upper.tri(mat)] <- expdat$z[upper.tri(expdat$z)]
    expdat$z <- mat
  }
  
  loclim <- c(chr$V2, chr$V3)
  locbreaks <- scales::breaks_pretty()(loclim)
  complim <- lapply(expnames, function(i){comp_dat[, eval(as.symbol(i))]})
  complim <- range(do.call(c, complim))
  complim <- scales::expand_range(complim, 0.1)
  compbreaks <- scales::breaks_width(1)(complim)
  
  if (length(explist) == 2) {
    lay <- matrix(1:4, 2)
    layout(lay, widths = c(0.2, 1), heights = c(0.2, 1), respect = TRUE)
    par(mar = c(0,0,0,0))
    plot.new()
    par(mar = c(3, rep(0.2, 3)))
    plot(
      x = comp_dat[, eval(as.symbol(expnames[2]))],
      y = comp_dat$end,
      col = attr(exp2, "colour"),
      type = 'l',
      ylim = rev(loclim), xlim = rev(complim), 
      xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n"
    )
    axis(1, compbreaks, compbreaks, cex.axis = 1)
  } else {
    lay <- matrix(1:2, 2)
    layout(lay, widths = 1, heights = c(0.2, 1), respect = TRUE)
  }
  
  par(mar = c(rep(0.2, 3), 3))
  plot(
    x = comp_dat$end,
    y = comp_dat[, eval(as.symbol(expnames[1]))],
    col = attr(exp1, "colour"),
    type = 'l',
    xlim = loclim, ylim = complim,
    xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n"
  )
  axis(2, compbreaks, compbreaks, cex.axis = 1, las = 1)
  par(mar = c(3, 0.2, 0.2, 3))
  
  if (div_scale) {
    colours <- bezier_corrected_divergent
  } else {
    colours <- .choose_palette()
  }
  colours <- colorRampPalette(colours)
  
  if (is.null(colour_lim)) {
    nondiag <- expdat$z[!as.logical(diag(ncol(expdat$z)))]
    colour_lim <- switch(
      metric,
      contacts = c(0, quantile(nondiag, 0.99, na.rm = TRUE)),
      obsexp = {
        q <- quantile(nondiag, c(0.01, 0.99), na.rm = TRUE)
        (max(abs(q - 1)) * c(-1, 1)) + 1
      },
      correlation = c(-1, 1)
    )
    message(paste0("Limits set to [", format(colour_lim[1], digits = 2), ",", 
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

}



