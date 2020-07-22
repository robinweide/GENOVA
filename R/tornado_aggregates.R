#' Tornado plot for aggregated discoveries
#' 
#' Takes the bins parallel to the Hi-C diagonal through the anchor locations
#' and perform k-means clustering on these.
#'
#' @param discovery One of the following \code{discovery} objects: \itemize{
#'   \item{An \code{APA_discovery} object}
#'   \item{A \code{CSCAn_discovery} object}
#'   \item{A \code{PESCAn_discovery} object}
#' }
#' @param znorm A \code{logical} of length one deciding whether to z-score
#'   normalize every feature within samples (\code{TRUE}) or not (\code{FALSE}).
#' @param K An \code{integer} of length or \code{NULL}, setting the number of
#'   clusters that \code{kmeans()} is instructed to find. When set to 
#'   \code{NULL}, a search for an appropriate \code{K} is conducted to explain
#'   a fraction of the variance controlled by the \code{min_var_exp} argument.
#' @param sort An \code{integer} setting sample indices on which the 
#'   row-ordering must be computed. When \code{NULL}, orders on the mean of all 
#'   samples.
#' @param colour_lim A \code{numeric} of length two setting the minimum and
#'   maximum values of the colour scale respectively.
#' @param min_var_exp A \code{numeric} between 0-1 that controls how much 
#'   variance the \code{kmeans()} result should explain before choosing a K. 
#'   When the \code{K} argument is not \code{NULL}, this is ignored.
#' @param print A \code{logical} whether to immediately render the plot.
#'
#' @return A \code{list} with two elements; first a plot and second a 
#'   \code{data.frame} with cluster annotations per feature.
#' @export
#'
#' @examples \dontrun{
#' apa <- APA(explist, bed)
#' tornado_aggregate(apa)
#' }
tornado_aggregate <- function(
  discovery, znorm = TRUE, K = NULL,
  sort = NULL, colour_lim = NULL, min_var_exp = 0.7,
  print = TRUE
) {
  # Input checks
  if (!("signal_raw" %in% names(discovery))) {
    stop("Please use a discovery with raw data.",
         call. = FALSE)
  }
  if (!inherits(discovery, 
                c("APA_discovery", "CSCAn_discovery", "PESCAn_discovery"))) {
    stop("The discovery must be one of the following: APA, CSCAn or PESCAn.",
         call. = FALSE)
  }
  
  # Grab metadata
  raw <- discovery$signal_raw
  dim <- dim(discovery$signal)
  expnames <- expnames(discovery)
  res <- resolution(discovery)
  
  # Seek features that are in all subsets
  shared_features <- lapply(raw, rownames)
  shared_features <- Reduce(intersect, shared_features)
  
  # Get (normalized) diagonal slices
  slices <- lapply(raw, function(x) {
    x <- x[shared_features, , , drop = FALSE]
    d <- dim(x)
    x <- x[slice.index(x, 2) == slice.index(x, 3)]
    x[is.na(x)] <- 0
    dim(x) <- d[1:2]
    if (znorm) {
      x[] <- t(apply(x, 1, scale))
    }
    x
  })
  
  # Combine slices
  d <- dim(slices[[1]])
  slices <- do.call(c, slices)
  dim(slices) <- c(d, length(raw))
  
  # Take midpoint values
  mid <- (ncol(slices) + 1) / 2
  midval <- slices[, mid, , drop = FALSE]
  dim(midval) <- dim(midval)[c(1, 3)]
  midval[!is.finite(midval)] <- 0
  
  KM <- NULL
  if (is.null(K)) {
    # Loop to find decent K for kmeans
    sumvar <- sum(apply(midval, 2, var, na.rm = TRUE)) * (d[1] - 1)
    K <- varexp <- 0
    while (varexp < min_var_exp) {
      K <- K + 1
      KM <- kmeans(midval, centers = K)
      varexp <- KM$betweenss / sumvar
    }
    K <- max(K, 1)
  }
  
  # Perform KM if K was provided
  if (is.null(KM)) {
    KM <- kmeans(midval, centers = K)
  }
  clust <- KM$cluster
  
  # Sort on cluster / mid value
  if (length(sort) == 0) {
    sort <- seq_len(ncol(midval))
  }
  sort <- rowMeans(midval[, sort, drop = FALSE])
  # Seek grand order of things
  ord <- order(clust, sort)
  # Seek within cluster order of things
  ord2 <- do.call(c, tapply(sort[ord], clust[ord], order))
  slices <- slices[ord, , , drop = FALSE]
  
  # Melt slices
  df <- data.frame(
    row = as.vector(slice.index(slices, 1)),
    col = (as.vector(slice.index(slices, 2)) - mid) * res,
    sample = expnames[as.vector(slice.index(slices, 3))],
    value = as.vector(slices)
  )
  df$row <- ord2[df$row]
  df$K <- clust[ord]
  
  # Setup plotting stuff
  if (is.null(colour_lim)) {
    colour_lim <- c(NA, NA)
  }
  sscale <- if (znorm) {
    scale_fill_GENOVA_div(name = "Z-Score", midpoint = 0, limits = colour_lim,
                          oob = scales::squish)
  } else {
    scale_fill_GENOVA(name = "Contacts", limits = colour_lim,
                      oob = scales::squish)
  }
  breaks <- scales::extended_breaks()(range(df$row))
  
  # Build plot
  g <- ggplot2::ggplot(df, ggplot2::aes(col, row, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_x_continuous(expand = c(0,0),
                                labels = scales::number_format(scale = 1e-3),
                                name = "Distance to anchor (kb)") +
    ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, NA),
                                name = "", breaks = breaks) +
    sscale +
    ggplot2::facet_grid(K ~ sample, space = "free_y", scales = "free_y") +
    GENOVA_THEME() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA, colour = "black"),
      panel.spacing.x = grid::unit(0.6, "strwidth", 
                                   data = as.character(max(df$col)/1000)),
      panel.spacing.y = grid::unit(0.9, "strheight", data = "1000")
    )
  
  features <- interpret_location_string(shared_features, feature_id = FALSE)
  features$cluster <- clust
  if (print) {
    print(g)
  }
  return(invisible(list(plot = g, feature_cluster = features)))
}