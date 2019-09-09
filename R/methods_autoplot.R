#' Visualise APA and PE-SCAn results
#'
#' Plot APA and PE-SCAn results, by default with comparisons of differences.
#'
#' @name RMT_visualisation
#' @aliases APA_visualisation PESCAn_visualisation
#' @param results A results object as returned by the \code{APA} function or the
#'   \code{PESCAn} function.
#' @param subtract An \code{integer} or \code{character} matching an experiment
#'   name of length 1, specifying which sample should be subtracted to make the
#'   difference panels. Alternatively, set to \code{NULL} for no difference
#'   panels.
#' @param mode (PE-SCAn only) A \code{character} vector of length 1 with
#'   \code{"obsexp"} for plotting the observed over expected or \code{"signal"}
#'   for plotting the signal only (no background correction).
#' @param raw A \code{logical} of length 1: should a bare mimimum plot be
#'   returned?
#'
#' @details Difference (\u0394) panels are created by subtracting the values of each
#'   sample by the values of the sample indicated by the \code{subtract}
#'   argument. Hence, a clear positive signal in the difference panels indicates
#'   enrichment in that sample versus the \code{subtract}-sample.
#'
#'   The \code{raw = TRUE} argument allows customisation of the plot: it returns
#'   a plot object with no position or fill scales and no theme settings.
#'
#'   When \code{raw = TRUE} and \code{subtract} is not \code{NULL}, the plot
#'   will not render unless a continuous fill scale is supplied with the
#'   argument \code{aesthetics = "fill2"}. This a small trade-off for having two
#'   independent fill scales in a single plot.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link[GENOVA]{APA}} and \code{\link[GENOVA]{PESCAn}}
#'
#' @examples
#' # APA
#' apa <- APA(list(WT = WT_40kb, KO = KO_40kb), WT_loops)
#' autoplot(apa)
#'
#' # PE-SCAn
#' pescan <- PESCAn(list(WT = WT_40kb, KO = KO_40kb), super_enhancers)
#' autoplot(pescan)
#'
#' # To plot PE-SCAn without background correction
#' autoplot(pescan, mode = "signal")
#'
#' # Handling 'raw' plots
#' autoplot(pescan, raw = TRUE) +
#'   ggplot2::scale_fill_gradient(aesthetics = "fill2")
#' @export
#' @rdname APA_PESCAn_visualisation
autoplot.APA_results <- function(results, subtract = 1, raw = FALSE) {
  # Name by order if names are missing
  if (is.null(names(results))) {
    names(results) <- seq_along(results)
  }
  # Melt matrices to data.frames
  mats <- lapply(results, `[[`, "signal")
  melts <- lapply(names(mats), function(i) {
    cbind(reshape2::melt(mats[[i]]), mode = "Individual", sample = i)
  })
  melts <- do.call(rbind, melts)

  # Calculate diff data if subtract is supplied
  if (!is.null(subtract)) {
    ref <- mats[[subtract]]
    diffs <- lapply(names(mats), function(i) {
      cbind(reshape2::melt(mats[[i]] - ref), mode = "\u0394", sample = i)
    })
    diffs <- do.call(rbind, diffs)
    melts <- rbind(melts, diffs)
  }

  # Setup base of plot
  g <- ggplot2::ggplot(melts, aes(Var1, Var2))

  if (!is.null(subtract)) {
    # Setup basics of the diff plots
    # Warnings are supressed because ggplot doesn't recognise
    # the fill2 aesthetic (YET!)
    suppressWarnings(
      g <- g + ggplot2::geom_raster(
        data = function(x) {
          x[x$mode == "\u0394", ]
        },
        aes(fill2 = value)
      ) +
        ggplot2::facet_grid(mode ~ sample, switch = "y")
    )

    if (!raw) {
      g <- g + ggplot2::scale_fill_gradientn(
        colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
        aesthetics = "fill2",
        name = "\u0394",
        limits = c(-1, 1) * max(abs(melts$value[melts$mode == "\u0394"]))
      )
      # Hack for scale
      g$scales$scales[[1]]$guide <- ggplot2::guide_colourbar()
      g$scales$scales[[1]]$guide$available_aes[[3]] <- "fill2"
      g$scales$scales[[1]]$guide$order <- 2
    } else {
      g <- g + guides("fill2" = guide_colourbar(available_aes = "fill2"))
    }

    # Good old ggnomics-style hack for different fill scales
    old_geom <- g$layers[[1]]$geom
    old_nahandle <- old_geom$handle_na
    new_nahandle <- function(self, data, params) {
      colnames(data)[colnames(data) %in% "fill2"] <- "fill"
      old_nahandle(data, params)
    }
    new_geom <- ggplot2::ggproto(paste0(sample(1e6, 1), class(old_geom)),
      old_geom,
      handle_na = new_nahandle
    )
    names(new_geom$default_aes)[1] <- "fill2"
    new_geom$non_missing_aes <- "fill2"
    g$layers[[1]]$geom <- new_geom
  } else {
    # No special treatment needed for just one colourscale
    g <- g + ggplot2::facet_grid(~sample)
  }

  # Add the non-diff plots
  g <- g + ggplot2::geom_raster(
    data = function(x) {
      x[x$mode == "Individual", ]
    },
    aes(fill = value)
  )

  if (raw) {
    return(g)
  }

  # Decorating of the plot
  g <- g + ggplot2::scale_fill_gradientn(
    colours = c("white", "orange", "red", "black"),
    guide = guide_colourbar(order = 1),
    name = "\u03bc Contacts"
  ) +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = function(x) {
        x <- scales::extended_breaks()(x)
        if (length(x) > 3) {
          head(tail(x, -1), -1)
        } else {
          x
        }
      },
      labels = function(x) {
        ifelse(x == 0, "3'", paste0(x / 1000, "kb"))
      }
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = function(x) {
        x <- scales::extended_breaks()(x)
        if (length(x) > 3) {
          head(tail(x, -1), -1)
        } else {
          x
        }
      },
      labels = function(x) {
        ifelse(x == 0, "5'", paste0(x / 1000, "kb"))
      }
    ) +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.placement = "outside",
      strip.background = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey30")
    )
  g
}

#' @export
#' @rdname APA_PESCAn_visualisation
autoplot.PESCAn_results <- function(results,
                                    subtract = 1,
                                    mode = c("obsexp", "signal"),
                                    raw = FALSE) {
  # Handle mode settings
  mode <- match.arg(mode)
  hasobsexp <- all(sapply(results, function(res) "obsexp" %in% names(res)))
  hasobsexp <- hasobsexp && mode == "obsexp"
  res <- lapply(results, function(res) {
    if (hasobsexp) {
      list(signal = res$obsexp)
    } else {
      list(signal = res$signal)
    }
  })

  # Cannibalise the APA autoplot function
  g <- autoplot.APA_results(res, subtract = subtract, raw = raw)

  # Adjust for PESCAn
  if (hasobsexp) {
    # Strip scale
    scales <- g$scales$scales
    if (length(scales) == 4) {
      scales[[2]] <- NULL
      scales[[1]]$palette <- scales::gradient_n_pal(
        colours = c("#7b3294", "#c2a5cf", "#f7f7f7", "#a6dba0", "#008837")
      )
    } else if (length(scales) == 3) {
      scales[[1]] <- NULL
    }
    g$scales$scales <- scales

    # Add new scales if non-raw output is desired
    if (!raw) {
      g <- g + ggplot2::scale_fill_gradientn(
        colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
        name = expression(frac("Observed", "Expected")),
        limits = c(-1, 1) * max(abs(g$data$value[g$data$mode == "Individual"] - 1)) + 1,
        guide = guide_colourbar(order = 1)
      )
    }
  }

  g
}
