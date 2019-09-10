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
#' @details Difference panels are created by subtracting the values of each
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
visualise.APA_discovery <- function(discovery, subtract = 1, raw = FALSE) {

  mats <- discovery$signal

  # Format dimnames if missing
  dims <- dim(mats)
  names_dims <- dimnames(mats)
  names_dims <- lapply(seq_along(dims), function(i) {
    if (is.null(names_dims[[i]])) {
      return(seq_len(dims[i]))
    } else {
      return(names_dims[[i]])
    }
  })

  # Setup coordinates
  coords <- data.frame(
    x = as.numeric(names_dims[[2]][as.vector(col(mats[,,1]))]),
    y = as.numeric(names_dims[[1]][as.vector(row(mats[,,1]))])
  )

  # Fill coordinates with experiment data
  df <- lapply(seq_len(dims[3]), function(i) {
    cbind.data.frame(coords,
                     name = names_dims[[3]][i],
                     mode = "Individual",
                     value = as.vector(mats[,,i]))
  })
  df <- do.call(rbind, df)

  # Calculate diff data if subtract is supplied
  if (!is.null(subtract)) {
    diff <- as.vector(mats[,,subtract])
    diff <- as.vector(mats) - rep(diff, dims[3])
    dim(diff) <- dim(mats)
    diff <- lapply(seq_len(dims[3]), function(i) {
      cbind.data.frame(coords,
                       name = names_dims[[3]][i],
                       mode = "Difference",
                       value = as.vector(diff[,,i]))
    })
    diff <- do.call(rbind, diff)
    df <- rbind(df, diff)
  }

  # Setup base of plot
  g <- ggplot2::ggplot(df, ggplot2::aes(x, y))

  if (!is.null(subtract)) {
    # Setup basics of the diff plots
    # Warnings are supressed because ggplot doesn't recognise
    # the fill2 aesthetic (YET!)
    suppressWarnings(
      g <- g + ggplot2::geom_raster(
        data = function(x) {
          x[x$mode == "Difference", ]
        },
        ggplot2::aes(fill2 = value)
      ) +
        ggplot2::facet_grid(mode ~ name, switch = "y")
    )

    if (!raw) {
      g <- g + ggplot2::scale_fill_gradientn(
        colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
        aesthetics = "fill2",
        name = "Difference",
        limits = c(-1, 1) * max(abs(df$value[df$mode == "Difference"]))
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
    g <- g + ggplot2::facet_grid(~ name)
  }

  # Add the non-diff plots
  g <- g + ggplot2::geom_raster(
    data = function(x) {
      x[x$mode == "Individual", ]
    },
    ggplot2::aes(fill = value)
  )

  if (raw) {
    return(g)
  }

  # Decorating of the plot
  g <- g + ggplot2::scale_fill_gradientn(
    colours = c("white", "orange", "red", "black"),
    guide = ggplot2::guide_colourbar(order = 1),
    name = expression(mu*" Contacts")
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
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey30")
    )
  g
}

#' @export
#' @rdname APA_PESCAn_visualisation
visualise.PESCAn_discovery <- function(discovery,
                                       subtract = 1,
                                       mode = c("obsexp", "signal"),
                                       raw = FALSE) {
  # Handle mode settings
  mode <- match.arg(mode)
  hasobsexp <- "obsexp" %in% names(discovery)
  if (mode == "obsexp" & !hasobsexp) {
    warning("Mode was set to 'obsexp' but no such result was found")
  }
  hasobsexp <- hasobsexp && mode == "obsexp"
  res <- if (hasobsexp) {
    setNames(discovery[names(discovery) %in% "obsexp"], "signal")
  } else {
    discovery
  }

  # Cannibalise the APA autoplot function
  g <- visualise.APA_discovery(res, subtract = subtract, raw = raw)

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
        guide = ggplot2::guide_colourbar(order = 1)
      )
    }
  }

  g
}
