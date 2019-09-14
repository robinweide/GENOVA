# Documentation -----------------------------------------------------------

#' @name visualise
#' @title Visualise discoveries
#'
#' @description Plot the results of \code{discovery} objects. By default
#'   contrasts one sample with the others.
#'
#' @param discovery A \code{discovery} object as returned by GENOVA analysis
#'   functions.
#' @param contrast An \code{integer} or \code{character} matching an experiment
#'   (name) of length 1, specifying which sample should be used to make a
#'   contrast with all other samples. Alternatively, set to \code{NULL} to not
#'   plot contrast panels.
#' @param metric A \code{character} of length 1: what contrast metric should be
#'   used? \code{"diff"} for difference by subtraction or \code{"lfc"} for
#'   \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} fold changes.
#' @param raw A \code{logical} of length 1: should a bare bones plot be
#'   returned?
#' @param mode (PESCAn_discovery only) What result slot should be used for
#'   visualisation? \code{"obsexp"} for the observed over expected metric or
#'   \code{"signal"} for mean contacts at unshifted anchors.
#'
#' @details The \code{"diff"} \code{metric} value creates contrast panels by
#'   subtracting the values of each sample by the values of the sample indicated
#'   by the '\code{contrast}' argument. The \code{"lfc"} \code{metric} value
#'   creates contrast panels by dividing every samples' values by the sample
#'   indicated by the '\code{contrast}' argument and then taking the base 2
#'   logarithm of that division.
#'
#'   The '\code{raw = TRUE}' argument allows custimisation of the plots. The
#'   returned plot will not have position- or fill-scales and no theme settings,
#'   which can be set freely afterwards to match personal aesthetic tastes. When
#'   '\code{raw = TRUE}' and '\code{subtract}' is not \code{NULL}, the fill
#'   scale of the contrast panels can be manipulated by setting the
#'   '\code{aesthetics = "altfill"}' inside ggplot2's fill scale functions.
#'
#' @note For \code{ATA_discovery} and \code{ARA_discovery} objects which include
#'   the matrix's diagonal, the upper limit for the contacts fill scale is set to
#'   the 95th percentile of the data to increase the dynamic range of colours.
#'
#' @examples
#' # APA
#' apa <- APA(list(WT = WT_40kb, KO = KO_40kb), loops)
#' visualise(apa)
#'
#' # PE-SCAn
#' pescan <- PESCAn(list(WT = WT_40kb, KO = KO_40kb), super_enhancers)
#' visualise(pescan)
#'
#' # To plot PE-SCAn without background correction
#' visualise(pescan, mode = "signal")
#'
#' # ATA
#' ata <- APA(list(WT = WT_10kb, KO = KO_10kb), tads)
#' visualise(ata)
#'
#' # ARA
#' ara <- ARA(list(WT = WT_20kb, KO = KO_20kb), ctcf_sites)
#' visualise(ara)
#'
#' # Handling 'raw' plots
#' visualise(pescan, raw = TRUE) +
#'   ggplot2::scale_fill_gradient(aesthetics = "altfill")
NULL

# Default -----------------------------------------------------------------

# If you fail, it is best to fail graciously

#' @export
#' @rdname visualise
#' @usage NULL
visualise.default <- function(discovery, ...) {
  stop("No visualise method for class '", class(discovery),
       "' has been implemented.", call. = FALSE)
}


# Common elements ---------------------------------------------------------

# Common ancestor for aggregate repeated matrix lookup analysis plots
visualise.ARMLA <- function(discovery, contrast = 1,
                            metric = c("diff", "lfc"),
                            raw = FALSE, altfillscale) {
  metric <- match.arg(metric)
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

  # Calculate metric if contrast is supplied
  if (!is.null(contrast)) {
    contrast <- as.vector(mats[,,contrast])
    contrast <- switch(
      metric,
      "diff" = as.vector(mats) - rep(contrast, dims[3]),
      "lfc" = log2(as.vector(mats) / rep(contrast, dims[3]))
    )
    dim(contrast) <- dim(mats)
    contrast_name <- switch(metric,
                            "diff" = "Difference",
                            "lfc" = "Change")
    contrast <- lapply(seq_len(dims[3]), function(i) {
      cbind.data.frame(coords,
                       name = names_dims[[3]][i],
                       mode = contrast_name,
                       value = as.vector(contrast[,,i]))
    })
    contrast <- do.call(rbind, contrast)
    outcontrast <<- contrast
    df <- rbind(df, contrast)
  }

  # Setup base of plot
  g <- ggplot2::ggplot(df, ggplot2::aes(x, y))

  if (!is.null(contrast)) {
    # Setup basics of the diff plots
    # Warnings are supressed because ggplot doesn't recognise
    # the altfill aesthetic (YET!)
    suppressWarnings(
      g <- g + ggplot2::geom_raster(
        data = function(x) {
          x[x$mode == contrast_name, ]
        },
        ggplot2::aes(altfill = value)
      ) +
        ggplot2::facet_grid(mode ~ name, switch = "y")
    )

    if (!raw) {
      g <- g + altfillscale
      # Hack for scale
      g$scales$scales[[1]]$guide <- ggplot2::guide_colourbar()
      g$scales$scales[[1]]$guide$available_aes[[3]] <- "altfill"
      g$scales$scales[[1]]$guide$order <- 2
    } else {
      g <- g + ggplot2::guides(
        "altfill" = ggplot2::guide_colourbar(available_aes = "altfill")
      )
    }

    # Good old ggnomics-style hack for different fill scales
    old_geom <- g$layers[[1]]$geom
    old_nahandle <- old_geom$handle_na
    new_nahandle <- function(self, data, params) {
      colnames(data)[colnames(data) %in% "altfill"] <- "fill"
      old_nahandle(data, params)
    }
    new_geom <- ggplot2::ggproto(paste0(sample(1e6, 1), class(old_geom)),
                                 old_geom,
                                 handle_na = new_nahandle
    )
    names(new_geom$default_aes)[1] <- "altfill"
    new_geom$non_missing_aes <- "altfill"
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
  } else {
    g <- g + ggplot2::theme(
      aspect.ratio = 1,
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey30")
    )
    return(g)
  }
}

# Discovery objects -------------------------------------------------------

#' @rdname visualise
#' @export
visualise.APA_discovery <- function(discovery, contrast = 1,
                                    metric = c("diff", "lfc"),
                                    raw = FALSE) {
  metric <- match.arg(metric)

  altfillscale <- ggplot2::scale_fill_gradientn(
    colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
    aesthetics = "altfill",
    name = switch (metric,
      "diff" = "Difference",
      "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    limits = function(x){c(-1, 1) * max(abs(x))}
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = discovery,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale
  )
  if (raw) {
    return(g)
  }

  pos_breaks <- function(x) {
    x <- scales::extended_breaks()(x)
    if (length(x) > 3) head(tail(x, -1), -1) else x
  }

  g <- g + ggplot2::scale_fill_gradientn(
    colours = c('white', '#f5a623', '#d0021b', 'black'),
    guide = ggplot2::guide_colourbar(order = 1),
    name = expression(mu*" Contacts")
  ) +
  ggplot2::scale_x_continuous(
    name = "",
    expand = c(0, 0),
    breaks = pos_breaks,
    labels = function(x) {
      ifelse(x == 0, "3'", paste0(x / 1000, "kb"))
    }
  ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        ifelse(x == 0, "5'", paste0(x / 1000, "kb"))
      }
    )

  g
}

#' @rdname visualise
#' @export
visualise.PESCAn_discovery <- function(discovery, contrast = 1,
                                       metric = c("diff", "lfc"),
                                       mode = c("obsexp", "signal"),
                                       raw = FALSE) {
  metric <- match.arg(metric)
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

  altcols <- if (hasobsexp) {
    c("#23AA17", "#90D48E", "#FFFFFF", "#F0A9F1", "#DA64DC")
  } else {
    c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49")
  }
  altfillscale <- ggplot2::scale_fill_gradientn(
    colours = altcols,
    aesthetics = "altfill",
    name = switch (metric,
                   "diff" = "Difference",
                   "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    limits = function(x){c(-1, 1) * max(abs(x))}
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = res,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale
  )
  if (raw) {
    return(g)
  }

  fillscale <- if (hasobsexp) {
    ggplot2::scale_fill_gradientn(
      colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(frac("Observed", "Expected")),
      limits = function(x) {
        c(-1, 1) * max(abs(x - 1)) + 1
      }
    )
  } else {
    ggplot2::scale_fill_gradientn(
      colours = c('white', '#f5a623', '#d0021b', 'black'),
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(mu*" Contacts")
    )
  }

  pos_breaks <- function(x) {
    x <- scales::extended_breaks()(x)
    if (length(x) > 3) head(tail(x, -1), -1) else x
  }

  g <- g + fillscale +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        ifelse(x == 0, "3'", paste0(x / 1000, "kb"))
      }
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        ifelse(x == 0, "5'", paste0(x / 1000, "kb"))
      }
    )

  g
}

#' @rdname visualise
#' @export
visualise.ATA_discovery <- function(discovery, contrast = 1,
                                    metric = c("diff", "lfc"),
                                    raw = FALSE) {
  metric <- match.arg(metric)

  altfillscale <- ggplot2::scale_fill_gradientn(
    colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
    aesthetics = "altfill",
    name = switch (metric,
                   "diff" = "Difference",
                   "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    limits = function(x){c(-1, 1) * max(abs(x))}
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = discovery,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale
  )
  if (raw) {
    return(g)
  }

  haspad <- "padding" %in% names(attributes(discovery))
  if (haspad) {
    padding <- attr(discovery, "padding")
    pos_breaks <- function(x) {
      cumsum(c(padding - 0.5, 1)) / (2 * padding) * diff(x) + min(x)
    }
  } else {
    pos_breaks <- function(x) {
      seq(x[1], x[2], length.out = 5)[c(2,4)]
    }
  }

  upperq <- quantile(discovery$signal, 0.95)

  g <- g +
    ggplot2::scale_fill_gradientn(
      colours = c('white', '#f5a623', '#d0021b', 'black'),
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(mu*" Contacts"),
      limits = c(NA, upperq),
      oob = scales::squish
    ) +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = c("5' border", "3' border")
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = c("5' border", "3' border")
    )

  g
}

#' @rdname visualise
#' @export
visualise.ARA_discovery <- function(discovery, contrast = 1,
                                    metric = c("diff", "lfc"),
                                    raw = FALSE) {
  metric <- match.arg(metric)

  altfillscale <- ggplot2::scale_fill_gradientn(
    colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
    aesthetics = "altfill",
    name = switch (metric,
                   "diff" = "Difference",
                   "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    limits = function(x){c(-1, 1) * max(abs(x))}
  )

  # Get a default plot
  g <- GENOVA:::visualise.ARMLA(
    discovery = discovery,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale
  )
  if (raw) {
    return(g)
  }

  pos_breaks <- function(x) {
    x <- scales::extended_breaks()(x)
    if (length(x) > 3) head(tail(x, -1), -1) else x
  }

  upperq <- quantile(discovery$signal, 0.95)

  g <- g + ggplot2::scale_fill_gradientn(
    colours = c('white', '#f5a623', '#d0021b', 'black'),
    guide = ggplot2::guide_colourbar(order = 1),
    name = expression(mu*" Contacts"),
    limits = c(NA, upperq),
    oob = scales::squish
  ) +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        ifelse(x == 0, "3'", paste0(x / 1000, "kb"))
      }
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        ifelse(x == 0, "5'", paste0(x / 1000, "kb"))
      }
    )

  g
}

# Utilities ---------------------------------------------------------------

# Makes sure no errors are returned when visualise(..., raw = TRUE)
# Unfortunately has to be exported, but is better than throwing errors
#' @export
#' @keywords internal
scale_altfill_continuous <- function(...) {
  ggplot2::scale_fill_gradient2(..., aesthetics = "altfill")
}
