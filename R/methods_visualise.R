# Documentation -----------------------------------------------------------

#' @name visualise
#' @aliases visualise visualize
#' @title Visualise discoveries
#'
#' @description Plot the results of \code{discovery} objects. By default
#'   contrasts one sample with the others.
#'
#' @param discovery A \code{discovery} object as returned by GENOVA analysis
#'   functions.
#'
#' @param contrast An \code{integer} or \code{character} matching an experiment
#'   (name) of length 1, specifying which sample should be used to make a
#'   contrast with all other samples. Alternatively, set to \code{NULL} to not
#'   plot contrast panels. See also the \code{show_single_contrast} argument.
#'
#' @param colour_lim,colour_lim_contrast
#'   Indication of limits for the primary and secondary continuous colour scales 
#'   respectively. One of: \itemize{ 
#'   \item \code{NULL} to have an educated quess for the scale. 
#'   \item A \code{numeric} vector of length two providing the limits of the 
#'   scale. Use \code{NA} to refer to existing minima or maxima. 
#'   \item A \code{function} that accepts the existing (automatic) limits and 
#'   returns new limits.}
#'
#' @param raw A \code{logical} of length 1: should a bare bones plot be
#'   returned?
#'
#' @param title add a title
#' @param ... Further arguments specific to the discovery class. See the section 
#' extended arguments below.
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
#' @section Extended arguments:
#'
#' @note For \code{ATA_discovery} objects which include the matrix's diagonal,
#'   the upper limit for the contacts fill scale is set to the 95th percentile
#'   of the data to increase the dynamic range of colours. The same is true for
#'   \code{ARA_discovery}, when '\code{mode = "signal"}'.
#'   
#' @export
#'
#' @examples
#' \dontrun{
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
#' # Compartment score
#' cs <- compartment_score(list(WT = WT_100kb, KO = KO_100kb), H3K4me1_peaks)
#' visualise(cs, chr = "chr1")
#'
#' # Saddle function
#' sadl <- saddle(list(WT = WT_100kb, KO = KO_100kb), cs) # see example above
#' visualise(sadl)
#'
#' # Handling 'raw' plots
#' visualise(pescan, raw = TRUE) +
#'   ggplot2::scale_fill_gradient(aesthetics = "altfill")
#' }
visualise <- function(discovery, ...) {
  UseMethod("visualise", discovery)
}



# Default -----------------------------------------------------------------

# If you fail, it is best to fail graciously

#' @export
#' @rdname visualise
visualise.default <- function(discovery, contrast, raw, title, 
                              colour_lim, colour_lim_contrast, ...) {
  stop("No visualise method for class '", class(discovery),
       "' has been implemented.", call. = FALSE)
}

# Common elements ---------------------------------------------------------

# Common ancestor for aggregate repeated matrix lookup analysis plots
visualise.ARMLA <- function(discovery, contrast = 1,
                            metric = c("diff", "lfc"),
                            raw = FALSE, altfillscale, 
                            show_single_contrast = FALSE,
                            ...) {
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
  showcontrast <- !is.null(contrast) && 
    (dims[3] > 1L || literalTRUE(show_single_contrast))
  
  if (showcontrast) {
    contrast <- as.vector(mats[,,contrast])
    contrast <- switch(
      metric,
      "diff" = as.vector(mats) - rep(contrast, dims[3]),
      "lfc"  = log2(as.vector(mats) / rep(contrast, dims[3]))
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
    df <- rbind(df, contrast)
  }

  # Setup base of plot
  g <- ggplot2::ggplot(df, ggplot2::aes(x, y))

  if (showcontrast) {
    # Setup basics of the diff plots
    # Warnings are supressed because ggplot doesn't recognise
    # the altfill aesthetic (YET!)
    suppressWarnings(
      g <- g + ggplot2::geom_raster(
        data = datafilter(mode == contrast_name),
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
    data = datafilter(mode == "Individual"),
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

# Aggregate matrices ------------------------------------------------------

#' @rdname visualise
#' @section Extended arguments:\subsection{APA, PE-SCAn, ATA, ARA & C-SCAn}{
#' \describe{
#'  \item{\code{metric}}{A \code{character} of length one: \code{"diff"} for 
#'   difference by subtraction or \code{"lfc"} for 
#'   \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} fold changes.}
#'  \item{\code{colour_lim}, \code{colour_lim_contrast}}{
#'   Indication of limits for the primary and secondary colour scale 
#'   respectively. One of: \itemize{ 
#'   \item \code{NULL} to have an educated quess for the scale. 
#'   \item A \code{numeric} vector of length two providing the limits of the 
#'   scale. Use \code{NA} to refer to existing minima or maxima. 
#'   \item A \code{function} that accepts the existing (automatic) limits and 
#'   returns new limits.}}
#'  \item{\code{mode}}{A \code{character} of length one indicating what type of
#'  result to plot. Either \code{"signal"} or \code{"obsexp"}, referring to the
#'  slot in the discovery object. Applicatble to PE-SCAn, C-SCAn and ARA.}
#'  \item{\code{show_single_contrast}}{A \code{logical} of length 1; 
#'   if \code{FALSE} (default), does not show contrasts when \code{discovery} 
#'   describes one experiment. If \code{TRUE}, plots empty panel.}
#' }
#' }
#' @export
#' @usage NULL
visualise.APA_discovery <- function(discovery, contrast = 1,
                                    metric = c("lfc", "diff"),
                                    raw = FALSE, title = NULL,
                                    colour_lim = NULL,
                                    colour_lim_contrast = NULL, 
                                    show_single_contrast = FALSE,
                                    ...) {
  metric <- match.arg(metric)
  
  # Decide on limits
  if (is.null(colour_lim)) {
    colour_lim <- c(NA, NA)
  }

  altfillscale <- scale_fill_GENOVA_div(
    aesthetics = "altfill",
    name = switch (metric,
      "diff" = "Difference",
      "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    limits = colour_lim_contrast,
    oob = scales::squish,
    midpoint = 0
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = discovery,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale,
    show_single_contrast = show_single_contrast
  )
  if (raw) {
    return(g)
  }

  g <- g + scale_fill_GENOVA(
    guide = ggplot2::guide_colourbar(order = 1),
    name = expression(mu*" Contacts"),
    limits = colour_lim,
    oob = scales::squish
  ) +
  ggplot2::scale_x_continuous(
    name = "",
    expand = c(0, 0),
    breaks = breaks_trim_outer(),
    labels = label_kilobase_relative("3'")
  ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = breaks_trim_outer(),
      labels = label_kilobase_relative("5'")
    )

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  g
}

#' @rdname visualise
#' @export
#' @usage NULL
visualise.CSCAn_discovery <- function(discovery, mode = c("obsexp", "signal"),
                                      raw = FALSE, title = NULL, 
                                      colour_lim = NULL,
                                      ...) {
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
  
  df <- c(lapply(seq_along(dim(res$signal)), function(i) {
    dnm <- dimnames(res$signal)[[i]]
    if (is.null(dnm)) {
      return(as.vector(slice.index(res$signal, i)))
    } else {
      return(dnm[as.vector(slice.index(res$signal, i))])
    }
  }), list(as.vector(res$signal)))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- setNames(df, c(paste0("Var", seq_along(dim(res$signal))), "value"))
  if (ncol(df) == 5) {
    groups <- strsplit(as.character(df$Var3), "-")
    df$left  <- vapply(groups, `[`, character(1), 1)
    df$right <- vapply(groups, `[`, character(1), 2)
  }
  df$Var1 <- as.integer(df$Var1)
  df$Var2 <- as.integer(df$Var2)
  
  g <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = value)) +
    ggplot2::geom_raster()
  
  if (all(c("left", "right") %in% names(df))) {
    g <- g + ggplot2::facet_grid(left ~ Var4 + right, switch = "y")
  } else if ("Var4" %in% names(df)) {
    g <- g + ggplot2::facet_grid(~ Var4, switch = "y")
  }
  
  if (raw) {
    return(g)
  }

  panel_spac <- c(rep(5.5, pmax(length(unique(df$right)) - 1L, 0)))
  panel_spac <- c(panel_spac, rep(c(11, panel_spac), 
                                  pmax(length(unique(df$Var4)) - 1, 0)))
  panel_spac <- if (length(panel_spac)) panel_spac else 5.5
  panel_spac <- ggplot2::unit(panel_spac, "points")
  
  fillscale <- if (hasobsexp) {
    scale_fill_GENOVA_div(
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(frac("Observed", "Expected")),
      oob = scales::squish,
      midpoint = 1,
      limits = colour_lim
    )
  } else {
    scale_fill_GENOVA(
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(mu*" Contacts"),
      oob = scales::squish,
      limits = colour_lim
    )
  }

  g <- g + fillscale + 
    ggplot2::scale_x_continuous(breaks = breaks_trim_outer(), 
                                labels = label_kilobase_relative("3'"), 
                                expand = c(0,0),
                                name = "") +
    ggplot2::scale_y_continuous(breaks = breaks_trim_outer(),
                                labels = label_kilobase_relative("5'"), 
                                expand = c(0,0), name = "") +
    GENOVA_THEME() +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.placement = "outside",
      panel.spacing.x = panel_spac
    )
  
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }
  
  g
}

#' @rdname visualise
#' @export
#' @usage NULL
visualise.PESCAn_discovery <- function(discovery, contrast = 1,
                                       metric = c("diff", "lfc"),
                                       mode = c("obsexp", "signal"),
                                       raw = FALSE, title = NULL,
                                       colour_lim = NULL,
                                       colour_lim_contrast = NULL, 
                                       show_single_contrast = FALSE,
                                       ...) {
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
  
  # Decide on limits
  # if (is.null(colour_lim)) {
  if (hasobsexp) {
    midpoint <- 1
  } else {
    midpoint <- NA
  }
  # }

  altcols <- if (hasobsexp) {
    "greenpink"
  } else {
    "divergent"
  }
  altfillscale <- scale_fill_GENOVA_div(
    palette = altcols,
    aesthetics = "altfill",
    name = switch (metric,
                   "diff" = "Difference",
                   "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    midpoint = 0,
    limits = colour_lim_contrast,
    oob = scales::squish
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = res,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale,
    show_single_contrast = show_single_contrast
  )
  if (raw) {
    return(g)
  }

  fillscale <- if (hasobsexp) {
    scale_fill_GENOVA_div(
      palette = "divergent",
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(frac("Observed", "Expected")),
      oob = scales::squish,
      limits = colour_lim,
      midpoint = 1
    )
  } else {
    scale_fill_GENOVA(
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(mu*" Contacts"),
      oob = scales::squish,
      limits = colour_lim
    )
  }

  g <- g + fillscale +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = breaks_trim_outer(),
      labels = label_kilobase_relative("3'")
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = breaks_trim_outer(),
      labels = label_kilobase_relative("5'")
    )

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  
  g
}

#' @rdname visualise
#' @export
#' @usage NULL
visualise.ATA_discovery <- function(discovery, contrast = 1,
                                    metric = c("lfc", "diff"),
                                    raw = FALSE, title = NULL,
                                    colour_lim = NULL,
                                    colour_lim_contrast = NULL, 
                                    show_single_contrast = FALSE,
                                    ...) {
  metric <- match.arg(metric)
  
  # Decide on limits
  # if (is.null(colour_lim_contrast)) {
  #   colour_lim_contrast <- centered_limits()
  # }

  altfillscale <- scale_fill_GENOVA_div(
    aesthetics = "altfill",
    name = switch (metric,
                   "diff" = "Difference",
                   "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    oob = scales::squish,
    midpoint = 0,
    limits = colour_lim_contrast
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = discovery,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale,
    show_single_contrast = show_single_contrast
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
  
  if (is.null(colour_lim)) {
    colour_lim <- c(NA, quantile(discovery$signal, 0.95))
  }

  g <- g +
    scale_fill_GENOVA(
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(mu*" Contacts"),
      limits = colour_lim,
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
      labels = c("3' border", "5' border")
    )

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  
  g
}

#' @rdname visualise
#' @export
#' @usage NULL
visualise.ARA_discovery <- function(discovery, contrast = 1,
                                    metric = c("diff", "lfc"),
                                    mode = c("obsexp", "signal"),
                                    raw = FALSE, title = NULL,
                                    colour_lim = NULL,
                                    colour_lim_contrast = NULL, 
                                    show_single_contrast = FALSE,
                                    ...) {
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
    pal <- "greenpink"
  } else {
    pal <- "divergent"
  }
  altfillscale <- scale_fill_GENOVA_div(
    palette = pal,
    aesthetics = "altfill",
    name = switch (metric,
                   "diff" = "Difference",
                   "lfc" = expression(atop("Log"[2]*" Fold", "Change"))
    ),
    limits = colour_lim_contrast,
    midpoint = 0,
    oob = scales::squish
  )

  # Get a default plot
  g <- visualise.ARMLA(
    discovery = res,
    contrast = contrast,
    metric = metric, raw = raw,
    altfillscale = altfillscale,
    show_single_contrast = show_single_contrast
  )
  if (raw) {
    return(g)
  }

  # pos_breaks <- function(x) {
  #   x <- scales::extended_breaks()(x)
  #   if (length(x) > 3) head(tail(x, -1), -1) else x
  # }
  
  # Decide on limits
  if (is.null(colour_lim)) {
    if (!hasobsexp) {
      colour_lim <- c(NA, quantile(discovery$signal, 0.95))
    }
  }

  fillscale <- if (hasobsexp) {
    scale_fill_GENOVA_div(
      palette = "divergent",
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(frac("Observed", "Expected")),
      oob = scales::squish,
      limits = colour_lim,
      midpoint = 1
    )
  } else {
    scale_fill_GENOVA(
      guide = ggplot2::guide_colourbar(order = 1),
      name = expression(mu*" Contacts"),
      limits = colour_lim,
      oob = scales::squish
    )
  }

  g <- g + fillscale +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = breaks_trim_outer(),
      labels = scales::label_number(scale = 1e-3, suffix = " kb")
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = breaks_trim_outer(),
      labels = scales::label_number(scale = 1e-3, suffix = " kb")
    )

  if(!is.null(title)){
    g <- g + ggplot2::ggtitle(title)
  }
  
  g
}

# Genome wide scores ------------------------------------------------------

#' @rdname visualise
#' @section Extended arguments:\subsection{Compartment-, Insulation score & 
#' Directionality Index}{
#' \describe{
#'  \item{\code{chr}}{A \code{character} of length one with the chromosome 
#'  name. Defaults to \code{"chr1"}.}
#'  \item{\code{start}, \code{end}}{An \code{integer} of length one setting the 
#'  start or end of the region to plot. If \code{NULL} (default), is set to 
#'  \code{-Inf} and \code{Inf} respectively}
#' }
#' }
#' @export
#' @usage NULL
visualise.CS_discovery <- function(discovery, contrast = NULL,
                                   chr = "chr1", start = NULL, end = NULL,
                                   raw = FALSE, show_single_contrast = FALSE,
                                   ...) {
  start <- if (is.null(start)) -Inf else start
  end <- if (is.null(end)) Inf else end
  loc <- standardise_location(chr, start, end, singular = TRUE)
  df <- as.data.table(discovery$compart_scores)
  ii <- which(df[["chrom"]] == loc$chrom & df[["start"]] >= loc$start & 
                df[["end"]] <= loc$end)
  df <- df[ii,]
  expnames <- colnames(df)[5:ncol(df)]
  
  df <- data.frame(mid = (df[["start"]] + df[["end"]])/2,
                   df[, ..expnames])
  
 showcontrast <- !is.null(contrast) && (length(expnames) > 1L || literalTRUE(show_single_contrast))
  if (showcontrast) {
    cdf <- as.matrix(df[,expnames])
    cdf <- cdf - cdf[, contrast]
    cdf <- cbind.data.frame(mid = df$mid, cdf)
  } else {
    cdf <- NULL
  }
  
  yname <- if (attr(discovery, "signed")) "Compartment Score" else "Score"
  
  # Melt df
  df <- data.frame(mid = rep(df$mid, length(expnames)),
                   score = unlist(df[2:ncol(df)]),
                   exp = rep(expnames, each = nrow(df)))
  
  if (showcontrast) {
    # Melt cdf
    cdf <- data.frame(mid = rep(cdf$mid, length(expnames)),
                      score = unlist(cdf[2:ncol(cdf)]),
                      exp = rep(expnames, each = nrow(cdf)),
                      panel = "Difference")
    df$panel <- yname
  }
  
  
  rownames(df) <- NULL
  
  g <- ggplot2::ggplot(df, ggplot2::aes(mid, score, colour = exp)) +
    ggplot2::geom_line()
  
  if (!is.null(cdf)) {
    g <- g + ggplot2::geom_line(data = cdf) +
      ggplot2::facet_grid(panel ~ ., scales = "free_y", switch = "y")
  }
  
  if (raw) {
    return(g)
  }
  cols <- attr(discovery, "colours")
  g <- g + ggplot2::scale_x_continuous(
    name = paste0("Location ", chr),
    expand = c(0.01,0),
    labels = scales::label_number(scale = 1e-6, suffix = "Mb")
    ) +
    ggplot2::scale_y_continuous(
      name = yname,
      # limits = c(-1.25, 1.25),
      oob = scales::squish,
      limits = centered_limits()
    ) +
    ggplot2::scale_colour_manual(
      name = "Sample",
      breaks = expnames,
      limits = expnames,
      values = cols
    )
  g <- g + ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    legend.key = ggplot2::element_blank(),
    aspect.ratio = 2 / (1 + sqrt(5)),
    text = ggplot2::element_text(color = 'black'),
    axis.text = ggplot2::element_text(color = 'black'),
  )
  if (!is.null(cdf)) {
    g <- g + ggplot2::theme(
      strip.placement = "outside",
      axis.title.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = ggplot2::rel(1))
    )
  }
  g
}

#' @rdname visualise
#' @export
#' @usage NULL
visualise.IS_discovery <- function(discovery, contrast = NULL, chr = "chr1",
                                   start = NULL, end = NULL, raw = FALSE, 
                                   show_single_contrast = FALSE,
                                   ...) {
  start <- if (is.null(start)) -Inf else start
  end <- if (is.null(end)) Inf else end
  loc <- standardise_location(chr, start, end, singular = TRUE)
  df <- as.data.table(discovery$insula_score)
  ii <- which(df[["chrom"]] == loc$chrom & df[["start"]] >= loc$start & 
                df[["end"]] <= loc$end)
  df <- df[ii,]
  expnames <- colnames(df)[5:ncol(df)]
  
  df <- data.table(mid = (df[["start"]] + df[["end"]])/2,
                   df[, ..expnames])

  showcontrast <- !is.null(contrast) && 
    (length(expnames) > 1L || literalTRUE(show_single_contrast))
  
  if (showcontrast) {
    cdf <- as.matrix(df[,..expnames])
    cdf <- cdf - cdf[, contrast]
    cdf <- cbind.data.frame(mid = df$mid, cdf)
  } else {
    cdf <- NULL
  }
  
  yname <- "Insulation Score"
  
  # Melt df
  df <- data.frame(mid = rep(df$mid, length(expnames)),
                   score = unlist(as.list(df)[2:ncol(df)]),
                   exp = rep(expnames, each = nrow(df)))
  
  if (showcontrast) {
    # Melt cdf
    cdf <- data.frame(mid = rep(cdf$mid, length(expnames)),
                      score = unlist(cdf[2:ncol(cdf)]),
                      exp = rep(expnames, each = nrow(cdf)),
                      panel = "Difference")
    df$panel <- yname
  }
  
  
  rownames(df) <- NULL
  
  g <- ggplot2::ggplot(df, ggplot2::aes(mid, score, colour = exp)) +
    ggplot2::geom_line()
  
  if (!is.null(cdf)) {
    g <- g + ggplot2::geom_line(data = cdf) +
      ggplot2::facet_grid(panel ~ ., scales = "free_y", switch = "y")
  }
  
  if (raw) {
    return(g)
  }
  cols <- attr(discovery, "colours")
  g <- g + ggplot2::scale_x_continuous(
    name = paste0("Location ", chr),
    expand = c(0.01,0),
    labels = scales::label_number(scale = 1e-6, suffix = "Mb")
  ) +
    ggplot2::scale_y_continuous(
      name = yname,
      oob = scales::squish,
      limits = symmetric_quantiles(df$score, c(0.05, 0.95))
    ) +
    ggplot2::scale_colour_manual(
      name = "Sample",
      breaks = expnames,
      limits = expnames,
      values = cols
    )
  g <- g + ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    legend.key = ggplot2::element_blank(),
    aspect.ratio = 2 / (1 + sqrt(5)),
    text = ggplot2::element_text(color = 'black'),
    axis.text = ggplot2::element_text(color = 'black'),
  )
  if (!is.null(cdf)) {
    g <- g + ggplot2::theme(
      strip.placement = "outside",
      axis.title.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = ggplot2::rel(1))
    )
  }
  g
}

#' @rdname visualise
#' @export
#' @usage NULL
visualise.DI_discovery <-  function(discovery, contrast = NULL, chr = "chr1",
                                    start = NULL, end = NULL, raw = FALSE, 
                                    show_single_contrast = FALSE,
                                    ...) {
  start <- if (is.null(start)) -Inf else start
  end <- if (is.null(end)) Inf else end
  loc <- standardise_location(chr, start, end, singular = TRUE)
  df <- discovery$DI
  ii <- which(df[["chrom"]] == loc$chrom & df[["start"]] >= loc$start & 
                df[["end"]] <= loc$end)
  df <- df[ii,]
  expnames <- tail(colnames(df), -4)
  
  df <- cbind.data.frame(mid = (df[["start"]] + df[["end"]]/2),
                         df[, expnames])
  
  yname <- "Directionality Index"
  
  showcontrast <- !is.null(contrast) && (length(expnames) > 1L || literalTRUE(show_single_contrast))
  
  if (showcontrast) {
    
    cdf <- do.call(cbind.data.frame, c(list(mid = df$mid), lapply(df[, -1], function(x){
      x - df[[contrast + 1]]
    })))
    setDT(cdf)
    cdf <- melt.data.table(cdf, value.name = "dir_index", id.vars = "mid")
    cdf[, panel := factor("Difference", levels = c(yname, "Difference"))]
    setDT(df)
    df <- melt.data.table(df, value.name = "dir_index", id.vars = "mid")
    df$panel <- factor(yname, levels = c(yname, "Difference"))
  } else {
    setDT(df)
    df <- melt.data.table(df, value.name = "dir_index", id.vars = "mid")
    cdf <- NULL
  }

  rownames(df) <- NULL
  
  g <- ggplot2::ggplot(df, ggplot2::aes(mid, dir_index, colour = variable)) +
    ggplot2::geom_line()
  
  if (!is.null(cdf)) {
    g <- g + ggplot2::geom_line(data = cdf) +
      ggplot2::facet_grid(panel ~ ., scales = "free_y", switch = "y")
  }
  
  if (raw) {
    return(g)
  }
  cols <- attr(discovery, "colours")
  g <- g + ggplot2::scale_x_continuous(
    name = paste0("Location ", chr),
    expand = c(0.01,0),
    labels = scales::label_number(scale = 1e-6, suffix = "Mb")
  ) +
    ggplot2::scale_y_continuous(
      name = yname,
      # limits = c(-1.25, 1.25),
      oob = scales::squish,
      limits = symmetric_quantiles(df$dir_index, c(0.005, 0.995))
    ) +
    ggplot2::scale_colour_manual(
      name = "Sample",
      breaks = expnames,
      limits = expnames,
      values = cols
    )
  g <- g + ggplot2::theme(
    panel.background = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black"),
    legend.key = ggplot2::element_blank(),
    aspect.ratio = 2 / (1 + sqrt(5)),
    text = ggplot2::element_text(color = 'black'),
    axis.text = ggplot2::element_text(color = 'black'),
  )
  if (!is.null(cdf)) {
    g <- g + ggplot2::theme(
      strip.placement = "outside",
      axis.title.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = ggplot2::rel(1))
    )
  }
  g
}

# Miscellaneous discoveries ----------------------------------------------------

#' @rdname visualise
#' @section Extended arguments:\subsection{Relative Contact Probability}{
#' \describe{
#'  \item{\code{metric}}{A \code{character} of length one. The choices are:
#'  \describe{
#'   \item{\code{"smooth"}}{A 
#'   \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}}-smoothed line 
#'   (default)}
#'   \item{\code{"both"}}{Like \code{"smooth"}, but also adds the raw binned 
#'   distances as points.}
#'   \item{\code{"lfc"}}{Displays 
#'   \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} fold change compared to 
#'   the sample specified in the \code{contrast} argument.}
#'  }}
#'  \item{\code{flipFacet}}{A \code{logical} of length one. If the 
#'  \code{bedlist} argument was provided to \code{RCP()}, combine all regions 
#'   in one panel? Defaults to \code{FALSE}.}
#' }
#' }
#' @export
#' @usage NULL
visualise.RCP_discovery = function(discovery, contrast = 1, 
                                   metric = c("smooth","both","lfc"), raw = FALSE, 
                                   title = NULL, flipFacet = FALSE, ...){
  
  metric <- match.arg(metric)
  
  nregion = length(unique(discovery$smooth$region))
  
  # colours
  smplCols = unique(discovery$smooth[,c('samplename', 'colour')])
  
  smplCols = lapply(split(smplCols,smplCols$colour), function(x){
    NL = nrow(x)
    if(NL > 1){
      
      cols = (col2rgb(unique(x$colour))+1)/255
      cols = sapply(1:NL, function(i){
        
        factor = ((1/(NL+2))*i)
        CF = cols+factor
        
        if(any(CF > 1)){
          factor = ((1/(NL+2))*i)
          CF = cols-factor
        }
        
        rgb(t(CF), maxColorValue = 1)
        
      })
      x$colour =  cols
    } else {
      x$colour =  rgb(t(col2rgb(unique(x$colour))/255))
    }
    x
  })
  
  smplCols = data.table::rbindlist(smplCols)
  smplCols = setNames(smplCols$colour,smplCols$samplename)
  
  D3cols <-c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", 
             "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78", 
             "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2", "#C7C7C7", 
             "#DBDB8D", "#9EDAE5")
  regCols = D3cols[1:nregion]
  
  
  ####################################################################### LFC
  
  if(metric == 'lfc'){
    smplevels = levels(discovery$smooth$samplename)
    if(is.numeric(contrast)){
      contrast = smplevels[contrast]
    } else if(!contrast %in% smplevels){
      stop('The contrast given is not found as samplename.')
    }
    
    # breaks
    logRange = round(log10(unique(discovery$smooth$distance)))
    logRange = range(logRange[is.finite(logRange)])
    logRange[2] = logRange[2] +1
    
    breaks <- 10**seq(from = logRange[1],
                      to =  logRange[2], length.out = 101)
    breaks = c(0,breaks)
    
    
    lfcDT = lapply(unique(discovery$smooth$region), function(REGION){
      
      tmp = RCPlfc(discovery$raw[discovery$raw$region == REGION,], contrast, breaks )
      tmp$region = REGION
      tmp
    })
    lfcDT = rbindlist(lfcDT)
    
    
    
    
    GG = NULL
    if(nregion == 1){
      GG = ggplot2::ggplot(lfcDT, ggplot2::aes(x = log10(distance), y = P ,col = samplename)) +
        ggplot2::labs(x = 'distance (Mb)', col = 'sample', y = expression("log2(P"[sample]*"/P"[contrast]*")")) +
        ggplot2::scale_color_manual(values = smplCols)
    } else if(flipFacet){
      GG =   ggplot2::ggplot(lfcDT, ggplot2::aes(x = log10(distance), y = P ,col = region)) +
        ggplot2::facet_grid(. ~ samplename)+
        ggplot2::labs(x = 'distance (Mb)', col = 'region', y = expression("log2(P"[sample]*"/P"[contrast]*")")) +
        ggplot2::scale_color_manual(values = regCols)
    } else {
      GG =ggplot2::ggplot(lfcDT, ggplot2::aes(x = log10(distance), y = P ,col = samplename)) +
        ggplot2::facet_grid(. ~ region)+
        ggplot2::labs(x = 'distance (Mb)', col = 'sample', y = expression("log2(P"[sample]*"/P"[contrast]*")"))+
        ggplot2::scale_color_manual(values = smplCols)
    }
    RAUW = GG +   
      ggplot2::geom_line() 
    
    # GG: misc
    breaks = unique(round(log10(unique(discovery$smooth$distance))))
    breaks = breaks[is.finite(breaks)]
    
    GG = GG + ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(breaks = breaks, labels = paste0((10**breaks)/1e6)) + 
      ggplot2::geom_line() +
      ggplot2::coord_fixed(ylim = range(lfcDT$P))
    
    GG = GG +  ggplot2::theme(panel.background = ggplot2::element_blank(),
                              aspect.ratio = 1,
                              strip.background = ggplot2::element_rect(fill = NA, colour = NA),
                              panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                              text = ggplot2::element_text(color = 'black'),
                              axis.line = ggplot2::element_blank(),
                              axis.text = ggplot2::element_text(colour = 'black'),
                              strip.text = ggplot2::element_text(colour = 'black') )
    
    if(!is.null(title)){
      GG = GG + ggplot2::ggtitle(title)
    }
    
    if(raw){
      suppressWarnings(RAUW)
    } else {
      suppressWarnings(GG)
    }
  } else {
    
    ####################################################################### \ LFC
    # GG: in
    GG = NULL
    if(nregion == 1){
      GG = ggplot2::ggplot(discovery$smooth, ggplot2::aes(x = log10(distance), y = P ,col = samplename)) +
        ggplot2::labs(x = 'distance (Mb)', col = 'sample') +
        ggplot2::scale_color_manual(values = smplCols)
    } else if(flipFacet){
      GG =ggplot2::ggplot(discovery$smooth, ggplot2::aes(x = log10(distance), y = P ,col = region)) +
        ggplot2::facet_grid(. ~ samplename)+
        ggplot2::labs(x = 'distance (Mb)', col = 'region') +
        ggplot2::scale_color_manual(values = regCols)
    } else {
      GG =ggplot2::ggplot(discovery$smooth, ggplot2::aes(x = log10(distance), y = P ,col = samplename)) +
        ggplot2::facet_grid(. ~ region)+
        ggplot2::labs(x = 'distance (Mb)', col = 'sample')+
        ggplot2::scale_color_manual(values = smplCols)
    }
    RAUW = GG +   
      ggplot2::scale_y_log10() + 
      ggplot2::geom_line() 
    
    # GG: cloud
    if(metric == 'both'){
      GG = GG + ggplot2::geom_point(data = discovery$raw, pch = '.', cex = 1, alpha = 0.01) 
      RAUW = RAUW + ggplot2::geom_point(data = discovery$raw, pch = '.', cex = 1, alpha = 0.01) 
    }
    
    # GG: misc
    breaks = unique(round(log10(unique(discovery$smooth$distance))))
    breaks = breaks[is.finite(breaks)]
    
    GG = GG + ggplot2::theme_classic() +
      ggplot2::scale_x_continuous(breaks = breaks, labels = paste0((10**breaks)/1e6)) + 
      ggplot2::scale_y_log10() + 
      ggplot2::geom_line() +
      ggplot2::coord_fixed(ylim = range(discovery$smooth$P))
    
    GG = GG +  ggplot2::theme(panel.background = ggplot2::element_blank(),
                              aspect.ratio = 1,
                              strip.background = ggplot2::element_rect(fill = NA, colour = NA),
                              panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                              text = ggplot2::element_text(color = 'black'),
                              axis.line = ggplot2::element_blank(),
                              axis.text = ggplot2::element_text(colour = 'black'),
                              strip.text = ggplot2::element_text(colour = 'black') )
    
    if(!is.null(title)){
      GG = GG + ggplot2::ggtitle(title)
    }
    
    if(raw){
      suppressMessages(RAUW)
    } else {
      suppressWarnings(GG)
    }
    
  }
}

#' @rdname visualise
#' @section Extended arguments: \subsection{Virtual 4C}{
#' \describe{
#'  \item{\code{bins}}{An \code{integer} of length 1 with the number of bins to
#'   aggregate signal in. If \code{NULL} (the default), this is set to the 
#'   number of Hi-C bins in the viewpoint's chromosome.}
#'  \item{\code{bedlist}}{Either a BED-formatted \code{data.frame} or a 
#'  \code{list} thereof, indicating genomic intervals to annotate in the bottom 
#'  margin of the plot.}
#'  \item{\code{extend_viewpoint}}{An \code{integer} of length one in basepairs
#'  indicating by how much to widen the viewpoint censor-box. Affects the 
#'  scaling of the y-axis.}
#' }
#' }
#' @export
#' @usage NULL
visualise.virtual4C_discovery <- function(discovery, bins = NULL, 
                                          bedlist = NULL, bed_colours = "black", 
                                          extend_viewpoint = NULL, ...){
  data <- as.data.table(discovery$data)
  VP   <- attr(discovery,"viewpoint")
  data <- data[chromosome == VP[1, 1]]

  expnames <- tail(colnames(data), -2)
  
  if(!is.null(extend_viewpoint)){
    VP[, 2] <- VP[, 2] - extend_viewpoint
    VP[, 3] <- VP[, 3] + extend_viewpoint
  }
  
  if( is.null(bins) ) {
    bins = nrow(data[!duplicated(mid)])
  }
  
  VP_mid <- rowMeans(attr(discovery, 'viewpoint')[, 2:3])
  
  blackout_up   <- VP[, 2]
  blackout_down <- VP[, 3]
  
  data <- melt.data.table(data, id.vars = c("chromosome", "mid"))
  
  data_blackout <- data
  for (i in seq_len(nrow(VP))) {
    if (nrow(VP) > 1) {
      expcheck1 <- data_blackout$variable != expnames[i]
      expcheck2 <- data$variable != expnames[i]
    } else {
      expcheck1 <- expcheck2 <- FALSE
    }
    data_blackout <- data_blackout[(mid >= blackout_up[i] & 
                                      mid <= blackout_down[i]) |
                                     expcheck1]
    data <- data[!(mid >= blackout_up[i] & mid <= blackout_down[i]) |
                   expcheck2]
  }

  if( !is.null(bedlist)) {
    if (inherits(bedlist, "list")) {
      bedlist <- lapply(seq_along(bedlist), function(i) {
        bed <- bedlist[[i]]
        if (!inherits(bed, "data.frame") | is.data.table(bed)) {
          bed <- as.data.frame(bed)
        }
        bed <- bed[bed[,1] == VP[1, 1], 2:3]
        colnames(bed) <- c("start", "end")
        bed$entry <- i
        bed
      })
      bed <- do.call(rbind, bedlist)
    } else {
      if (!inherits(bedlist, "data.frame") | is.data.table(bedlist)) {
        bedlist <- as.data.frame(bedlist)
      }
      bed <- bedlist
      bed <- bed[bed[,1] == VP[1, 1], 2:3]
      colnames(bed) <- c("start", "end")
      bed$entry <- 1
    }
  } else {
    bed <- NULL
  }

  breaks   <- seq(min(data$mid), max(data$mid), length.out = bins)
  bin_size <- median(diff(breaks))
  data <- data[is.finite(value)]
  smooth   <- data[, mean(value), 
                   by = list(findInterval(data$mid, breaks), variable)]
  smooth$mid = breaks[smooth[[1]]] + (bin_size/2)
  smooth[,1] = NULL
  colnames(smooth) = c("variable", "value","mid")
  
  smooth$experiment <- factor(smooth$experiment, levels = expnames)
  data$experiment <- factor(data$experiment, levels = expnames)
  p = ggplot2::ggplot(data, ggplot2::aes(x= mid, y = value)) +
    ggplot2::geom_col(data = smooth, fill = 'black', width = bin_size,
                      colour = NA) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_continuous(name = paste0("Location ", VP[1, 1], " (Mb)"),
                                labels = scales::label_number(scale = 1e-6)) +
    ggplot2::scale_y_continuous(name = "Signal",
                                breaks = function(x) {
                                  scales::extended_breaks()(pmax(x, 0))
                                }) +
    ggplot2::facet_grid(variable ~ .)
  
  # draw_blackout ===+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ymax <- ceiling(max(smooth$value))
  data_blackout$value <- ymax

  if (nrow(data_blackout) > 0) {
    blackout <- data_blackout[, list(min = min(mid), max = max(mid)), 
                              by = "variable"]
    blackout[, "mid" := (min + max) / 2]
    p <- p + ggplot2::geom_rect(
      data = blackout, fill = "#D8D8D8", colour = "#D8D8D8",
      ggplot2::aes(xmin = min, xmax = max, ymin = 0, ymax = ymax),
      inherit.aes = FALSE
    ) +
      ggplot2::geom_text(
        data = blackout, label = "\u2693", vjust = 1,
        ggplot2::aes(mid, ymax * 0.9)
      )
  }
  bed_size <- ymax/50
  bed_padding <- ymax/200
  
  if( !is.null(bed) ){
    bed$colour <- bed_colours[pmin(bed$entry, length(bed_colours))]
    bed$ymin <- -1 * bed_size * bed$entry - bed_padding * (bed$entry)
    bed$ymax <- bed$ymin + bed_size
    n_entries <- length(unique(bed$entry))
    
    p <- p + ggplot2::geom_rect(
      data = bed, inherit.aes = FALSE,
      ggplot2::aes(xmin = start, xmax = end, 
                   ymin = ymin,
                   ymax = ymax),
      fill = rep(bed$colour, length(expnames))
    ) +
      ggplot2::coord_cartesian(
        expand = F,
        ylim = c(-1 * (bed_size * n_entries) - bed_padding * (n_entries + 1),
                 ymax)
      )
  } else {
    p = p + ggplot2::coord_cartesian(expand = F, 
                                     ylim = c(0, ymax))
  }

  p <- p + ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black'),
                          axis.text = ggplot2::element_text(colour = 'black'),
                          strip.background = ggplot2::element_blank())
  suppressWarnings(p)
}

#' @rdname visualise
#' @section Extended arguments:\subsection{Saddle}{
#' \describe{
#'  \item{\code{chr}}{A \code{character} of length one with the chromosome 
#'   name. Defaults to \code{"all"}.}
#'  \item{\code{show_single_contrast}}{A \code{logical} of length 1; 
#'   if \code{FALSE} (default), does not show contrasts when \code{discovery} 
#'   describes one experiment. If \code{TRUE}, plots empty panel.}
#' }
#' }
#' @export
#' @usage NULL
visualise.saddle_discovery <- function(discovery, contrast = 1,
                                       chr = "all",
                                       raw = FALSE, title = NULL,
                                       colour_lim = NULL,
                                       colour_lim_contrast = NULL, 
                                       show_single_contrast = FALSE,
                                       ...) {
  df <- discovery$saddle
  df <- df[!is.na(mean) & !is.na(q1),]
  df$exp <- factor(df$exp, levels = unique(df$exp))
  
  if (chr != "all") {
    if (endsWith(chr, "p") || endsWith(chr, "q")) {
      chromo <- chr
      df <- df[chr == chromo,]
    } else {
      chromo <- df[["chr"]]
      chromo <- substr(chromo, 1, nchar(chromo) - 1)
      keep <- chromo == chr
      df <- df[keep,]
    }
  }
  
  expnames <- df[, unique(exp)]
  
  # Aggregate across chromosomes
  df <- df[, mean(mean), by = list(exp, q1, q2)]
  setnames(df, 4, "obsexp")
  
  # Mirror along diagonal
  comp <- df[q1 != q2, ]
  setnames(comp, 2:3, names(comp)[3:2])
  df <- rbindlist(list(df, comp), use.names = TRUE)

  showcontrast <- !is.null(contrast) && (length(expnames) > 1L || show_single_contrast)
  if (showcontrast) {
    contrast <- df[exp == expnames[contrast],]
    m <- matrix(NA_real_, max(contrast$q1), max(contrast$q2))
    i <- as.matrix(contrast[,list(q1, q2)])
    m[i] <- contrast[["obsexp"]]
    
    contrast <- lapply(expnames, function(j) {
      xx <- df[exp == j]
      mm <- m
      mm[as.matrix(xx[,list(q1, q2)])] <- xx[["obsexp"]]
      mm <- mm - m
      data.table(exp = j, q1 = row(mm)[TRUE], q2 = col(mm)[TRUE], 
                 diff = mm[TRUE])
    })
    contrast <- rbindlist(contrast)
    contrast <- as.data.frame(contrast)
    contrast$panel <- factor("Difference", 
                             levels = c("Individual", "Difference"))
    contrast$q1 <- contrast$q1 / max(contrast$q1)
    contrast$q2 <- contrast$q2 / max(contrast$q2)
  }
  
  df$panel <- factor("Individual",
                     levels = c("Individual", "Difference"))
  df <- as.data.frame(df)
  df$q1 <- df$q1 / max(df$q1)
  df$q2 <- df$q2 / max(df$q2)

  g <- ggplot2::ggplot(df, ggplot2::aes(q1, q2)) +
    ggplot2::geom_raster(hjust = 0, vjust = 0, ggplot2::aes(fill = obsexp))
  
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }
  
  if (showcontrast) {
    
    suppressWarnings(
      g <- g + ggplot2::geom_raster(data = contrast, 
                                    ggplot2::aes(q1, q2, altfill = diff),
                                    inherit.aes = FALSE,
                                    hjust = 0, vjust = 0) +
        ggplot2::facet_grid(panel ~ exp, switch = "y")
    )
    
    if (!raw) {
      g <- g + scale_fill_GENOVA_div(
        palette = "greenpink",
        aesthetics = "altfill",
        oob = scales::squish,
        name = "Difference",
        limits = colour_lim_contrast,
        midpoint = 0
      )
      g$scales$scales[[1]]$guide <- ggplot2::guide_colourbar()
      g$scales$scales[[1]]$guide$available_aes[[3]] <- "altfill"
      g$scales$scales[[1]]$guide$order <- 2
    } else {
      g <- g + ggplot2::guides(
        "altfill" = ggplot2::guide_colourbar(available_aes = "altfill")
      )
    }
    
    # ggnomics style hack
    old_geom <- g$layers[[2]]$geom
    old_nahandle <- old_geom$handle_na
    new_nahandle <- function(self, data, params) {
      colnames(data)[colnames(data) %in% "altfill"] <- "fill"
      old_nahandle(data, params)
    }
    new_geom <- ggplot2::ggproto(paste0(sample(1e6, 1), class(old_geom)),
                                 old_geom,
                                 handle_na = new_nahandle)
    names(new_geom$default_aes)[1] <- "altfill"
    new_geom$non_missing_aes <- "altfill"
    g$layers[[2]]$geom <- new_geom
  } else if (length(unique(df$exp)) > 1) {
    g <- g + ggplot2::facet_grid(~ exp)
  }
  
  if (raw) {
    return(g)
  } else {
    g <- g + ggplot2::theme(
      aspect.ratio = 1,
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = "black",
                                           size = 0.25),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black", size = 0.25),
      strip.text = ggplot2::element_text(colour = "black"),
      panel.spacing.x = grid::unit(0.8 * 0.5, "strwidth", data = c("1.00.0")),
      panel.spacing.y = grid::unit(0.8 * 1.5, "strheight", data = c("1.00.0"))
    ) +
      ggplot2::scale_x_continuous(name = "Quantile", expand = c(0,0),
                                  breaks = c(0, 1), labels = c("B", "A")) +
      ggplot2::scale_y_continuous(name = "Quantile", expand = c(0,0),
                                  breaks = c(0, 1), labels = c("B", "A")) +
      scale_fill_GENOVA_div(
        palette = "divergent",
        midpoint = 1,
        limits = colour_lim,
        oob = scales::squish,
        name = expression(frac("Observed", "Expected")),
        guide = ggplot2::guide_colourbar(order = 1)) +
      ggplot2::coord_cartesian(clip = "off")
  }
  return(g)
}

#' @rdname visualise
#' @usage NULL
#' @export
visualise.domainogram_discovery <- function(discovery, 
                                            colour_lim = c(-1, 1),
                                            title = NULL,
                                            raw = FALSE, ...) {
  df <- discovery$scores
  df <- as.data.table(df)
  df <- melt.data.table(df, id.vars = c("window", "position"), 
                        value.name = "insulation")
  setnames(df, 3, "experiment")
  
  g <- ggplot2::ggplot(df, ggplot2::aes(position, window, fill = insulation)) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(experiment ~ .)
  
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }
  
  if (!raw) {
    g <- g + 
      ggplot2::scale_x_continuous(name = paste0("Position ", 
                                                attr(df, "chrom"), " (Mb)"),
                                  labels = scales::label_number(scale = 1e-6),
                                  expand = c(0,0)) +
      ggplot2::scale_y_continuous(name = "Window Size", expand = c(0, 0)) +
      scale_fill_GENOVA_div(midpoint = 0,
                            trans = "reverse",
                            limits = rev(colour_lim), 
                            oob = scales::squish,
                            name = "Insulation\nScore") +
      ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
                     strip.background = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_line(colour = "black"),
                     axis.line  = ggplot2::element_line(colour = "black"),
                     panel.grid = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
  }
  
  g
}

#' @rdname visualise
#' @section Extended arguments:\subsection{Intra-inter TAD}{
#' \describe{
#'  \item{\code{geom}}{One of the following length one \code{character}s:
#'    \describe{
#'      \item{\code{"boxplot"}}{Display as boxplots.}
#'      \item{\code{"violin"}}{Display as violin plots.}
#'      \item{\code{"jitter"}}{Display as jittered points plot.}
#'    }
#'  }
#'  \item{\code{censor_contrast}}{A \code{logical} of length 1 deciding whether 
#'  the contrasting experiment itself should be censored (\code{TRUE}) or
#'  included (\code{FALSE}).}
#' }
#' }
#' @export
#' @usage NULL
visualise.IIT_discovery <- function(discovery, contrast = 1, raw = FALSE,
                                    geom = c("boxplot", "violin", "jitter"),
                                    censor_contrast = TRUE, title = NULL, 
                                    show_single_contrast = FALSE,
                                    ...) {
  geom <- match.arg(geom)
  dat <- as.data.table(discovery$results)
  cols <- attr(discovery, "colours")
  
  expnames <- tail(colnames(dat), -2)
  
  if (!is.null(contrast) & length(expnames) < 2) {
    if (show_single_contrast) {
      message("Cannot compute a contrast for one sample. Reverting to ",
              "visualising plain values.")
    }
    contrast <- NULL
  } else if (!is.null(contrast)){
    contrast <- expnames[contrast]
    
    trans <- lapply(setNames(expnames, expnames), function(i) {
      log2(dat[[i]] / dat[[contrast]])
    })
    
    for (i in expnames) {
      dat[, as.character(i) := trans[[i]]]
    }
    
    if (censor_contrast & !is.null(contrast)) {
      dat <- dat[, -..contrast]
      cols <- cols[which(expnames != contrast)]
    }
  }
  
  df <- melt.data.table(
    dat, id.vars = c("x", "y"), 
    measure.vars = intersect(colnames(dat), expnames)
  )
  
  df$diff <- as.factor(df$y - df$x)

  g <- ggplot2::ggplot(df, ggplot2::aes(diff, value))
  
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }
  
  if (geom == "boxplot") {
    g <- g + ggplot2::geom_boxplot(ggplot2::aes(fill = variable),
                                   key_glyph = "polygon")
  } else if (geom == "violin") {
    g <- g + ggplot2::geom_violin(ggplot2::aes(fill = variable))
  } else if (geom == "jitter") {
    g <- g + ggplot2::geom_point(ggplot2::aes(colour = variable),
                                 position = ggplot2::position_jitterdodge(),
                                 size = 1, alpha = 0.3, shape = 16)
  }
  
  if (raw) {
    return(g)
  }
  
  if ("fill" %in% names(g$layers[[1]]$mapping)) {
    g <- g + ggplot2::scale_fill_manual(
      values = cols, name = "Experiment"
    )
  } else {
    g <- g + ggplot2::scale_colour_manual(
      values = cols, name = "Experiment",
      guide = ggplot2::guide_legend(
        override.aes = list(size = 2, alpha = 1)
      )
    )
  }
  
  g <- g + ggplot2::scale_x_discrete(
    name = "TAD Distance",
    labels = scales::math_format(n + .x)
    # labels = scales::label_number(prefix = "n +")
  )
  
  if (is.null(contrast)) {
    g <- g + ggplot2::scale_y_continuous(
      name = "Contacts",
      trans = "log10",
      labels = scales::math_format(format = log10)
    )
  } else {
    ytitle <- bquote("Log"[2]~" (Experiment contacts /" ~ 
                       paste(.(contrast)) ~ "contacts)")
    g <- g + ggplot2::scale_y_continuous(
      name = ytitle,
    )
  }
  
  g <- g + ggplot2::theme(
    text = ggplot2::element_text(colour = "black"),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.background = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  )
  
  return(g)
}

#' @rdname visualise
#' @usage NULL
#' @export
visualise.chrommat_discovery <- function(discovery, raw = FALSE, title = NULL,
                                         colour_lim = NULL, ...) {
  obsexp <- discovery$obs
  dim <- dim(obsexp)
  if (attr(discovery, "mode") %in% c("trans", "regress")) {
    tmp <- obsexp
    tmp[slice.index(obsexp, 1) == slice.index(obsexp, 2)] <- 0
  } else {
    tmp <- obsexp
  }
  sums <- apply(tmp, 3, sum)
  obsexp[] <- log2((obsexp / rep(sums, each = prod(dim[1:2]))) / discovery$exp)
  
  df <- data.frame(
    row = as.vector(slice.index(obsexp, 1)),
    col = as.vector(slice.index(obsexp, 2)),
    exp = as.vector(slice.index(obsexp, 3))
  )
  df[] <- mapply(function(x, i){x[i]}, x = dimnames(obsexp), i = df)
  df$value <- as.vector(obsexp)
  
  if (is.null(colour_lim)) {
    colour_lim <- range(with(df, value[row != col]))
    colour_lim[1] <- min(colour_lim[1], -1)
    colour_lim[2] <- max(colour_lim[2], 1) 
  }

  
  if (utils::packageVersion("ggplot2") > "3.2.1") {
    guide_x <- ggplot2::guide_axis(check.overlap = TRUE, angle = 90)
    guide_y <- ggplot2::guide_axis(check.overlap = TRUE)
  } else {
    guide_x <- guide_y <- ggplot2::waiver()
  }
  
  g <- ggplot2::ggplot(df, ggplot2::aes(row, col, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::facet_grid(~ exp) +
    ggplot2::scale_x_discrete(limits = dimnames(obsexp)[[1]], guide = guide_x,
                              expand = c(0,0), name = "") +
    ggplot2::scale_y_discrete(limits = dimnames(obsexp)[[2]], guide = guide_y,
                              expand = c(0,0), name = "")
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }
  if (raw) {
    return(g)
  }
  
  g <- g + scale_fill_GENOVA_div(
    name = expression(Log[2]*frac("Observed", "Expected")),
    limits = colour_lim, oob = scales::squish,
    midpoint = 0
  ) +
    ggplot2::coord_equal() +
    GENOVA_THEME()
  g
}

# Utilities ---------------------------------------------------------------

#' scale_altfill_continuous Makes sure no errors are returned when
#' visualise(..., raw = TRUE) Unfortunately has to be exported, but that is
#' better than throwing errors anytime an attempt is mode to print a raw plot.
#'
#' @param ... Arguments to \code{\link[ggplot2]{scale_fill_gradient2}}.
#'
#' @export
scale_altfill_continuous <- function(...) {
  guide <- ggplot2::guide_colourbar()
  guide$available_aes <- "altfill"
  ggplot2::scale_fill_gradient2(..., guide = guide,
                                aesthetics = "altfill")
}

# Function factory for centering limits around a value
centered_limits <- function(around = 0) {
  function(input) {
    c(-1, 1) * max(abs(input - around)) + around
  }
}

symmetric_quantiles <- function(x, q = c(0.05, 0.95)) {
  c(-1, 1) * diff(quantile(x[is.finite(x)], q))
}

# Function factory labelling
label_kilobase_relative <- function(zero = "0") {
  fun <- scales::label_number(scale = 1e-3, suffix = " kb")
  function(x) {
    ifelse(x == 0, zero, fun(x))
  }
}

# Break function factory
breaks_trim_outer <- function(min = 3) {
  function(x) {
    x <- scales::extended_breaks()(x)
    if (length(x) > min) head(tail(x, -1), -1) else x
  }
}

# Convenience function for filtering layer data when inheriting from
# main ggplot2 call
datafilter <- function(expr = NULL) {
  expr <- substitute(expr)
  if (is.null(expr)) {
    expr <- substitute(TRUE)
  }
  f <- function(x) subset.data.frame(x, eval(expr))
  parent.env(environment(f)) <- parent.frame()
  f
}
