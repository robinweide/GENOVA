# User function -----------------------------------------------------------

#' @rdname pyramid
#' @title Plot a region as a pyramid
#'
#' @description Plots a region of the Hi-C map with a 45 degree rotation. The
#'   rotation makes it such that the diagonal coincides with the x-axis and the
#'   y-axis indicates distance. \code{pyramid_difference} does the same, but
#'   first subtracts one sample from the other.
#'
#' @param exp A GENOVA \code{contacts} or \code{contacts_matrix} object.
#' @param exp1,exp2 Like \code{exp} but only used in
#'   \code{pyramid_difference()}.
#' @param chrom A \code{character} of length one indicating a chromosome name.
#' @param start A \code{numeric} of length one for the start position in
#'   basepairs.
#' @param end A \code{numeric} of length one for the end position in basepairs.
#' @param crop_x A \code{numeric} of length two indicating positions in
#'   basepairs of where to crop the x-axis.
#' @param crop_y A \code{numeric} of length two indicating distances in
#'   basepairs of where to crop the y-axis.
#' @param colour_scale A continuous ggplot2 scale with the \code{"fill"}
#'   aesthetic.
#' @param display_yaxis A \code{logical} of length 1: should the y-axis be
#'   displayed?
#' @param raw A \code{logical} of length 1: should a bare bones plot be
#'   returned? When \code{TRUE}, no position- or colour-scales or theme elements
#'   are added.
#' @param ... Arguments passed to downstream functions.
#'
#' @details Some \code{colour_scale} settings are adjusted. If no \code{limits}
#'   were set, new limits are set as \code{c(0, quantile(x, 0.975))}, wherein
#'   \code{x} are the Hi-C values. Also, the \code{oob} parameter is replaced by
#'   \code{scales::squish()} and the name is set by default to
#'   \code{"Contacts"}.
#'
#' @return A \code{ggplot} object
#' @export
#'
#' @section Annotations: For annotations along the linear genome, see
#'   \code{\link[GENOVA]{pyramidtracks}}. For annotations in the Hi-C map, see
#'   \code{\link[GENOVA]{add_tads}} and \code{\link[GENOVA]{add_loops}}.
#'
#' @examples
#' \dontrun{
#' pyramid(exp, "chr2", 25e6, 30e6)
#' }
pyramid <- function(exp, chrom, start, end, crop_x, crop_y, 
                    colour_scale, display_yaxis, raw, ...) {
  UseMethod("pyramid", exp)
}

#' @method pyramid default
#' @export
pyramid.default <- function(exp, ...) {
  stop("Don't know how to make a pyramid plot out of this data.")
}

#' @method pyramid contacts
#' @export
pyramid.contacts <- function(exp, chrom = "chr1", start = 0, end = 25e6, ...) {
  y <- select_subset(exp, chrom, start, end)
  if (attr(exp, "znorm")) {
    pyramid(y, colour_scale = scale_fill_GENOVA_div(
      name = "Z-score",
      midpoint = 0, limits = c(NA, NA))
    )
  } else {
    pyramid(y, ...)
  }
}

#' @method pyramid contacts_matrix
#' @export
pyramid.contacts_matrix <- function(
  exp,
  ...
) {
  rowcoord <- exp$x
  colcoord <- exp$y
  res <- attr(exp, "resolution")
  location <- range(rowcoord) + c(-0.5, 0.5) * res
  location <- list(attr(exp, "chrom"), location[1], location[2])
  
  .core_pyramid(
    x = rowcoord[as.vector(row(exp$z))],
    y = colcoord[as.vector(col(exp$z))],
    z = as.vector(exp$z),
    location = location,
    resolution = res,
    ...
  )
}

#' @method pyramid matrix
#' @export
pyramid.matrix <- function(exp, ...) {
  .core_pyramid(
    x = as.vector(row(exp)),
    y = as.vector(col(exp)),
    z = as.vector(exp),
    location = list(NULL, 0.5, dim(exp)[1] + 0.5),
    resolution = 1,
    ...
  )
}

#' @export
#' @rdname pyramid
pyramid_difference <- function(exp1, exp2, chrom, 
                               start, end, ...) {
  a <- select_subset(exp1, chrom, start, end)
  b <- select_subset(exp2, chrom, start, end)
  a$z <- a$z - b$z
  
  pyramid(a, colour_scale = scale_fill_GENOVA_div(
    name = "Difference",
    midpoint = 0,
    limits = quantile(a$z, c(0.05, 0.95)))
  )
}

.core_pyramid <- function(x, y, z, 
                         crop_x = c(-Inf, Inf), 
                         crop_y = c(-Inf, Inf), 
                         location,
                         resolution,
                         colour_scale = scale_fill_GENOVA(),
                         display_yaxis = FALSE,
                         raw = FALSE) {
  # Build data
  df <- data.table(
    x = x, y = y, contacts = z
  )
  
  # Polygonise rectangles
  xres <- c(-1, -1, 1, 1) * resolution * 0.5
  yres <- c(-1, 1, 1, -1) * resolution * 0.5
  df <- df[, list(x = x + xres, 
                  y = y + yres, 
                  contacts = rep(contacts, 4)),
           by = list(id = seq_len(NROW(df)))]
  # Do 45 degree rotation
  df[, c("x", "y") := .transform_xy_coords(x, y)]
  # Trim sawtooth
  df <- df[y > -1]
  
  # Apply cropping
  if (!is.null(crop_x) && length(crop_x) == 2L) {
    crop_x <- sort(crop_x)
    location[[2]] <- max(crop_x[[1]], location[[2]])
    location[[3]] <- min(crop_x[[2]], location[[3]])
    df <- df[x >= crop_x[1] & x <= crop_x[2]]
  }
  if (!is.null(crop_y) && length(crop_y) == 2L) {
    crop_y <- sort(crop_y)
    df <- df[y >= crop_y[1] & y <= crop_y[2]]
    crop_y <- ifelse(is.finite(crop_y), crop_y, range(df$y))
  } else {
    crop_y <- range(df$y)
  }
  
  df$facet <- factor(".pyramid", levels = ".pyramid")
  asp_ratio <- diff(range(df$x)) / diff(range(df$y))
  
  # Make layer
  triangle <- ggplot2::geom_polygon(
    data = as.data.frame(df),
    ggplot2::aes(x, y, group = id, fill = contacts)
  )
  
  # Setup main plot
  p <- structure(
    list(
      data = ggplot2::waiver(),
      layers = list(triangle),
      scales = ggplot2::ggplot()$scales,
      mapping = ggplot2::aes(),
      theme = list(),
      coordinates = ggplot2::coord_cartesian(default = TRUE),
      facet = facet_pyramid(ggplot2::vars(facet), col_size = 2 * asp_ratio),
      plot_env = parent.frame(),
      location = location,
      labels = list(x = location[[1]], y = "distance", fill = "contacts")
    ), class = c("ggpyramid", "gg", "ggplot")
  )
  
  # Tweak plot
  if (!raw) {
    # Setup colour scale
    # Always replace `oob`
    colour_scale$oob <- scales::squish
    # Only set limits when none are given
    if (is.null(colour_scale$limits)) {
      colour_scale$limits <- c(0, quantile(df$contacts, 0.975))
    }
    # Only set name when none is given
    if (inherits(colour_scale$name, "waiver")) {
      if (!is.null(location[[1]])) {
        colour_scale$name <- "Contacts"
      } else {
        colour_scale$name <- "Value"
      }
    }
    
    # Add y-scale
    if (display_yaxis) {
      yscale <- ggplot2::scale_y_continuous(
        name = "Distance", expand = c(0, 0),
        labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
        limits = crop_y
      )
    } else {
      yscale <- ggplot2::scale_y_continuous(
        # breaks = NULL,
        guide = ggplot2::guide_none(),
        name = "", expand = c(0, 0),
        limits = crop_y
      )
    }
    
    p$facet$add_scale(yscale, ".pyramid")
    
    p <- p + GENOVA_THEME() +
      ggplot2::theme(panel.border = ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_line(),
                     strip.placement = "outside",
                     strip.text.y.left = ggplot2::element_text(angle = 0))
    if (!is.null(location[[1]])) {
      p <- p + ggplot2::scale_x_continuous(
        expand = c(0, 0),
        limits = c(location[[2]], location[[3]]),
        labels = scales::label_number(1, scale = 1e-6, suffix = " Mb")
      )
    } else {
      p <- p + ggplot2::scale_x_continuous(
        expand = c(0,0),
        limits = c(location[[2]], location[[3]]),
        breaks = scales::breaks_extended(5, Q = c(1,5,2,4,3))
      )
    }
    p <- p + colour_scale
  }
  ggplot2::set_last_plot(p)
  return(p)
}

# Facets ---------------------------------------------------------------

# In order to flexibly add new annotation layers, we need a facet system that
# allows this.

# This is an unwieldy facet function that protests against anything you would
# ideally like to do. I don't recommend using it interactively.
facet_pyramid <- 
  function(rows = NULL, scales = "free_y", space = "free_y",
           shrink = TRUE, labeller = ggplot2::label_value, as.table = TRUE, 
           switch = "y", drop = TRUE, margins = FALSE,
           col_size = 2, row_size = 1
  ) {
    scales <- match.arg(scales, c("fixed", "free_x", "free_y", "free"))
    space  <- match.arg(scales, c("fixed", "free_x", "free_y", "free"))
    free <- list(x = any(scales %in% c("free_x", "free")),
                 y = any(scales %in% c("free_y", "free")))
    space <- list(x = any(space %in% c("free_x", "free")),
                  y = any(space %in% c("free_x", "free")))
    if (!is.null(switch)) {
      switch <- match.arg(switch, c("x", "y", "both"))
    }
    
    # Don't print the pyramid facet label
    labfun <- function(x){ifelse(x == ".pyramid", "", as.character(x))}
    labfun <- ggplot2::as_labeller(labfun)
    
    names(rows) <- "facets"

    ggplot2::ggproto(
      NULL, FacetPyramid, shrink = shrink, 
      ratios = list(rows = row_size, cols = col_size),
      params = list(rows = rows, cols = ggplot2::vars(),
                    margins = margins, free = free, space_free = space,
                    labeller = labfun, as.table = as.table, switch = switch,
                    drop = drop)
    )
  }

# Basically this is a combination of ggh4x::facetted_pos_scales() and 
# ggh4x::force_panelsizes() with a few extras. These extra's are:
#  * Adding scales dynamically instead of statically.
#  * A clone function to deal with the fact that ggproto's are based on
#    environments, so we create a child that inherits from the parent 
#    environment. Necessary if plots are forked.
# Should work fine except I don't know how to allow transformative position 
# scales (e.g. scale_(x/y)_log10()) to transform data prior to 
# calculating stats. Hence, we should pre-transform stat layers at the data or 
# aes() level.
FacetPyramid <- ggplot2::ggproto(
  "FacetPyramid", ggplot2::FacetGrid,
  add_scale = function(scale, facet_name, self) {
    current <- self$y_scales
    if (facet_name %in% names(current)) {
      current[[facet_name]] <- scale
    } else {
      len <- length(current)
      current <- c(current, scale)
      names(current)[len + 1] <- facet_name
    }
    self$y_scales <- current
    return(invisible())
  },
  y_scales = list(),
  ratios = list(rows = 1, cols = 2),
  init_scales = function(layout, x_scale = NULL, y_scale = NULL, params, self) {
    # As in ggh4x::facetted_pos_scales()
    # We initialise a new y-scale for every facet
    scales <- list()
    # Initialise x-scales as per usual
    if (!is.null(x_scale)) {
      scales$x <- lapply(seq_len(max(layout$SCALE_X)), 
                         function(i) x_scale$clone())
    }
    if (!is.null(y_scale)) {
      new_y <- self$y_scales
      yidx <- seq_len(max(layout$SCALE_Y))
      
      # Check if there are any additional y-scales
      if (length(new_y) == 0) {
        scales$y <- lapply(yidx, function(i) {
          y_scale$clone()
        })
      } else {
        scales$y <- lapply(yidx, function(i) {
          if (length(new_y) >= i) {
            if (!is.null(new_y[[i]])) {
              new <- new_y[[i]]$clone()
              new$oob <- function(x, ...) x
              return(new)
            }
          }
          return(y_scale$clone())
        })
      }
    }
    return(scales)
  },
  train_scales = function(x_scales, y_scales, layout, data, params, self) {
    # As in ggh4x::facetted_pos_scales()
    # Scale limits are trained on transformed data
    data <- lapply(data, function(layer_data) {
      self$finish_data(layer_data, layout,
                       x_scales, y_scales, params)
    })
    ggplot2::ggproto_parent(ggplot2::Facet, self)$train_scales(
      x_scales, y_scales, layout, data, params
      )
  },
  finish_data = function(data, layout, x_scales, y_scales, params) {
    # As in ggh4x::facetted_pos_scales()
    # Performs the data transformation. Note that this step occurs after (!)
    # stat calculation, so stat layers need to be pre-transformed.
    panels <- split(data, data$PANEL, drop = FALSE)
    panels <- lapply(names(panels), function(i) {
      dat <- panels[[i]]
      panel_id <- match(as.numeric(i), layout$PANEL)
      xidx <- layout[panel_id, "SCALE_X"]
      yidx <- layout[panel_id, "SCALE_Y"]
      
      y_vars <- intersect(y_scales[[yidx]]$aesthetics, names(dat))
      x_vars <- intersect(x_scales[[xidx]]$aesthetics, names(dat))
      
      for (j in y_vars) {
        dat[, j] <- y_scales[[yidx]]$transform(dat[, j])
      }
      for(j in x_vars) {
        dat[, j] <- x_scales[[xidx]]$transform(dat[, j])
      }
      dat
    })
    data <- unsplit(panels, data$PANEL)
    data
  },
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                          data, theme, params, self) {
    # As in ggh4x::force_panelsizes()
    # Basically, execute FacetGrid drawing, then force sizes
    gt <- ggplot2::ggproto_parent(ggplot2::FacetGrid, self)$draw_panels(
      panels, layout, x_scales, y_scales, ranges, coord, data, theme, params
    )
    
    prows <- ggplot2::panel_rows(gt)
    pcols <- ggplot2::panel_cols(gt)
    
    ratios <- self$ratios
    xsize <- ratios$cols
    ysize <- ratios$rows
    
    if (length(ysize) > 0) {
      rowheights <- rep(ysize, length.out = nrow(prows))
      gt$heights[prows$t] <- ggplot2::unit(rowheights, "null")
    }
    if (length(xsize) > 0) {
      colwidths <- rep(xsize, length.out = nrow(pcols))
      gt$widths[pcols$l] <- ggplot2::unit(colwidths, "null")
    }
    
    gt$respect <- TRUE
    gt
  },
  clone = function(self) {
    # Whenever we add a new facet layer we should clone the Facet ggproto
    # such that forked plot construction still works.
    ggplot2::ggproto(NULL, self)
  }
)

# Adding tracks -----------------------------------------------------------

#' @name pyramidtracks
#' @title Adding tracks to a pyramid plot.
#'
#' @description These functions add new facets to a pyramid plot with additional
#'   annotations.
#'
#' @param layer A ggplot layer instance
#' @param size A \code{numeric} of length one noting the relative size to the
#'   main facet panel.
#' @param name A \code{character} of length one title to display as facet strip
#'   label.
#' @param scale_y A ggplot position scale with the y-aesthetic to control
#'   parameters of the track's y-scale.
#'
#' @details Currently, \code{as_track()} only accepts layers wherein the
#'   \code{data} argument is explicitly defined before evaluation. For example,
#'   \code{geom_point(aes(x = 1, y = 1))} is not an appropriate track, while
#'   \code{geom_point(aes(x + 1, y = y), data = data.frame(x = 1, y = 1))} is an
#'   appropriate track.
#'
#' @return A \code{track_layer} object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialising a pyramid
#' p <- pyramid(exp, "chr2", 25e6, 30e6)
#'
#' # Adding a bed annotation
#' p + add_bed_track(bed, fill = "blue", name = "bed")
#'
#' # Adding a bedgraph annotation
#' p <- p + add_bed_graph(bedgraph, linetype = 2, name = "Some statistic")
#' p
#'
#' # Adding custom annotation
#' # Use the name to display in the same panel as previous annotation
#' p + as_track(ggplot2::geom_vline(xintercept = 27.5e6),
#'              name = "Some statistic")
#' }
as_track <- function(layer, 
                     size = 0.1, 
                     name,
                     scale_y = NULL) {
  if (is.null(layer$data$facet)) {
    layer$data$facet <- factor(name)
  }
  name <- force(name)
  if (is.null(scale_y)) {
    scale_y <- ggplot2::scale_y_continuous(n.breaks = 3)
  }
  structure(
    list(
      layer = layer,
      scale = scale_y,
      size = size
    ),
    class = "track_layer"
  )
}

#' @export
#' @rdname pyramidtracks
#' @param ... Additional arguments to pass to \code{ggplot2::layer()}.
#' @inheritParams pyramidannotations
add_bed_track <- function(bed, size = 0.1, name, ...) {
  name <- force(name)
  df <- as.data.frame(bed[,1:3])
  names(df)[1:3] <- c("chrom", "xmin", "xmax")
  df$facet <- factor(name)
  
  rct <- ggplot2::geom_rect(ggplot2::aes(xmin = xmin, ymin = 0,
                                         xmax = xmax, ymax = 1),
                            data = df, ...)
  scale <- ggplot2::scale_y_continuous(breaks = NULL)
  as_track(rct, size = size, name = name, scale_y = scale)
}

#' @export
#' @rdname pyramidtracks
#' @param bedgraph A \code{data.frame} with four columns noting the chromosome
#'   name, start- and end-positions of genomic intervals, and a value to
#'   display.
add_bed_graph <- function(bedgraph, size = 0.1, name, ...) {
  name <- force(name)
  df <- data.table(bedgraph[,1:4])
  setnames(df, 1:4, c("chrom", "start", "end", "y"))
  df <- df[, list(chrom, x = (start + end)/2, y)]
  setDF(df)
  df$facet <- factor(name)
  pth <- ggplot2::geom_path(ggplot2::aes(x, y), data = df, ...)
  as_track(pth, size = size, name = name)
}

#' @export
#' @rdname pyramidtracks
#' @usage NULL
#' @keywords internal
ggplot_add.track_layer <- function(object, plot, object_name) {
  lay <- object$layer
  if (inherits(plot, "ggpyramid")) {
    layerdata <- lay$data
    loca <- plot$location
    if ("chrom" %in% names(layerdata)) {
      layerdata <- layerdata[layerdata$chrom == loca[[1]],]
    }
    # variable names from ggplot2:::ggplot_global$x_aes
    xvar <- intersect(names(layerdata), 
                      c("x", "xmin", "xmax", "xend", "xintercept",
                        "xmin_final", "xmax_final", "xlower", "xmiddle",
                        "xupper", "x0"))
    validpos <- lapply(layerdata[xvar], function(x) {
      x >= loca[[2]] & x <= loca[[3]]
    })
    validpos <- rowSums(do.call(cbind, validpos)) >= 1
    lay$data <- layerdata[validpos,]
  }
  if (NROW(lay$data) == 0) {
    warning("No appropriate data was found for this position.", call. = FALSE)
  }
  plot <- ggplot2::ggplot_add(lay, plot = plot, 
                              object_name = object_name)
  if (inherits(plot$facet, 'FacetPyramid')) {
    plot$facet <- plot$facet$clone()
    plot$facet$add_scale(object$scale, object$layer$data$facet[[1]])
    plot$facet$ratios$rows <- c(plot$facet$ratios$rows, object$size)
    plot <- .set_facet_factorlevel(plot)
  }
  plot
}

# Adding marks ------------------------------------------------------------

# Main difference with tracks is that info would need to be rotated

#' @name pyramidannotations
#' @aliases add_tads add_loops
#' @title Rotated pyramid annotations
#'
#' @description Add annotations to a pyramid plot that also require the 45
#'   degree rotation.
#'
#' @param bed A BED-formatted \code{data.frame} with the following 3 columns:
#'   \enumerate{ \item A \code{character} giving the chromosome names. \item An
#'   \code{integer} with start positions. \item An \code{integer} with end
#'   positions. }
#' @param bedpe A BEDPE-formatted \code{data.frame} with the following 6
#'   columns: \enumerate{ \item A \code{character} giving the chromosome names
#'   of the first coordinate. \item An \code{integer} giving the start positions
#'   of the first coordinate. \item An \code{integer} giving the end positions
#'   of the first coordinate. \item A \code{character} giving the chromosome
#'   names of the second coordinate. \item An \code{integer} giving the start
#'   positions of the second coordinate. \item An \code{integer} giving the end
#'   positions of the second coordinate. }
#' @param colour,shape,... Additional parameters to pass to
#'   \code{ggplot2::layer()}.
#'
#' @return A \code{mark_layer} object.
#'
#' @examples
#' \dontrun{
#' p <- pyramid(exp, "chr2", 25e6, 30e6)
#' p + add_tads(tads)
#' p + add_loops(loops)
#' }
NULL

#' @export
#' @rdname pyramidannotations
add_tads <- function(bed, colour = "limegreen", ...) {
  df <- as.data.frame(bed[, 1:3])
  n <- NROW(df)
  seq <- rep(seq_len(n), each = 2)
  df <- .list_as_df(list(
    chrom = df[[1]],
    x = c(df[[2]][[1]], df[[2]][tail(seq, -1)], df[[3]][[n]]),
    y = c(df[[2]][[1]], df[[3]][head(seq, -1)], df[[3]][[n]])
  ))
  layer <- ggplot2::geom_path(data = df, 
                              ggplot2::aes(x = x, y = y), 
                              colour = colour, ...)
  
  structure(list(layer = layer), class = "mark_layer")
}

#' @export
#' @rdname pyramidannotations
add_loops <- function(bedpe, colour = "limegreen", shape = 1, ...) {
  bedpe <- bedpe[bedpe[[1]] == bedpe[[4]],] # drop trans
  x = (bedpe[[2]] + bedpe[[3]]) / 2
  y = (bedpe[[5]] + bedpe[[6]]) / 2
  df <- data.frame(
    chrom = bedpe[[1]],
    x = pmin(x, y),
    y = pmax(x, y)
  )
  layer <- ggplot2::geom_point(data = df,
                               ggplot2::aes(x = x, y = y),
                               colour = colour, shape = shape, ...)
  structure(list(layer = layer), class = "mark_layer")
}

#' @export
#' @rdname pyramidannotations
#' @usage NULL
#' @keywords internal
ggplot_add.mark_layer <- function(object, plot, object_name) {
  lay <- object$layer
  
  if (inherits(plot, "ggpyramid")) {
    data <- lay$data
    loca <- plot[["location"]]
    data <- subset(data, chrom == loca[[1]])
    if (inherits(lay$geom, "GeomPoint")) {
      validpos <- data$x >= loca[[2]] & data$x <= loca[[3]] &
        data$y >= loca[[2]] & data$y <= loca[[3]]
      data <- data[validpos, ]
    } else {
      validpos <- with(data, (x > loca[[2]] & x < loca[[3]]) | 
                         (y >= loca[[2]] & y <= loca[[3]]))
      data <- data[validpos, ]
      data <- .list_as_df(.clip_poly(data$x, data$y,
                                     xrange = c(loca[[2]], loca[[3]]),
                                     yrange = c(loca[[2]], loca[[3]])))
    }
    data$facet <- plot$layers[[1]]$data$facet[1]
    data[c("x", "y")] <- .transform_xy_coords(data$x, data$y)
    lay$data <- data
  }
  plot <- ggplot2::ggplot_add(lay, plot = plot, 
                              object_name = object_name)
  plot
}


# Adding discoveries ------------------------------------------------------

# Not satisfied nor finished with this yet

discovery_as_layer <- function(discovery, location, ...) {
  UseMethod("discovery_as_layer")
}

discovery_as_layer.genomescore_discovery <- function(discovery, location, ...) {
  df <- .main_data(discovery)
  expnames <- expnames(discovery)
  setDT(df)
  df <- df[chrom == location[[1]] & start >= location[[2]] &
             end <= location[[3]]]
  
  df <- .list_as_df(list(
    mid = (df[["start"]] + df[["end"]])/2,df[, ..expnames])
  )
  df <-  .list_as_df(
    list(x = rep(df$mid, length(expnames)),
         y = unlist(df[2:ncol(df)]),
         colour = factor(rep(expnames, each = nrow(df)), 
                         levels = expnames),
         facet = factor("Compartment\nScore"))
  )
  layer <- ggplot2::geom_line(ggplot2::aes(x, y, colour = colour), data = df, ...)
  as_track(layer, name = gsub(" ", "\n", .metric_name(discovery)))
}

ggplot_add.genomescore_discovery <- function(object, plot, object_name) {
  if (inherits(plot, "ggpyramid")) {
    plot + discovery_as_layer(object, plot$location)
  } else {
    NextMethod()
  }
}

# Helpers -----------------------------------------------------------------

.set_facet_factorlevel <- function(plot) {
  all_levels <- lapply(plot$layers, function(layer) {
    levels(layer$data$facet)
  })
  all_levels <- Reduce(union, all_levels)
  plot$layers <- lapply(plot$layers, function(layer) {
    layer$data$facet <- factor(as.character(layer$data$facet), all_levels)
    layer
  })
  plot
}

# Angle is counterclockwise and in degrees
# Scale/share is for x- and y-direction respectively
# Order of operations is scale, shear, rotate
.transform_xy_coords <- function(x, y, 
                                 scale = sqrt(2) * c(0.5, 1),
                                 shear = c(0, 0),
                                 angle = 45) {
  # Quick shortcut for pyramid case
  if (angle == 45 && all(scale == sqrt(2) * c(0.5, 1)) && 
      all(shear == 0)) {
    
    rotmat <- matrix(c(0.5, -1, 0.5, 1), ncol = 2)
    
  } else {
    
    rotmat <- diag(2)
    # Apply scale
    rotmat <- rotmat * scale
    # Apply shear
    rotmat <- rotmat %*% matrix(c(1, shear, 1), ncol = 2)
    # Apply rotation
    angle <- -angle * pi / 180
    angle <- matrix(
      c(cos(angle),  sin(angle),
        -sin(angle), cos(angle)),
      ncol = 2
    )
    rotmat <- rotmat %*% angle
    
  }
  
  coords <- matrix(c(x, y), ncol = 2)
  # Transform
  coords <- t(rotmat %*% t(coords))
  .list_as_df(list(x = coords[, 1], y = coords[, 2]))
}

# From R4.0.0, cheapo data.frame constructor that doesn't care about rownames,
# column classes coersion or any of that stuff
.list_as_df <- function (x = list(), nrow = NULL) 
{
  stopifnot(is.list(x), is.null(nrow) || nrow >= 0L)
  if (n <- length(x)) {
    if (is.null(nrow)) 
      nrow <- max(lengths(x), 0L)
    x <- lapply(x, rep_len, nrow)
  }
  else {
    if (is.null(nrow)) 
      nrow <- 0L
  }
  if (is.null(names(x))) 
    names(x) <- character(n)
  class(x) <- "data.frame"
  attr(x, "row.names") <- .set_row_names(nrow)
  x
}

# Wrapper for .clip_edge for clipping polygons/tads
.clip_poly <- function(x, y, id = 1, 
                       xrange = c(-Inf, Inf), yrange = c(-Inf, Inf)) {
  # Clip right
  xy <- .clip_edge(x, y, xrange[2])
  # Clip left, flip x as edge is clipped on the right
  xy <- .clip_edge(-xy$x, xy$y, -xrange[1])
  # Clip top, switch x/y
  xy <- .clip_edge(xy$y, xy$x, yrange[2])
  # Clip bottom, flip x, which is the old y
  xy <- .clip_edge(-xy$x, xy$y, -yrange[1])
  # Switch and flip back x/y
  return(list(x = -xy$y, y = -xy$x))
}

# Only clips right edges, transform data prior to calling .clip_edge to make
# it seem like you've clipped a different edge
.clip_edge <- function(x, y, clip) {
  # Reparameterise as stretches of out of boundness
  rle <- rle(x > clip)
  # If nothing is out of bounds, skip clipping
  if (any(rle$values)) {
    # Get runlength parameters
    starts <- {ends <-  cumsum(rle$lengths)} - rle$lengths + 1
    is_outside <- which(rle$values)
    
    # Match inside/outside points
    inside  <- rbind(starts - 1, ends + 1)[, is_outside, drop = FALSE]
    outside <- rbind(starts, ends)[, is_outside, drop = FALSE]
    
    # Drop global polygon starts and ends
    keep <- !(inside > sum(rle$lengths) | inside < 1)
    inside  <- inside[keep]
    outside <- outside[keep]
    
    # Track IDs
    vrtx_idx <- split(cumsum(keep)[keep], col(keep)[keep])
    
    # Calculate weights for y
    weight <- x[inside]
    weight <- (clip - weight) / (x[outside] - weight)

    # Weight y
    weighted <- y[inside]
    weighted <- weighted + weight * (y[outside] - weighted)

    # Setup new vectors
    rle$lengths[is_outside] <- lengths(vrtx_idx)
    new_x <- new_y <- numeric(sum(rle$lengths))
    new_starts <- {new_ends <- cumsum(rle$lengths)} - rle$lengths + 1
    new_vrtx <- {idx <- seq_along(new_ends)} %in% is_outside
    
    # Fill in new coordinates
    j <- 1
    for (i in idx) {
      new_i <- new_starts[i]:new_ends[i]
      if (new_vrtx[i]) {
        # Place new vertices
        new_x[new_i] <- clip
        new_y[new_i] <- weighted[vrtx_idx[[j]]]
        j <- j + 1
      } else {
        # Keep old vertices
        old_i <- starts[i]:ends[i]
        new_x[new_i] <- x[old_i]
        new_y[new_i] <- y[old_i]
      }
    }
    x <- new_x
    y <- new_y
  }
  return(list(x = x, y = y))
} 