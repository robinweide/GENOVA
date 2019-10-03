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
#' @param metric A \code{character} of length 1:
#'
#'   \describe{ \item{A*A}{\code{"diff"} for difference by subtraction or
#'   \code{"lfc"} for \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} fold
#'   changes.} \item{RCP}{\code{"smooth"} for a \ifelse{html}{\out{log<sub>10</sub>}}{\eqn{log_10}}-smoothed line, \code{both}
#'   for adding the raw distance-bins as points and \code{lfc} for
#'   \ifelse{html}{\out{log<sub>2</sub>}}{\eqn{log_2}} fold changes.} }
#'
#' @param colour_lim,colour_lim_contrast One of: \itemize{
#'   \item \code{NULL} to have the middle colour fall on a sensible mid-point.
#'   \item A \code{numeric} vector of length two providing the limits of the scale. Use \code{NA} to refer to existing minima or maxima.
#'   \item A \code{function} that accepts the existing (automatic) limits and returns new limits.
#' }
#'
#' @param raw A \code{logical} of length 1: should a bare bones plot be
#'   returned?
#'
#' @param bins \code{[virtual4C]} Set the number of histogram-bins
#'
#' @param bed \code{[virtual4C]} A data.frame of bed-entries to plot underneath.
#'
#' @param extend_viewpoint \code{[virtual4C]} Add a set bp to both sides of the
#'   viewpoint. Makes the viewpoint-box broader.
#'
#' @param mode \code{[PESCAn & ARA]} What result slot should be used for
#'   visualisation? \code{"obsexp"} for the observed over expected metric or
#'   \code{"signal"} for mean contacts at unshifted anchors.
#'
#' @param flipFacet \code{[RCP]} Do you want to have RCP's of different
#'   regions in one plot, instead of facets? (default : \code{FALSE})
#' @param chr \code{[CS]} A \code{character} of length 1 indicating a
#'   chromosome name.
#' @param start,end \code{[CS]} A \code{numeric} of length 1 with start-
#'   and end-positions for the region to plot. If \code{NULL}, is set to
#'   \code{-Inf} and \code{Inf} respectively.
#'
#'
#' @param flipFacet \code{[RCP]} Do you want to have RCP's of different regions
#'   in one plot, instead of facets? (default : \code{FALSE})
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
#' @note For \code{ATA_discovery} objects which include the matrix's diagonal,
#'   the upper limit for the contacts fill scale is set to the 95th percentile
#'   of the data to increase the dynamic range of colours. The same is true for
#'   \code{ARA_discovery}, when '\code{mode = "signal"}'.
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
                                    raw = FALSE, title = NULL) {
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

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  g
}

#' @rdname visualise
#' @export
visualise.PESCAn_discovery <- function(discovery, contrast = 1,
                                       metric = c("diff", "lfc"),
                                       mode = c("obsexp", "signal"),
                                       raw = FALSE, title = NULL) {
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

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  
  g
}

#' @rdname visualise
#' @export
visualise.ATA_discovery <- function(discovery, contrast = 1,
                                    metric = c("diff", "lfc"),
                                    raw = FALSE, title = NULL) {
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

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  
  g
}

#' @rdname visualise
#' @export
visualise.ARA_discovery <- function(discovery, contrast = 1,
                                    metric = c("diff", "lfc"),
                                    mode = c("obsexp", "signal"),
                                    raw = FALSE, title = NULL) {
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
  g <- GENOVA:::visualise.ARMLA(
    discovery = res,
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
      name = expression(mu*" Contacts"),
      limits = c(NA, upperq),
      oob = scales::squish
    )
  }

  g <- g + fillscale +
    ggplot2::scale_x_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        paste0(x / 1000, "kb")
      }
    ) +
    ggplot2::scale_y_continuous(
      name = "",
      expand = c(0, 0),
      breaks = pos_breaks,
      labels = function(x) {
        paste0(x / 1000, "kb")
      }
    )

  if(!is.null(title)){
    g = g + ggplot2::ggtitle(title)
  }
  
  g
}

#' @rdname visualise
#' @export
visualise.RCP_discovery = function(discovery, contrast = 1, metric = c("smooth","both","lfc"), raw = F, title = NULL, flipFacet = F){
  
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
#' @export
visualise.CS_discovery <- function(discovery, contrast = NULL,
                                   chr = "chr1", start = NULL, end = NULL,
                                   raw = FALSE) {
  start <- if (is.null(start)) -Inf else start
  end <- if (is.null(end)) Inf else end
  df <- discovery$compart_scores
  ii <- which(df[["chrom"]] == chr & df[["start"]] >= start & 
                df[["end"]] <= end)
  df <- df[ii,]
  expnames <- colnames(df)[5:ncol(df)]
  
  df <- data.frame(mid = (df[["start"]] + df[["end"]])/2,
                   df[, ..expnames])
  
  if (!is.null(contrast)) {
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
  
  if (!is.null(contrast)) {
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
    labels = function(x){paste0(x/1e6, " Mb")}
    ) +
    ggplot2::scale_y_continuous(
      name = yname,
      # limits = c(-1.25, 1.25),
      oob = scales::squish,
      limits = function(x){c(-1, 1) * max(abs(x))}
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
visualise.virtual4C_discovery <- function(discovery, bins = NULL, bed = NULL, extend_viewpoint = NULL){
  # ! someday: allow mulitple samples
  data <- discovery$data
  VP   <- attr(discovery,"viewpoint")
  
  if(!is.null(extend_viewpoint)){
    VP[1,2] <- VP[1,2] - extend_viewpoint
    VP[1,3] <- VP[1,3] + extend_viewpoint
  }
  
  if( is.null(bins) ) {
    bins = nrow(data)
  }
  
  VP_mid <- rowMeans(attr(discovery, 'viewpoint')[1,2:3])
  
  blackout_up   <- VP[1,2]
  blackout_down <- VP[1,3]
  data_blackout <- data[(data$mid >= blackout_up & data$mid <= blackout_down)]
  
  if( !is.null(bed)) {
    bed = bed[bed[,1] == attr(discovery, 'viewpoint')[1,1],2:3]
  }
  
  data <- data[!(data$mid >= blackout_up & data$mid <= blackout_down)]
  
  breaks   <- seq(min(data$mid), max(data$mid), length.out = bins)
  bin_size <- median(diff(breaks))
  smooth   <- data[, mean(signal),by = findInterval(data$mid, breaks)]
  smooth$mid = breaks[unlist(smooth[,1])] +(bin_size/2)
  smooth[,1] = NULL
  colnames(smooth) = c("signal","mid")
  
  p = ggplot2::ggplot(data, ggplot2::aes(x= mid/1e6, y = signal)) +
    ggplot2::geom_col(data = smooth, fill = 'black', width = bin_size/1e6) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = attr(discovery, 'viewpoint')[1,1])
  
  # draw_blackout ===+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ymax <- ceiling(max(smooth$signal))
  data_blackout$signal <- ymax

  p = p + ggplot2::annotate('rect', 
                            fill = "#D8D8D8", 
                            xmin = (min(data_blackout$mid)/1e6)-(bin_size/1e6), 
                            xmax = (max(data_blackout$mid)/1e6)+(bin_size/1e6), 
                            ymin = 0, 
                            ymax = max(data_blackout$signal) )+
    ggplot2::annotate(geom = 'text',
                      vjust = 1,
                      x = rowMeans(VP[,2:3])/1e6,
                      y =  ymax*0.9,
                      label = '\u2693')

  if( !is.null(bed)){
    p = p + ggplot2::annotate('rect', 
                              fill = "black",
                              xmin = bed[,1]/1e6, 
                              xmax = bed[,2]/1e6,
                              ymin = -ceiling(max(smooth$signal))/100,
                              ymax = 0)  
  }
  p
  bed_track_size <- ceiling(max(smooth$signal))/50
  bed_track_padding <- ceiling(max(smooth$signal))/200
  
  if( !is.null(bed) ){
    bed <- bed/1e6
    
    p = p+ ggplot2::annotate(geom = 'rect', 
                             xmin = bed[,1], 
                             xmax = bed[,2], 
                             ymin = -bed_track_size, 
                             ymax = -bed_track_padding, 
                             fill = 'black') +
      ggplot2::coord_cartesian(expand = F, 
                               ylim = c(-1*(bed_track_size+
                                              bed_track_padding+
                                              bed_track_padding), 
                                        ymax))
  } else {
    p = p + ggplot2::coord_cartesian(expand = F, 
                                     ylim = c(0, ymax))
  }

  p <- p + ggplot2::theme(axis.line = ggplot2::element_line(colour = 'black'),
                          axis.text = ggplot2::element_text(colour = 'black'))
  suppressWarnings(p)
}

#' @rdname visualise
#' @export
visualise.saddle_discovery <- function(discovery, contrast = 1,
                                       raw = FALSE, title = NULL,
                                       colour_lim = NULL,
                                       colour_lim_contrast = NULL) {
  df <- discovery$saddle
  df <- df[!is.na(mean) & !is.na(q1),]
  
  expnames <- df[, unique(exp)]
  
  # Aggregate across chromosomes
  df <- df[, mean(mean), by = list(exp, q1, q2)]
  setnames(df, 4, "obsexp")
  
  # Mirror along diagonal
  comp <- df[q1 != q2, ]
  setnames(comp, 2:3, names(comp)[3:2])
  df <- rbindlist(list(df, comp), use.names = TRUE)
  
  if (is.null(colour_lim)) {
    colour_lim <- function(x) {c(-1, 1) * max(abs(x - 1)) + 1}
  }
  if (is.null(colour_lim_contrast)) {
    colour_lim_contrast <- function(x){c(-1, 1) * max(abs(x))}
  }

  if (!is.null(contrast)) {
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
  
  if (!is.null(contrast)) {
    
    suppressWarnings(
      g <- g + ggplot2::geom_raster(data = contrast, 
                                    ggplot2::aes(q1, q2, altfill = diff),
                                    inherit.aes = FALSE,
                                    hjust = 0, vjust = 0) +
        ggplot2::facet_grid(panel ~ exp, switch = "y")
    )
    
    if (!raw) {
      g <- g + ggplot2::scale_fill_gradientn(
        colours = c("#23AA17", "#90D48E", "#FFFFFF", "#F0A9F1", "#DA64DC"),
        aesthetics = "altfill",
        oob = scales::squish,
        name = "Difference",
        limits = colour_lim_contrast
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
                                breaks = c(0, 0.5, 1)) +
      ggplot2::scale_y_continuous(name = "Quantile", expand = c(0,0),
                                  breaks = c(0, 0.5, 1)) +
      ggplot2::scale_fill_gradientn(
        colours = c("#009BEF", "#7FCDF7", "#FFFFFF", "#FFADA3", "#FF5C49"),
        limits = colour_lim,
        oob = scales::squish,
        name = expression(frac("Observed", "Expected")),
        guide = ggplot2::guide_colourbar(order = 1)) +
      ggplot2::coord_cartesian(clip = "off")
  }
  return(g)
}

# Utilities ---------------------------------------------------------------

# Makes sure no errors are returned when visualise(..., raw = TRUE)
# Unfortunately has to be exported, but is better than throwing errors
#' @export
#' @keywords internal
scale_altfill_continuous <- function(...) {
  ggplot2::scale_fill_gradient2(..., aesthetics = "altfill")
}
