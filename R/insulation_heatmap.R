#' Plot a heatmap or profile of insulation scores
#'
#' Takes pre-calculated insulation scores around a set of genomic locations and
#' plots these in a heatmap and/or profile.
#'
#' @inheritParams call_TAD_insulation
#' @param bed A \code{data.frame} with 3 columns in BED format, containing
#'   genomic locations for which the insulation score is taken.
#' @param bed_pos A \code{character} vector of length 1 containing either
#'   \code{"start"}, \code{"end"} or \code{"center"} indicating what position of
#'   the '\code{bed}' argument to center insulation scores around.
#' @param region_width A \code{numeric} vector of length 1 for the size in
#'   basepairs around the feature to retrieve insulation scores for.
#' @param sort_bins A \code{integer} of length 1 indicating how many central
#'   bins should be used to sort the heatmap.
#' @param sort_on A \code{integer} of length 1 noting an experiment on which to
#'   sort the heatmap.
#' @param mode A \code{character} of length 1 containing either
#'   \code{"heatmap"}, \code{"profile"} or \code{"both"} for how the data should
#'   be plotted.
#' @param colour_lim A \code{numeric} vector of length 2 indicating the limits
#'   of the heatmap colour scale in insulation scores.
#' @param colours A vector of valid colours for the heatmap gradient. Internally
#'   defaults to a blue-white-red gradient.
#' @param profile_fun A \code{function} with which to calculate the profile.
#'   Should return 1 value per call.
#' @param custom_fill_scale A call to one of ggplot2's \code{scale_fill_*}
#'   family of functions to use as the heatmap's scale. If this argument is not
#'   \code{NULL}, it overrides the '\code{colour_lim}' and '\code{colour}'
#'   arguments.
#' @param title A \code{character} or \code{expression} of length 1 to use as
#'   text for the plot title.
#'
#' @details The profile colours are automatically assumed from the colours given
#'   to the \code{contacts} object during the \code{load_contacts} operation.
#'   Entries in the '\code{bed}' argument that contain > 1 percent missing data
#'   are dropped.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @seealso \code{\link[GENOVA]{insulation_score}} for calculating insulation
#'   scores. \code{\link[GENOVA]{call_TAD_insulation}} for calling TADs on the
#'   insulation score.
#'
#' @examples
#' \dontrun{
#' # Calculating insulation scores
#' ins <- insulation_score(list(WT_20kb, KO_20kb), window = 25)
#'
#' # Calling TADs from the insulation score
#' tadlist <- call_TAD_insulation(ins)
#'
#' # Plotting a heatmap
#' heatmap_insulation(ins, tadlist$WT_20kb)
#' }
heatmap_insulation <- function(IS_discovery, 
                               bed, bed_pos = "end", 
                               region_width = 1e6,
                               sort_bins = 11,
                               sort_on = 1,
                               mode = c("both", "heatmap", "profile"),
                               colour_lim = c(-1, 1),
                               colours = NULL,
                               profile_fun = mean,
                               custom_fill_scale = NULL,
                               title = NULL) {
  # Setting up parameters
  mode <- match.arg(mode, c("both", "heatmap", "profile"))
  
  # Experiment parameters
  ins <- IS_discovery$insula_score
  expcols <- attr(IS_discovery, "colours")
  expnames <- tail(colnames(ins), -4)
  res <- attr(IS_discovery, "resolution")
  
  # Profile function parse name
  funname <- as.character(substitute(profile_fun))[[1]]
  funname <- strsplit(funname, "")[[1]]
  funname[1] <- toupper(funname[1])
  funname <- paste0(funname, collapse = "")
  
  # Grab relevant data
  binwidth <- region_width / res + 1
  
  # Translate bed to indices
  features <- data.table(feat = bed2idx(ins[,1:4], bed, mode = bed_pos))
  features[, i := seq_len(.N)]
  
  # Determine region around bed entries
  idx <- seq_len(binwidth) - ceiling(0.5 * binwidth)
  idx <- features[, list("bin" = feat + idx), by = i]
  idx[, x := seq_len(.N), by = i]
  
  # Grab that region
  dat <- ins[idx, on = c("bin" = "bin")]
  dat <- melt(dat, id.vars = c("i", "x"), measure.vars = expnames)

  # Censor missing data
  missing <- dat[, mean(is.na(value) | !is.finite(value)), by = i]
  missing <- missing$i[missing$V1 > 0.01]
  dat <- dat[!(i %in% missing)]
  dat$value[is.na(dat$value)] <- 0
  
  # Find order of rows
  sorter <- (binwidth / 2) + c(-0.5, 0.5) * sort_bins
  sorter <- dat[between(x, sorter[1], sorter[2])]
  sorter <- sorter[variable == expnames[sort_on]]
  sorter <- sorter[, list(mu = mean(pmax(value, - 2), na.rm = TRUE)), by = i]
  sorter[, c("ord", "mu") := list(order(order(mu, decreasing = TRUE)), NULL)]
  
  # Combine heatmap with row orders
  dat <- dat[sorter, on = c("i" = "i")]
  dat$facet <- "Heatmap\n(Location Rank)"
  
  # Calculate profile
  prof <- dat[, list(value = profile_fun(value)), by = c("x", "variable")]
  prof$facet <- paste0("Profile\n(", funname, " Insulation)")
  
  # Expansion factor for profile so the x-axis is not touched
  efactor <- diff(range(prof$value)) * c(-0.6, 0.6) + mean(range(prof$value))
  
  # Take common factors
  g <- ggplot2::ggplot(dat, ggplot2::aes((x - 0.5 * binwidth) * res)) +
    ggplot2::scale_x_continuous(name = "Position relative to feature (kb)",
                                labels = function(x){x/1000},
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), name = "")
  
  # Take conditional factors
  if (mode %in% c("both", "heatmap")) {
    g <- g + ggplot2::geom_raster(ggplot2::aes(y = ord, fill = value)) +
      ggplot2::geom_blank(data = data.frame(facet = dat$facet[1],
                                            x = 0, y = 0),
                          ggplot2::aes(x, y))
    if (!is.null(custom_fill_scale)) {
      g <- g + custom_fill_scale
    } else if (!is.null(colours)) {
      g <- g + ggplot2::scale_fill_gradientn(colours = colours,
                                             name = "Insulation\nScore",
                                             limits = colour_lim,
                                             oob = scales::squish)
    } else {
      g <- g +
        ggplot2::scale_fill_gradient2(low = "#ff5c49", high = "#009bef",
                                      name = "Insulation\nScore",
                                      limits = colour_lim, oob = scales::squish)
    }
  }
  if (mode %in% c("both", "profile")) {
    g <- g + ggplot2::geom_line(data = prof, 
                                ggplot2::aes(y = value, colour = variable)) +
      ggplot2::geom_blank(data = data.frame(y = efactor, facet = prof$facet[1]), 
                          ggplot2::aes(x = 0, y = y)) +
      ggplot2::scale_colour_manual(values = expcols, name = "Experiment")
  }
  if (mode == "both") {
    g <- g + ggplot2::facet_grid(facet ~ variable, 
                                 scales = "free_y", switch = "y") +
      ggplot2::theme(strip.placement = "outside")
  } else if (mode == "heatmap") {
    g <- g + ggplot2::facet_grid(facet ~ variable, switch = "y")
  } else if (mode == "profile") {
    suppressMessages(
      g <- g + ggplot2::scale_y_continuous(name = prof$facet[1])
    )
    
  }
  if (!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }
  
  # Take common factors again
  g <- g + ggplot2::theme(
    strip.placement = "outside",
    axis.ticks = ggplot2::element_line(colour = "black"),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.line = ggplot2::element_line(colour = "black"),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    panel.spacing.x = grid::unit(0.7, "strwidth", "-500"),
    legend.key = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(size = ggplot2::rel(1.25))
  )
  g
}
