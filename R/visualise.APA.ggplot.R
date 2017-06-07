#' Plot the APA-results
#'
#' @param APAlist A list of results from `APA`.
#' @param title Text to plot
#' @param Focus Which sample does need to be the to-compare sample?
#' @param zTop The min and max colorscale-values for the first row of plots.
#' @param zBottom The min and max colorscale-values for the first row of plots.
#' @return A grid object, containing two ggplot-objects.
#' @export
visualise.APA.ggplot <- function(APAlist, resolution, title = 'APA', zTop = NULL, zBottom = NULL, focus = 1,...){
  size <- dim(as.data.frame(APAlist[[1]][1]))[1]
  size.banks <- (size - 1)/2
  tickLabelDownstream <- as.character(1 * ((size.banks/2 *
                                              resolution)/1000))
  tickLabelUpstream <- as.character(-1 * ((size.banks/2 * resolution)/1000))
  abovePlots <- data.frame(Var1 = integer(), Var2 = integer(),
                           value = numeric(), sample = factor())
  belowPlots <- abovePlots
  list.len <- length(APAlist)
  belownames <- vector()
  for (i in 1:list.len) {
    firstMat <- as.data.frame(APAlist[[i]][1]); colnames(firstMat) <- 1:size ; rownames(firstMat) <- 1:size
    secondMat <- as.data.frame(APAlist[[focus]][1]); colnames(secondMat) <- 1:size ; rownames(secondMat) <- 1:size
    a <- reshape2::melt(as.matrix(firstMat))
    a$sample <- factor(rep(names(APAlist)[i], length(a[, 1])))
    abovePlots <- rbind(abovePlots, a)

    e <-  melt(as.matrix(firstMat - secondMat))
    e$sample <- rep(paste0(names(APAlist)[focus], " vs ",
                           names(APAlist)[i]), length(e[, 1]))
    belowPlots <- rbind(belowPlots, e)
    belownames <- c(belownames, paste0(names(APAlist)[focus],
                                       " vs ", names(APAlist)[i]))
  }

  z <- NULL
  if(is.null(zTop) ){
    z <- c(quantile(na.exclude(abovePlots$value), 0.01),
           quantile(na.exclude(abovePlots$value), 0.99))
  } else {
    z <- zTop
  }

  abovePlots$value[abovePlots$value > z[2]] <- z[2]
  abovePlots$value[abovePlots$value < z[1]] <- z[1]
  abovePlots$sample <- factor(abovePlots$sample, levels = c(levels(abovePlots$sample)[focus],
                                                            levels(abovePlots$sample)[!levels(abovePlots$sample) %in%
                                                                                        levels(abovePlots$sample)[focus]]))
  volgorde <- match(names(APAlist), levels(abovePlots$sample))
  belowPlots$sample <- factor(belowPlots$sample, levels = belownames[volgorde])

  if(require('viridis') == TRUE){
    plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) +
      ggplot2::geom_raster(ggplot2::aes(fill = value)) +
      ggplot2::facet_grid(. ~ sample) + ggplot2::coord_fixed() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA)) +
      ggplot2::scale_x_continuous(breaks = c(size *  0.25, size * 0.5, size * 0.75), labels = c(paste0(tickLabelUpstream,  "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
      ggplot2::scale_y_continuous(breaks = c(size *  0.25, size * 0.5, size * 0.75), labels = c(paste0(tickLabelUpstream,  "kb"), "5'", paste0(tickLabelDownstream, "kb"))) +
      ggplot2::labs(title = title, x = "", y = "", fill = "Contacts ") + viridis::scale_fill_viridis(limits = z, option = "inferno", direction = -1)

  } else {
    plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) +
      ggplot2::geom_raster(ggplot2::aes(fill = value)) +
      ggplot2::facet_grid(. ~ sample) + ggplot2::coord_fixed() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA)) +
      ggplot2::scale_x_continuous(breaks = c(size *  0.25, size * 0.5, size * 0.75), labels = c(paste0(tickLabelUpstream,  "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
      ggplot2::scale_y_continuous(breaks = c(size *  0.25, size * 0.5, size * 0.75), labels = c(paste0(tickLabelUpstream,  "kb"), "5'", paste0(tickLabelDownstream, "kb"))) +
      ggplot2::labs(title = title, x = "", y = "", fill = "Contacts ") + ggplot2::scale_fill_distiller(limits = z, palette = "Spectral")
  }


  z2 <- NULL
  if(is.null(zBottom) ){
    z2 <- c(quantile(na.exclude(belowPlots$value), 0.01),
           quantile(na.exclude(belowPlots$value), 0.99))
  } else {
    z2 <- zBottom
  }


  belowPlots$value[belowPlots$value > z2[2]] <- z2[2]
  belowPlots$value[belowPlots$value < z2[1]] <- z2[1]


  plot2 <- ggplot2::ggplot(belowPlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::facet_grid(. ~ sample) + ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(limits = z2, midpoint = 0, low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA", colour = NA)) +
    ggplot2::scale_x_continuous(breaks = c(size *  0.25, size * 0.5, size * 0.75), labels = c(paste0(tickLabelUpstream,  "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(size *  0.25, size * 0.5, size * 0.75), labels = c(paste0(tickLabelUpstream, "kb"), "5'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::labs(x = "",   y = "", fill = "Difference")

  grid::grid.newpage()
  grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2),
                        size = "last"))
}
