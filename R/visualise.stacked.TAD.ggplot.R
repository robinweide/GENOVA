#' Plot the stacked TAD results.
#'
#' @param stackedlist List of results from `stacked.TAD`.
#' @param title Text to plot
#' @param Focus Wich sample does need to be the to-compare sample?
#' @return A grid object, containing two ggplot-objects.
#' @export
visualise.stacked.TAD.ggplot <- function(stackedlist, title, focus = 1){
  # Make two dataframes
  abovePlots <- data.frame(Var1 = integer(),
                           Var2 = integer(),
                           value = numeric(),
                           sample = factor())
  belowPlots <- abovePlots
  # Loop trough list and melt
  list.len <- length(stackedlist)
  belownames <- vector()
  for(i in 1:list.len){
    a <- reshape2::melt(stackedlist[[i]])
    a$sample <- factor(rep(names(stackedlist)[i], length(a[,1])))
    abovePlots <- rbind(abovePlots, a)

    e <- reshape2::melt(stackedlist[[i]]-stackedlist[[focus]])
    e$sample <- rep(paste0(names(stackedlist)[focus], ' vs ', names(stackedlist)[i] ) , length(e[,1]))
    belowPlots <- rbind(belowPlots, e)
    belownames <- c(belownames, paste0(names(stackedlist)[focus], ' vs ', names(stackedlist)[i] ))
  }
  colnames(abovePlots) <- c('Var1','Var2','value','sample')
  colnames(belowPlots) <- c('Var1','Var2','value','sample')
  # Find zlims
  zmaxAbove <- quantile(na.exclude(abovePlots$value), .85)
  zminAbove <- quantile(na.exclude(abovePlots$value), .15)
  abovePlots[abovePlots$value < zminAbove, 3  ] <- zminAbove
  abovePlots[abovePlots$value > zmaxAbove,  3 ] <- zmaxAbove
  z2 <- min( abs(c(  quantile(na.exclude(belowPlots$value), .035) , quantile(na.exclude(belowPlots$value), .965)        ) )   )
  belowPlots$value[belowPlots$value > z2] <- z2
  belowPlots$value[belowPlots$value < (z2 *-1)] <- (z2 *-1)
  # Set focus-sample to column 1
  abovePlots$sample <- factor(abovePlots$sample, levels = c(levels(abovePlots$sample)[focus], levels(abovePlots$sample)[!levels(abovePlots$sample) %in% levels(abovePlots$sample)[focus]] ))
  volgorde <- match(names(stackedlist),levels(abovePlots$sample))
  belowPlots$sample <- factor(belowPlots$sample, levels = belownames[volgorde] )
  # Plot first row
  plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) + ggplot2::geom_raster(ggplot2::aes(fill = value),interpolate= T) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_distiller(limits = c(zminAbove, zmaxAbove), palette = "Spectral")  +
    ggplot2::scale_x_continuous(breaks = c(26,77), trans = 'reverse',labels = c("3' border", "5' border"))+
    ggplot2::scale_y_continuous(breaks = c(26,77), labels = c("3' border", "5' border"))+
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA",colour = NA)) +ggplot2::facet_grid(.~ sample)+
    #geom_segment(aes(y = 26, x = 0, yend = 26, xend = 100, colour = "segment"))+
    #geom_segment(aes(y = 77, x = 0, yend = 77, xend = 100, colour = "segment"))
    ggplot2::labs(title = title,x = "", y = "", fill = "Contacts ")
  plot2 <- ggplot2::ggplot(belowPlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
    ggplot2::facet_grid(.~ sample) +
    ggplot2::scale_fill_gradient2(limits = c(z2*-1,z2), midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(breaks = c(26,77), trans = 'reverse', labels = c("3' border", "5' border"))+
    ggplot2::scale_y_continuous(breaks = c(26,77), labels = c("3' border", "5' border"))+
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA",colour = NA)) +
    ggplot2::labs(title = '',x = "", y = "", fill = "Difference")
  grid::grid.newpage()
  grid::grid.draw(rbind(grid::ggplotGrob(plot1), ggplotGrob(grid::plot2), size = "last"))
}
