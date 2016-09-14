#' Plot the APA-results
#'
#' @param APAlist A list of results from `APA`.
#' @param title Text to plot
#' @param Focus Wich sample does need to be the to-compare sample?
#' @return A grid object, containing two ggplot-objects.
#' @export
visualise.APA.ggplot <- function(APAlist, title, zCoef = 1, focus = 1, ...){
  # Make two dataframes
  abovePlots <- data.frame(Var1 = integer(),
                           Var2 = integer(),
                           value = numeric(),
                           sample = factor())
  belowPlots <- abovePlots
  # Loop trough list and melt
  list.len <- length(APAlist)
  belownames <- vector()
  for(i in 1:list.len){
    a <- reshape2::melt(APAlist[[i]])
    a$sample <- factor(rep(names(APAlist)[i], length(a[,1])))
    abovePlots <- rbind(abovePlots, a)
    
    e <- reshape2::melt(APAlist[[i]]-APAlist[[focus]])
    e$sample <- rep(paste0(names(APAlist)[focus], ' vs ', names(APAlist)[i] ) , length(e[,1]))
    belowPlots <- rbind(belowPlots, e)
    belownames <- c(belownames, paste0(names(APAlist)[focus], ' vs ', names(APAlist)[i] ))
  }
  # Find best color-limits
  z <- max( abs( c(  quantile(na.exclude(abovePlots$value), .01) , quantile(na.exclude(abovePlots$value), .99)        )  )   )*zCoef
  abovePlots$value[abovePlots$value > z] <- z
  abovePlots$value[abovePlots$value < (z *-1)] <- (z *-1)
  # Set focus-sample to column 1
  abovePlots$sample <- factor(abovePlots$sample, levels = c(levels(abovePlots$sample)[focus], levels(abovePlots$sample)[!levels(abovePlots$sample) %in% levels(abovePlots$sample)[focus]] ))
  volgorde <- match(names(APAlist),levels(abovePlots$sample))
  belowPlots$sample <- factor(belowPlots$sample, levels = belownames[volgorde] )
  # Plot per sample
  plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) + ggplot2::geom_raster(ggplot2::aes(fill = value)) + ggplot2::facet_grid(.~ sample) +
    ggplot2::coord_fixed() + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA",colour = NA)) +
    ggplot2::scale_x_continuous(breaks = c(5.5,11,16.5), labels = c("-50Kb", "3'", "+50Kb")) +
    ggplot2::scale_y_continuous(breaks = c(5.5,11,16.5), labels = c("-50Kb", "5'", "+50Kb")) +
    ggplot2::labs(title = title,x = "", y = "", fill = "Contacts ") +
    ggplot2::scale_fill_distiller(limits = c(0, z), palette = "Spectral")
  # Find best color-limits
  z2 <- max( abs( c(  quantile(na.exclude(belowPlots$value), .01) , quantile(na.exclude(belowPlots$value), .99)        )  )   )*zCoef
  belowPlots$value[belowPlots$value > z2] <- z2
  belowPlots$value[belowPlots$value < (z2 *-1)] <- (z2 *-1)
  # Plot difference-plots
  plot2 <- ggplot2::ggplot(belowPlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) + ggplot2::facet_grid(.~ sample) + ggplot2::coord_fixed()  +
    ggplot2::scale_fill_gradient2(limits = c(z2*-1,z2), midpoint=0, low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#FAFAFA",colour = NA))+
    ggplot2::scale_x_continuous(breaks = c(5.5,11,16.5), labels = c("-50Kb", "3'", "+50Kb"))+
    ggplot2::scale_y_continuous(breaks = c(5.5,11,16.5), labels = c("-50Kb", "5'", "+50Kb")) +
    ggplot2::labs(x = "", y = "", fill = "Difference")
  return(grid.arrange(plot1, plot2, ncol=1))
  #grid::grid.newpage()
  #grid::grid.draw(rbind(grid::ggplotGrob(plot1), grid::ggplotGrob(plot2), size = "last"))
}
