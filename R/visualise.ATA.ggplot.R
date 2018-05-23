#' visualise.ATA.ggplot
#'
#' Plot the stacked TAD results: takes in a named list of ATA-outputs and plots these next to eachother.
#' It also plots differentials.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param stackedlist Named list of results from \code{ATA}.
#' @param title Title-text to plot
#' @param Focus The index of the sample in \code{stackedlist} thats use to compare against.
#' @param zlim1 Vector of 2 values, denoting the min and max color-values of the first row
#' @param zlim2 Vector of 2 values, denoting the min and max color-values of the second row
#' @details
#' By substracting the values from each sample with the values from the \code{focus}-sample, we generate the differentials.
#' A positive value in the differential plots thus means an enrichment in that sample versus the \code{focus}-sample.
#' @return A \code{grid} object, containing two ggplot-objects.
#' @examples
#' # run ATA
#' ATA_RAOetal <- ATA(experiment = Rao_20k,tad.bed = TADs)
#' ATA_Sanbornetal <- ATA(experiment = Sanborn_20k,tad.bed = TADs)
#'
#' # make a named list
#' ATAlist <- list(Rao = ATA_RAOetal, Sanborn = ATA_Sanbornetal)
#'
#' # plot the results
#' visualise.ATA.ggplot(ATAlist, focus = 1, title = 'Rao vs Sanborn.')
#' @export
visualise.ATA.ggplot <- function(stackedlist, title = "ATA", focus = 1, zlim1 = NULL, zlim2 = NULL){

  #! check OUTLIERCORRECTIONSWITCH
  OL = c()
  for(i in 1:length(stackedlist)){
    OL = c(OL, stackedlist[[i]]$OUTLIERCORRECTIONSWITCH)
  }
  OL = unique(OL)
  if(length(OL) != 1){
    stop('There seem to be some ATAs with outlier-correction and some not...')
  }

  # issue 29: rename stackr -> ata
  # Make two dataframes
  higlassCol <- c('white', '#f5a623', '#d0021b', 'black')
  abovePlots <- data.frame(Var1 = integer(),
                           Var2 = integer(),
                           value = numeric(),
                           sample = factor())
  belowPlots <- abovePlots

  # Loop trough list and melt
  list.len <- length(stackedlist)
  belownames <- vector()
  m <- vector()
  for(i in 1:list.len){
    m <- c(m,max(stackedlist[[i]]$STACK[30:40,60:70]))
    a <- reshape2::melt(stackedlist[[i]]$STACK)
    a$sample <- factor(rep(names(stackedlist)[i], length(a[,1])))
    abovePlots <- rbind(abovePlots, a)

    e <- reshape2::melt(stackedlist[[i]]$STACK-stackedlist[[focus]]$STACK)
    e$sample <- rep(paste0(names(stackedlist)[focus], ' vs ', names(stackedlist)[i] ) , length(e[,1]))
    belowPlots <- rbind(belowPlots, e)
    belownames <- c(belownames, paste0(names(stackedlist)[focus], ' vs ', names(stackedlist)[i] ))
  }
  colnames(abovePlots) <- c('Var1','Var2','value','sample')
  colnames(belowPlots) <- c('Var1','Var2','value','sample')

  # set zlims
  # issue 28: GGatavis: manual zscale
  if(is.null(zlim1)){
    zminAbove <- min(na.exclude(abovePlots$value))
    zmaxAbove <- max(na.exclude(abovePlots$value))
  } else {
    zminAbove <- zlim1[1]
    zmaxAbove <- zlim1[2]
  }
  if(is.null(zlim2)){
    z2 <- min( abs(c(  quantile(na.exclude(belowPlots$value), 0.005) , quantile(na.exclude(belowPlots$value), .995)        ) )   )
    zminBelow <- z2 * -1
    zmaxBelow <- z2
  } else {
    zminBelow <- zlim2[1]
    zmaxBelow <- zlim2[2]
  }

  # issue 19: don't allow diff0 as diff-zlim
  if(diff(c(zminBelow,zmaxBelow))  == 0){
    warning("The Z-lim of the differentials is too close together.\nUse zlim2 to override this.")
    zminBelow <- -1
    zmaxBelow <- 1
  }

  # To avoid bleaching, set out-of-zlim values to min/max
  abovePlots[abovePlots$value < zminAbove, 3] <- zminAbove
  abovePlots[abovePlots$value > zmaxAbove, 3] <- zmaxAbove
  belowPlots[belowPlots$value < zminBelow, 3] <- zminBelow
  belowPlots[belowPlots$value > zmaxBelow, 3] <- zmaxBelow

  # Set focus-sample to column 1
  abovePlots$sample <- factor(abovePlots$sample, levels = levels(abovePlots$sample)[match(levels(abovePlots$sample) ,  names(stackedlist))])
  volgorde <- match(names(stackedlist),levels(abovePlots$sample))
  belowPlots$sample <- factor(belowPlots$sample, levels = belownames[volgorde] )

  # Plot first row
  upTick = ((1/99) *25)*100
  downTick = ((1/99) *75)*100
  shimmy = ((0.5/99) )*100

  plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value),interpolate= F) +
    ggplot2::coord_fixed(expand = F) +
    ggplot2::scale_fill_gradientn(colours = higlassCol, limits = c(zminAbove, zmaxAbove))+
    ggplot2::scale_x_continuous(expand = c(0,0),
                                breaks = c(upTick+shimmy,downTick-shimmy), trans = 'reverse',labels = c("3' border", "5' border"))+
    ggplot2::scale_y_continuous(expand = c(0,0),
                                breaks = c(upTick+shimmy,downTick-shimmy), labels = c("3' border", "5' border"))+
    GENOVA_THEME() +
    ggplot2::facet_grid(.~ sample)+
    ggplot2::labs(title = title,x = "", y = "", fill = "Contacts ")


  plot2 <- ggplot2::ggplot(belowPlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value), interpolate = F) +
    ggplot2::facet_grid(.~ sample) +
    ggplot2::scale_fill_gradient2(limits = c(zminBelow, zmaxBelow), low="#2166ac", mid="white", high="#b2182b") +
    ggplot2::coord_fixed(expand = F) +
    ggplot2::scale_x_continuous(expand = c(0,0),
                                breaks = c(upTick+shimmy,downTick-shimmy), trans = 'reverse', labels = c("3' border", "5' border"))+
    ggplot2::scale_y_continuous(expand = c(0,0),
                                breaks = c(upTick+shimmy,downTick-shimmy), labels = c("3' border", "5' border"))+
    GENOVA_THEME() +
    ggplot2::labs(title = '',x = "", y = "", fill = "Difference")

  grid::grid.newpage()
  grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2), size = "last"))
  }
