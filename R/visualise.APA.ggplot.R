#' visualise.APA.ggplot
#'
#' Plot the APA-results and the differential results.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param APAlist A list of results from the APA-function.
#' @param title Title-text to plot.
#' @param focus Which sample will be the to-compare sample in the differential row?
#' @param zTop The min and max values for the first row of plots.
#' @param zBottom The min and max values for the first row of plots.
#' @return A grid object, containing two ggplot-objects.
#' @details
#' By substracting the values from each sample with the values from the \code{focus}-sample, we generate the differentials.
#' A positive value in the differential plots thus means an enrichment in that sample versus the \code{focus}-sample.
#' @export
visualise.APA.ggplot <- function(APAlist, title = 'APA', zTop = NULL, zBottom = NULL, focus = 1,...){
  OL = c()
  for(i in 1:length(APAlist)){
    OL = c(OL, APAlist[[i]]$OUTLIERCORRECTIONSWITCH)
  }
  OL = unique(OL)
  if(length(OL) != 1){
    stop('There seem to be APAs with different outlier-corrections...')
  }

  higlassCol <- c('white', '#f5a623', '#d0021b', 'black')
  resolution = APAlist[[1]]$RES

  size <- dim(APAlist[[1]]$APA)[1]
  size.banks <- (size - 1)/2
  allTicks = seq(resolution*size.banks, -1*resolution*size.banks, length.out = size)/1e3

  tickLabelUpstream = 0
  tickLabelDownstream = 0

  tickPosUpstream =  (size.banks +1) - (size.banks/2)
  tickPosDownstream = (size.banks +1) + (size.banks/2)

  if(size.banks %% 2 == 0){ # like 21
    tickLabelDownstream = allTicks[tickPosDownstream]
    tickLabelUpstream = allTicks[tickPosUpstream]
  } else { # 19
    tickPosUpstream = tickPosUpstream - .5
    tickPosDownstream = tickPosDownstream + .5

    tickLabelDownstream = allTicks[tickPosDownstream]
    tickLabelUpstream = allTicks[tickPosUpstream]
  }




  abovePlots <- data.frame(Var1 = integer(), Var2 = integer(),
                           value = numeric(), sample = factor())
  belowPlots <- abovePlots
  list.len <- length(APAlist)
  belownames <- vector()
  for (i in 1:list.len) {
    firstMat <- as.data.frame(APAlist[[i]]$APA); colnames(firstMat) <- 1:size ; rownames(firstMat) <- 1:size
    secondMat <- as.data.frame(APAlist[[focus]]$APA); colnames(secondMat) <- 1:size ; rownames(secondMat) <- 1:size
    a <- reshape2::melt(as.matrix(firstMat))
    a$sample <- factor(rep(names(APAlist)[i], length(a[, 1])))
    abovePlots <- rbind(abovePlots, a)

    e <-  reshape2::melt(as.matrix(firstMat - secondMat))
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


  plot1 <- ggplot2::ggplot(abovePlots, ggplot2::aes(Var1, Var2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::facet_grid(. ~ sample) +
    ggplot2::coord_fixed() +

    GENOVA_THEME() +

    ggplot2::scale_x_continuous(breaks = c(tickPosDownstream,size.banks + 1, tickPosUpstream),
                                expand = c(0,0),
                                labels = c(paste0(tickLabelUpstream,  "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(tickPosDownstream,size.banks + 1, tickPosUpstream),
                                expand = c(0,0),
                                labels = c(paste0(tickLabelDownstream,  "kb"), "5'", paste0(tickLabelUpstream, "kb"))) +
    ggplot2::labs(title = title, x = "", y = "", fill = "Contacts ") +
    ggplot2::scale_fill_gradientn(colours = higlassCol, limits = z)




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
    ggplot2::facet_grid(. ~ sample) +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_gradient2(limits = z2, midpoint = 0, low = "#2166ac", mid = "white", high = "#b2182b") +
    GENOVA_THEME() +
    ggplot2::scale_x_continuous(breaks = c(tickPosDownstream,size.banks + 1, tickPosUpstream),
                                expand = c(0,0),
                                labels = c(paste0(tickLabelDownstream,  "kb"), "3'", paste0(tickLabelDownstream, "kb"))) +
    ggplot2::scale_y_continuous(breaks = c(tickPosDownstream,size.banks + 1, tickPosUpstream),
                                expand = c(0,0),
                                labels = c(paste0(tickLabelDownstream, "kb"), "5'", paste0(tickLabelUpstream, "kb"))) +
    ggplot2::labs(x = "",   y = "", fill = "Difference")

  grid::grid.newpage()
  grid::grid.draw(rbind(ggplot2::ggplotGrob(plot1), ggplot2::ggplotGrob(plot2),
                        size = "last"))
}
