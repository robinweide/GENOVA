#' Plot the RCP results
#'
#' Produces a line-plot per chromosome and (optionally) per-BED.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param RCPdata Output of the RCP-function.
#' @param smooth Plot a lowess-smoothed line?
#' @param combine Combine all chromosomes?
#' @param ylim An optional vector of two values, containing the y-axis limits.
#' @param xlim An optional vector of two values in Mb, containing the x-axis limits.
#' @param lineOnly Set to FALSE if you want to plot points and error-bars.
#' @param lineWidth Width of the RCP-lines
#' @param errWidth Width of SEM-lines
#' @param pointWidth Width of points
#' @examples
#' # Run RCP on the first five chromosomes
#' RCP_out <- RCP(list(WT,Mut1), chromsToUse = paste0('chr',1:5))
#'
#' # Plot the RCP-output between 5 and 100Mb per chromosome.
#' visualise.RCP.ggplot(RCP_out,combine = F,smooth = F, lineOnly = F, xlim = c(5,110))
#'
#' # Plot the RCP-output with a loess-smoothed line without the actual observed points
#' # between 5 and 100Mb and combine all chromosomes in one plot.
#' visualise.RCP.ggplot(RCP_out,combine = T,smooth = T, lineOnly = T, xlim = c(5,110))
#' @return A ggplot-object
#' @export
visualise.RCP.ggplot <-function(RCPdata, smooth =F, combine = T, ylim = NULL, xlim = NULL,lineOnly = T, lineWidth = 1, errWidth = 1, pointWidth = 0.5){



  cols <- RCPdata$color
  p <- NULL
  colorOnBed = F
  if(nlevels(RCPdata$sample) > 1){
    names(cols) <- RCPdata$sample
  } else if(nlevels(RCPdata$BED) > 1){
    names(cols) <- RCPdata$BED
    colorOnBed = T
  } else {
    names(cols) <- RCPdata$sample
  }

  if (combine == T) {
    RCPdata <- dplyr::group_by(RCPdata, distance, sample, BED)
    RCPdata <- dplyr::summarise(RCPdata, prob = median(prob), SEM = sum(SEM), col = unique(color))
    RCPdata$chrom <- "All chromosomes"
  }

  RCPdata$BED = factor(RCPdata$BED)
  if(nlevels(RCPdata$BED) == 0){
    RCPdata$BED = factor(1)
  }


  if (smooth == F) {

    if(colorOnBed){
      p <- ggplot2::ggplot(RCPdata, ggplot2::aes(col = BED, x = distance/1000000,
                                                 y = prob)) +
        ggplot2::geom_line(size = lineWidth, mapping = ggplot2::aes( lty = BED)) +
        ggplot2::facet_wrap(~chrom, nrow = floor(sqrt(length(unique(RCPdata$chrom))))) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::labs(title = "Relative contact probability",
                      x = "Distance (Mbp)", y = "RCP", col = "", lty = "") +
        GENOVA_THEME() +
        ggplot2::theme(aspect.ratio = 1) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5, hjust = 1))
    } else {

      p <- ggplot2::ggplot(RCPdata, ggplot2::aes(col = sample, x = distance/1000000,
                                                 y = prob)) +
        ggplot2::geom_line(size = lineWidth, mapping = ggplot2::aes( lty = BED)) +
        ggplot2::facet_wrap(~chrom, nrow = floor(sqrt(length(unique(RCPdata$chrom))))) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::labs(title = "Relative contact probability",
                      x = "Distance (Mbp)", y = "RCP", col = "", lty = "") +
        ggplot2::theme(aspect.ratio = 1) +
        GENOVA_THEME() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5, hjust = 1))


    }

    if(lineOnly == F){
      p <- p + ggplot2::geom_pointrange(size = errWidth,
                                        shape = 20,
                                        fatten = pointWidth,
                                        mapping =  ggplot2::aes(col = sample,
                                                                ymin = prob-SEM,
                                                                ymax = prob+SEM))
    }





    if(!is.null(ylim)){
      p <- p + ggplot2::scale_y_continuous(limits = ylim, trans = 'log10')
    } else {
      p <- p + ggplot2::scale_y_log10()
    }

    if(!is.null(xlim)){
      p <- p + ggplot2::scale_x_continuous(limits = xlim, trans = 'log10')
    } else {
      p <- p + ggplot2::scale_x_log10()
    }



  }

  ##########################
  #    smooth    #
  ##########################
  else {
    if(colorOnBed){

    p <- ggplot2::ggplot(RCPdata, ggplot2::aes(col = BED, x = distance/1000000,
                                               y = prob)) +
      ggplot2::geom_smooth(mapping = ggplot2::aes( lty = BED),size = lineWidth,  span = 0.25, se = F) +
      ggplot2::facet_wrap(~chrom, nrow = floor(sqrt(length(unique(RCPdata$chrom))))) +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::labs(title = "Relative contact probability",
                    x = "Distance (Mbp)", y = "RCP", col = "", lty = "") +
      GENOVA_THEME() +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         vjust = 0.5, hjust = 1))
    } else {
      p <- ggplot2::ggplot(RCPdata, ggplot2::aes(col = sample, x = distance/1000000,
                                                 y = prob)) +
        ggplot2::geom_smooth(mapping = ggplot2::aes( lty = BED),size = lineWidth,  span = 0.25, se = F) +
        ggplot2::facet_wrap(~chrom, nrow = floor(sqrt(length(unique(RCPdata$chrom))))) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::labs(title = "Relative contact probability",
                      x = "Distance (Mbp)", y = "RCP", col = "", lty = "") +
        GENOVA_THEME() +
        ggplot2::theme(aspect.ratio = 1) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5, hjust = 1))
    }

    if(lineOnly == F){
      p <- p + ggplot2::geom_pointrange(size = errWidth, shape = 20, fatten = pointWidth, ggplot2::aes(col = sample,ymin = prob-SEM,ymax = prob+SEM))
    }

    if(!is.null(ylim)){
      p <- p + ggplot2::scale_y_continuous(limits = ylim, trans = 'log10')
    } else {
      p <- p + ggplot2::scale_y_log10()
    }

    if(!is.null(xlim)){
      p <- p + ggplot2::scale_x_continuous(limits = xlim, trans = 'log10')
    } else {
      p <- p + ggplot2::scale_x_log10()
    }

  }


  # if all beds are the same (or NA), remove guide altogether
  if(nlevels(RCPdata$BED) < 2){
    p = p+ggplot2::guides(lty=FALSE)
  }
  p = p+ggplot2::guides(size=FALSE)
  suppressMessages(p + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                                      panel.grid.major =  ggplot2::element_line(colour = '#333333')))
}

GENOVA_THEME = function(){
  require(ggplot2)
  p = ggplot2::theme(panel.background = ggplot2::element_blank(),
                     legend.key =  ggplot2::element_rect(fill = 'white'),
                     strip.background = ggplot2::element_rect(fill = NA, colour = NA),
                     panel.border = ggplot2::element_rect(fill = NA, colour = 'black'),
                     text = ggplot2::element_text(color = 'black'),
                     axis.text = ggplot2::element_text(colour = 'black'),
                     strip.text = ggplot2::element_text(colour = 'black') )
  return(p)
}
