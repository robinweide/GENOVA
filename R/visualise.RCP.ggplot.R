#' Plot the RCP results
#'
#' Produces a line-plot per chromosome.
#'
#' @param RCPdata Output of `RCP`
#' @return A ggplot-object
#' @export
visualise.RCP.ggplot <-function(RCPdata){
  cols <- RCPdata$color
  names(cols) <- RCPdata$sample

  ggplot2::ggplot(RCPdata, ggplot2::aes(col = sample,x = distance, y = prob)) +
    #geom_point(size = .25,alpha = .25) +
    ggplot2::geom_smooth(span = 0.25, se = F) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap( ~ chrom, ncol = floor(sqrt( length(unique(RCP_WTWAPL$chrom))  )) ) +
    ggplot2::theme_linedraw()+
    ggplot2::scale_color_manual(values = cols)+
    ggplot2::labs(title = "Relative contact probability", x = "Distance (bp)", y = "RCP", col = '')+
    ggplot2::coord_fixed(ratio = .75)
}
