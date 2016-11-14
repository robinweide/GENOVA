#' Plot the RCP results
#'
#' Produces a line-plot per chromosome.
#'
#' @param RCPdata Output of `RCP`
#' @param smooth Plot a lowess-smoothed line?
#' @param combine Combine all chromosomes?
#' @return A ggplot-object
#' @export
visualise.RCP.ggplot <-function(RCPdata, smooth =F, combine = T){
  cols <- RCPdata$color
  names(cols) <- RCPdata$sample
  
  if(combine == T){
    RCPdata <- group_by(RCPdata, distance, sample) 
    RCPdata <- summarise(RCPdata, prob = mean(prob), col = unique(color))
    RCPdata$chrom <- "All chromosomes"
  }
  
  if(smooth == F){
    ggplot2::ggplot(RCPdata, ggplot2::aes(col = sample, x = distance/1e6, 
                                          y = prob)) + ggplot2::geom_line() + 
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + 
      ggplot2::facet_wrap(~chrom, nrow = floor(sqrt(length(unique(RCPdata$chrom))))) + 
      ggplot2::theme_linedraw() + ggplot2::scale_color_manual(values = cols) + 
      ggplot2::labs(title = "Relative contact probability", 
                    x = "Distance (Mbp)", y = "RCP", col = "") + ggplot2::theme(aspect.ratio = 1) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5, hjust=1))
  } else {
    ggplot2::ggplot(RCPdata, ggplot2::aes(col = sample, x = distance/1e6, 
                                          y = prob)) + ggplot2::geom_smooth(span = 0.25, se = F)  + 
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + 
      ggplot2::facet_wrap(~chrom, nrow = floor(sqrt(length(unique(RCPdata$chrom))))) + 
      ggplot2::theme_linedraw() + ggplot2::scale_color_manual(values = cols) + 
      ggplot2::labs(title = "Relative contact probability", 
                    x = "Distance (Mbp)", y = "RCP", col = "") + ggplot2::theme(aspect.ratio = 1) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,vjust = 0.5, hjust=1))
  }
}
