#' Plot the virtual 4C results.
#'
#' @param data Output of `virtual.4C`.
#' @return A plot.
visualise.virtual.4C <- function(data){
  plot(data, type = 'l', ylim = c(min(data), max(data)*1.1), xlab = 'Normalised distance', ylab = 'Average score per loop')
  points(x = 30, y = max(data[25:35])*1.05,pch = 25, col = 4)
  points(x = 71, y = max(data[70:80])+(max(data[25:35])*.05),pch = 25, col = 2)
  legend(x = 80, y = max(data)*1.125,legend = c('Viewpoint', 'Anchor'), bty = 'n', col = c(4,2), pch = c(25,25))
}
