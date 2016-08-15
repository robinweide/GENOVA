#' Plot the virtual 4C results.
#'
#' @param data Output of `virtual.4C`.
#' @return A plot.
#' @export
visualise.virtual.4C <- function(data){
  plot(data, type = 'l', ylim = c(min(data), max(data)*1.1), xlab = 'Normalised distance', ylab = 'Score per loop')
  points(x = 31, y = max(data[25:35])+1,pch = 25, col = 4)
  points(x = 73, y = max(data[70:80])+1,pch = 25, col = 2)
  legend(x = 80, y = max(data)*1.125,legend = c('Viewpoint', 'Anchor'), bty = 'n', col = c(4,2), pch = c(25,25))
}
