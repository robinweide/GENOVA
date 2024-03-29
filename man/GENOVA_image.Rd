% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_image.R
\name{GENOVA_image}
\alias{GENOVA_image}
\alias{image.contacts_matrix}
\title{Display a part of a Hi-C matrix}
\usage{
\method{image}{contacts_matrix}(x, col = NULL, zlim = NULL, xlab = NULL, ylab = NULL, ...)
}
\arguments{
\item{x}{A \code{contacts_matrix} object, as returned by
\code{\link[GENOVA]{select_subset}()}.}

\item{col}{a list of colors such as that generated by
\code{\link[grDevices]{hcl.colors}}, \code{\link[grDevices]{gray.colors}}
or similar functions. Defaults to GENOVA palette.}

\item{zlim}{the minimum and maximum z values for which colors should be
plotted, defaulting to the range of finite values of off-diagonal elements
of \code{x}.}

\item{xlab, ylab}{each a character string giving the labels for the x and y
axis. Default is to use the chromosome name as axis title.}

\item{...}{
  Arguments passed on to \code{\link[graphics:image]{graphics::image.default}}
  \describe{
    \item{\code{xlim,ylim}}{ranges for the plotted \code{x} and \code{y} values,
    defaulting to the ranges of \code{x} and \code{y}.}
    \item{\code{add}}{logical; if \code{TRUE}, add to current plot (and disregard
    the following four arguments).  This is rarely useful because
    \code{image} \sQuote{paints} over existing graphics.}
    \item{\code{xaxs,yaxs}}{style of x and y axis.  The default \code{"i"} is
    appropriate for images.  See \code{\link[graphics]{par}}.}
    \item{\code{xlab,ylab}}{each a character string giving the labels for the x and
    y axis.  Default to the \sQuote{call names} of \code{x} or \code{y}, or to
    \code{""} if these were unspecified.}
    \item{\code{breaks}}{a set of finite numeric breakpoints for the colours:
    must have one more breakpoint than colour and be in increasing
    order.  Unsorted vectors will be sorted, with a warning.}
    \item{\code{oldstyle}}{logical. If true the midpoints of the colour intervals
    are equally spaced, and \code{zlim[1]} and \code{zlim[2]} were taken
    to be midpoints.  The default is to have colour intervals of equal
    lengths between the limits.}
    \item{\code{useRaster}}{logical; if \code{TRUE} a bitmap raster is used to
    plot the image instead of polygons. The grid must be regular in that
    case, otherwise an error is raised.   For the behaviour when this is
    not specified, see \sQuote{Details}.}
  }}
}
\description{
Display a part of a Hi-C matrix
}
\examples{
\dontrun{
x <- select_subset(exp, "chr1", 50e6, 60e6)
image(x)
}
}
