#' fields.interp.surface
#'
#' fields.interp.surface
#'
#' @param obj A surface or image  object like the list for contour, persp or image.
#' @param lof A matrix of 2 d locations -- new points to evaluate the surface.
#' @note
#' THIS FUNCTION IS TAKEN FROM FIELDS!
#'
#' fields is a package for analysis of spatial data written for
#' the R software environment .
#' Copyright (C) 2016
#' University Corporation for Atmospheric Research (UCAR)
#' Contact: Douglas Nychka, nychka@ucar.edu,
#' National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#'
#' This program is free software; you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation; either version 2 of the License, or
#' (at your option) any later version.
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#'
#' You should have received a copy of the GNU General Public License
#' along with the R software environment if not, write to the Free Software
#' Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#' or see http://www.r-project.org/Licenses/GPL-2
#' @return An interpolated matrix
#' @export
fields.interp.surface <- function(obj, loc) {

    # obj is a surface or image  object like the list for contour, persp or image.
    # loc a matrix of 2 d locations -- new points to evaluate the surface.
    x <- obj$x
    y <- obj$y
    z <- obj$z
    nx <- length(x)
    ny <- length(y)
    # this clever idea for finding the intermediate coordinates at the new points
    # is from J-O Irisson
    lx <- approx(x, 1:nx, loc[, 1])$y
    ly <- approx(y, 1:ny, loc[, 2])$y
    lx1 <- floor(lx)
    ly1 <- floor(ly)
    # x and y distances between each new point and the closest grid point in the lower left hand corner.
    ex <- lx - lx1
    ey <- ly - ly1
    # fix up weights to handle the case when loc are equal to
    # last grid point.  These have been set to NA above.
    ex[lx1 == nx] <- 1
    ey[ly1 == ny] <- 1
    lx1[lx1 == nx] <- nx - 1
    ly1[ly1 == ny] <- ny - 1
    # bilinear interpolation finds simple weights based on the
    # the four corners of the grid box containing the new
    # points.
    return(z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + z[cbind(lx1 +
        1, ly1)] * ex * (1 - ey) + z[cbind(lx1, ly1 + 1)] * (1 -
        ex) * ey + z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}
