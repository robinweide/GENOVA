#' @export
#' @keywords internal
persp.PESCAn_discovery <- function(
  x, ...
) {
  mat <- x$obsexp

  par(mfrow = seq_len(dim(mat)[3]))
  
  # Set default arguments
  args <- pairlist(
    x = rev(as.numeric(dimnames(mat)[[1]])/1000),
    y = as.numeric(dimnames(mat)[[2]])/1000,
    z = mat,
    ylab = "5' (kb)", xlab = "3' (kb)",
    zlab = "Observed / Expected", col = "dodgerblue", ticktype = "detailed",
    theta = 30, phi = 25, shade = 0.75, ltheta = -15,
    ... = NULL
  )
  args <- c(
    args,
    list(
      xlim = range(args$x),
      ylim = range(args$y),
      zlim = range(args$z, na.rm = TRUE)
    )
  )
  
  # Replace default arguments with ellipsis arguments
  dot_args <- list(...)
  
  xtra <- setdiff(names(args), names(dot_args))
  
  args <- c(dot_args, args[xtra])
  
  # Fill rest of argument list with the function's formals
  fun <- getFromNamespace("persp.default", "graphics")
  
  fun_formals <- formals(fun)
  
  xtra <- setdiff(names(fun_formals), names(args))
  
  args <- c(args, fun_formals[xtra])
  
  cols <- args$col
  
  # Call function with arguments
  for (i in seq_len(dim(mat)[3])) {
    m <- mat[, , i]
    m[cbind(rev(row(m)[T]), col(m)[T])] <- m
    args$z <- m
    args$col <- cols[[min(i, length(cols))]]
    do.call(fun, args)
  }
}

#' @export
#' @keywords internal
persp.APA_discovery <- function(
  x, ...
) {
  mat <- x$signal
  
  par(mfrow = seq_len(dim(mat)[3]))
  
  # Set default arguments
  args <- pairlist(
    x = rev(as.numeric(dimnames(mat)[[1]])/1000),
    y = as.numeric(dimnames(mat)[[2]])/1000,
    z = mat,
    ylab = "5' (kb)", xlab = "3' (kb)",
    zlab = "Mean Contacts", col = "dodgerblue", ticktype = "detailed",
    theta = 30, phi = 25, shade = 0.75, ltheta = -15,
    ... = NULL
  )
  args <- c(
    args,
    list(
      xlim = range(args$x),
      ylim = range(args$y),
      zlim = range(args$z, na.rm = TRUE)
    )
  )
  
  # Replace default arguments with ellipsis arguments
  dot_args <- list(...)
  
  xtra <- setdiff(names(args), names(dot_args))
  
  args <- c(dot_args, args[xtra])
  
  # Fill rest of argument list with the function's formals
  fun <- getFromNamespace("persp.default", "graphics")
  
  fun_formals <- formals(fun)
  
  xtra <- setdiff(names(fun_formals), names(args))
  
  args <- c(args, fun_formals[xtra])
  
  cols <- args$col
  
  # Call function with arguments
  for (i in seq_len(dim(mat)[3])) {
    m <- mat[, , i]
    m[cbind(rev(row(m)[T]), col(m)[T])] <- m
    args$z <- m
    args$col <- cols[[min(i, length(cols))]]
    do.call(fun, args)
  }
}
