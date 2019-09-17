# Contacts class ----------------------------------------------------------

#' @export
#' @keywords internal
print.contacts <- function(x) {
  res <- x$RES
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  string <- paste0("A ", attr(x, "package"), " ", class(x), " object named '",
                   x$NAME, "' with the following slots:\n")
  slots <- sapply(x, is.null)
  slots <- names(slots[!slots])

  ice <- paste0("- ICE\t:\tTriplet format matrix containing ", dim(x$ICE)[1],
                " informative bins.\n")
  abs <- paste0("- ABS\t:\t", dim(x$ABS)[1], " genomic indices in BED format.\n")
  cname <- paste0('- NAME\t:\t"', x$NAME, '"\n')
  pres <- paste0("- RES\t:\tResolution at ", res, ".\n")
  chroms <- paste0("- CHRS\t:\tA vector of ", length(x$CHRS), " chromosome names.\n")
  col <- paste0("- COL\t:\tThe colour '", x$COL, "'\n")
  mask <- paste0("- MASK\t:\t", sum(x$MASK), " bins are masked.\n")
  rmchrom <- paste0("- RMCHROM:\t", if(x$RMCHROM){"Some"} else {"No"},
                    " chromosomes have been removed.\n")
  zscore <- paste0("- ZSCORE:\tThe data have ", if (!x$ZSCORE) {"not"},
                   " been Z-score normalised.\n")
  comment <- paste0("- COMM\t:\tThe following comment: ",
                    if (!is.null(x$COMM)) {x$COMM} else {"NULL"},
                    "\n")
  centros <- if (is.null(x$CENTROMERES)) {
    paste0("- CENTROMERES:\tNo centromere information.\n")
  } else {
    paste0("- CENTROMERES:\tLocations of ", nrow(x$CENTROMERES), " centromeres.\n")
  }

  cat(string)
  cat(ice)
  cat(abs)
  cat(cname)
  cat(pres)
  cat(chroms)
  cat(col)
  cat(mask)
  cat(rmchrom)
  cat(zscore)
  cat(centros)
  cat(comment)
}


# Discovery classes -------------------------------------------------------

print.ARMLA_discovery <- function(x) {
  res <- attr(x, "resolution")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }

  myclass <- class(x)[[1]]
  miniclass <- strsplit(myclass, "_")[[1]][[1]]
  aan <- if(startsWith(myclass, "A")) {
    "An"
  } else {
    "A"
  }

  opening <- paste0(
    aan, " ", attr(x, "package"), " '", myclass, "' object involving the ",
    "following ", dim(x$signal)[3],' experiments:\n"',
    paste0(dimnames(x$signal)[[3]], collapse = '", "'),
    '" at a resolution of ', res, '.\n\n'
  )

  slots0 <- paste0("Contains the following slots:\n")
  slots1 <- paste0("- signal:\tAn ", paste0(dim(x$signal), collapse = " x "),
                   " array with summarised ", miniclass, " results at",
                   " anchor positions.\n")
  slots2 <- if ("obsexp" %in% names(x)) {
    paste0("- obsexp:\tAn ", paste0(dim(x$obsexp), collapse = " x "),
           " array with summarised and normalised ", miniclass, " results.\n")
  } else ""
  slots3 <- if ("shifted" %in% names(x)) {
    paste0("- shifted:\tAn ", paste0(dim(x$shifted), collapse = " x "),
           " array with summarised ", miniclass, " results at shifted ",
           "anchor positions.\n")
  } else ""
  slots4 <- if ("signal_raw" %in% names(x)) {
    paste0("- signal_raw:\tA list of length ", length(x$signal_raw),
          " containing raw results at each anchor position.\n")
  } else ""
  slots5 <- if("shifted_raw" %in% names(x)) {
    paste0("- shifted_raw:\tA list of length ", length(x$shifted_raw),
           " containing raw results at each shifted anchor position.\n")
  } else ""

  cat(opening)
  cat(slots0)
  cat(slots1)
  cat(slots2)
  cat(slots3)
  cat(slots4)
  cat(slots5)
}

# Other classes -----------------------------------------------------------

#' @export
#' @keywords internal
print.anchors <- function(x) {
  str1 <- paste0("An 'anchors' object of type '", attr(x, "type"),
                 "' of length ", nrow(x),":\n")
  class(x) <- "matrix"
  x <- x[T,]
  cat(str1)
  p <- print(head(x))
  outp <<- p
  n <- max(nchar(p))
  if (is.null(dimnames(p))) {
    if (nrow(x) > nrow(p)) {
      space <- "."
      str2 <- paste0("[", space, ",] ",
                     paste0(rep(".", n), collapse = ""), " ",
                     paste0(rep(".", n), collapse = ""))
      cat(str2)
    }
  } else {
    m <- max(nchar(dimnames(p)[[1]]))
    val <- tail(dimnames(p)[[1]], 1)
    val <- as.numeric(substr(val, 2, m - 2))
    if (nrow(x) > val) {
      space <- paste0(rep(".", m - 3), collapse = "")
      str2 <- paste0("[", space, ",] ",
                     paste0(rep(".", n), collapse = ""), " ",
                     paste0(rep(".", n), collapse = ""))
      cat(str2)
    }
  }
}

