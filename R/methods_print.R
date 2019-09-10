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

#' @export
#' @keywords internal
print.APA_discovery <- function(x) {
  res <- abs(median(diff(as.numeric(dimnames(x$signal)[[1]]))))
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  string <- paste0("A", attr(x, "package"), "APA_discovery object involving the following ",
                   dim(x$signal)[3],
                   ' experiment(s):\n"',
                   paste0(dimnames(x$signal)[[3]], collapse = '", "'),
                   '" at a resolution of ', res,'.\n\n')
  slots0 <- paste0("Contains the following slots:\n")
  slots1 <- paste0("- signal:\tAn ", dim(x$signal)[1], " x ", dim(x$signal)[2],
                   " x ", dim(x$signal)[3], " array with summarised APA results.\n")
  slots2 <- if("signal_raw" %in% names(x)) {
    paste0("- signal_raw:\tA list of length ", length(x$signal_raw),
           " containing arrays of unsummarised results.\n")
  } else {""}

  cat(string)
  cat(slots0)
  cat(slots1)
  cat(slots2)
}

#' @export
#' @keywords internal
print.PESCAn_discovery <- function(x) {
  res <- abs(median(diff(as.numeric(dimnames(x$signal)[[1]]))))
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  string <- paste0("A ", attr(x, "package"), " PESCAn_discovery object involving the following ",
                   dim(x$signal)[3], ' experiments(s):\n"',
                   paste0(dimnames(x$signal)[[3]], collapse = '", "'),
                   '" at a resolution of ', res, '.\n\n')
  slots0 <- paste0("Contains the following slots:\n")
  slots1 <- if ("obsexp" %in% names(x)) {
    paste0("- obsexp:\tAn ", dim(x$obsexp)[1], " x ", dim(x$obsexp)[2], " x ",
           dim(x$obsexp)[3], " array with summarised PESCAn results.\n")
  } else {""}
  slots2 <- if ("signal" %in% names(x)) {
    paste0("- signal:\tAn ", dim(x$signal)[1], " x ", dim(x$signal)[2], " x ",
           dim(x$signal)[3], " array with the summarised signal for the anchors.\n")
  } else {""}
  slots3 <- if ("shifted" %in% names(x)) {
    paste0("- shifted:\tAn ", dim(x$shifted)[1], " x ", dim(x$shifted)[2], " x ",
           dim(x$shifted)[3], " array with the summarised signal for shifted anchors.\n")
  }
  slots4 <- if ("signal_raw" %in% names(x)) {
    paste0("- signal_raw:\tA list of length ", length(x$signal_raw),
           " containing arrays of unsummarised results for the anchors.\n")
  } else {""}
  slots5 <- if ("shifted_raw" %in% names(x)) {
    paste0("- shifted_raw:\tA list of length ", length(x$shifted_raw),
           " containing arrays of unsummarised results for the shifted anchors.\n")
  }

  cat(string)
  cat(slots0)
  cat(slots1)
  cat(slots2)
  cat(slots3)
  cat(slots4)
  cat(slots5)
}
