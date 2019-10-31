# Contacts class ----------------------------------------------------------

#' @export
#' @keywords internal
print.contacts <- function(x, ...) {
  res <- attr(x, "res")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  string <- paste0("A ", attr(x, "package"), " ", class(x), " object named '",
                   attr(x, "samplename"), "' at a resolution of ", res, 
                   ".\nContains the following slots:\n")
  slots <- sapply(x, is.null)
  slots <- names(slots[!slots])

  mat <- paste0("- MAT\t:\tTriplet format matrix containing ", dim(x$MAT)[1],
                " informative bins.\n")
  idx <- paste0("- IDX\t:\t", dim(x$IDX)[1], " genomic indices in BED format.\n")
  chroms <- paste0("- CHRS\t:\tA vector of ", length(x$CHRS), " chromosome names.\n")
  
  
  col <- paste0("This object is assigned the colour '", attr(x, "col"), "'\n")
  mask <- paste0(sum(x$MASK), " bins are masked.\n")
  rmchrom <- paste0(if(attr(x, "rmChrom")){"Some"} else {"No"},
                    " chromosomes have been removed.\n")
  zscore <- paste0("The data have ", if (!attr(x, "znorm")) {"not"},
                   " been Z-score normalised.\n")
  centros <- if (is.null(x$CENTROMERES)) {
    paste0("- CENTROMERES:\tNo centromere information.\n")
  } else {
    paste0("- CENTROMERES:\tLocations of ", nrow(x$CENTROMERES), " centromeres.\n")
  }
  bal <- paste0("The orgininal data were loaded in as ", if (attr(x, "balanced")) {
    "balanced"
  } else { "raw" }, " data.")

  cat(string)
  cat(mat)
  cat(idx)
  cat(chroms)
  cat(centros)
  cat(col)
  cat(mask)
  cat(rmchrom)
  cat(zscore)
  cat(bal)
}


# Discovery classes -------------------------------------------------------

#' @export
#' @keywords internal
print.ARMLA_discovery <- function(x, ...) {
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

  opening <- paste0(
    "A ", attr(x, "package"), " '", myclass, "' object involving the ",
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

#' @export
#' @keywords internal
print.RCP_discovery <- function(x, ...) {
  
  string <- paste0("A ", attr(x, "package"), " ", 
                   'RCP_discovery',  " object with the following details:\n")
  
  smpls = unique(x$raw$samplename)
  smpls = paste(paste0(smpls[-length(smpls)], collapse = ', '), 
                smpls[length(smpls)], sep = ' & ')
  print_samples = paste0("- samples: ", smpls,"\n")
  
  regs = unique(x$raw$region)
  regs = paste(paste0(regs[-length(regs)], collapse = ', '), 
               regs[length(regs)], sep = ' & ')
  print_regions = paste0("- regions: ",regs,"\n")
  
  print_norm = paste0("- normalisation: ", attr(x, 'norm'),"\n")
  
  cat(string)
  cat(print_samples)
  cat(print_regions)
  cat(print_norm)
  
}

#' @export
#' @keywords internal
print.CS_discovery <- function(x, ...) {
  
  myclass <- class(x)[1]
  cols <- colnames(x$compart_scores)
  res <- attr(x, "resolution")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  part <- attr(x, "partitioning")
  
  string <- paste0("A ", attr(x, "PACKAGE"), " '", myclass, 
                   "' object involving the following ", 
                   length(cols) - 4, " experiments:\n'", 
                   paste0(cols[5:length(cols)], collapse = "', '"), "' at a ",
                   "resolution of ", res, ".\n")
  slot0 <- "Contains the following slots:\n\n"
  
  slot1 <- if (attr(x, "signed")) {
    paste0("- compart_scores:\tA data.frame containing ",
           sum(!is.na(x$compart_scores[[5]])),
           " compartment scores.\n\n")
  } else {
    paste0("- compart_scores:\tA data.frame containing ",
           sum(!is.na(x$compart_scores[[5]])),
           " scores.\n\n")
  }

  ncentro <- sum(grepl("centro$", part$values))
  centrobins <- sum(part$lengths[grepl("centro$", part$values)])
  
  details1 <- paste0(ncentro, " centromeres spanning a total of ", centrobins, 
                     " bins have been ignored.\n")
  signed <- paste0("The scores are signed.")
  unsigned <- paste0("The scores are unsigned.")
  
  cat(string)
  cat(slot0)
  cat(slot1)
  cat(details1)
  if (attr(x, "signed")) {
    cat(signed)
  } else {
    cat(unsigned)
  }
}

#' @export
#' @keywords internal
print.IS_discovery <- function(x, ...) {
  
  myclass <- class(x)[1]
  cols <- colnames(x$insula_score)
  res <- attr(x, "resolution")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  part <- attr(x, "partitioning")
  
  string <- paste0("A ", attr(x, "PACKAGE"), " '", myclass, 
                   "' object involving the following ", 
                   length(cols) - 4, " experiments:\n'", 
                   paste0(cols[5:length(cols)], collapse = "', '"), "' at a ",
                   "resolution of ", res, ".\n")
  slot0 <- "Contains the following slots:\n\n"
  
  slot1 <- paste0("- insula_score:\tA data.frame containing ",
                  sum(!is.na(x$insula_score[[5]])),
                  " insulation scores.\n\n")
  
  details1 <- paste0("The scores have been called using a ", attr(x, "window"),
                     " x ", attr(x, "window"), " sliding square.")

  
  cat(string)
  cat(slot0)
  cat(slot1)
  cat(details1)
}

#' @export
#' @keywords internal
print.saddle_discovery <- function(x, ...) {
  myclass <- class(x)
  n_bins <- max(c(x$saddle$q1, x$saddle$q2), na.rm = TRUE)
  expnames <- unique(x$saddle$exp)
  
  res <- attr(x, "resolution")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  
  n_arms <- length(unique(x$saddle$chr))
  
  string <- paste0("A ", attr(x, "package"), " '", myclass, 
                   "' object involving the following ", 
                   length(expnames), " experiments:\n'", 
                   paste0(expnames, collapse = "', '"), "' at a ",
                   "resolution of ", res, ".\n")
  slot0 <- "Contains the following slots:\n\n"
  slot1 <- paste0("- saddle:\tA data.frame containing quantile-quantile ",
                  "scores for ", n_arms, " chromosome arms.")
  
  cat(string)
  cat(slot0)
  cat(slot1)
}

#' @export
#' @keywords internal
print.domainogram_discovery <- function(x, ...) {
  myclass <- class(x)[1]
  pos <- format(range(x$position), scientific = FALSE)
  chrom <- attr(x, "chr")
  expnames <- unique(x$experiment)
  res <- attr(x, "resolution")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  
  string <- paste0("A ", attr(x, "package"), " '", myclass,
                   "' object involving the following ",
                   length(expnames), " experiments:\n'",
                   paste0(expnames, collapse = "', '"), "' at a ",
                   "resolution of ", res, ".\nPosition ",
                   chrom, ":", pos[1], "-", pos[2]," spanning ",
                  diff(range(x$window)), " window sizes.\n")
  cat(string)
  print(as.data.table(x))
}

#' @export
#' @keywords internal
print.virtual4C_discovery <- function(x, ...) {
  myclass <- class(x)[[1]]
  expnames <- unique(x$data$experiment)
  res <- attr(x, "resolution")
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  
  vp <- attr(x, "viewpoint")

  string <- paste0("A ", attr(x, "package"), " '", myclass,
                   "' object involving the following ",
                   length(expnames), " experiments:\n'",
                   paste0(expnames, collapse = "', '"), "' at a ",
                   "resolution of ", res, ".\n")
  string1 <- paste0("The viewpoint of this virtual 4C is located at ",
                    vp[1, 1], ":", format(vp[1, 2], scientific = FALSE), "-",
                    format(vp[1, 3], scientific = FALSE), ".")
  
  cat(string)
  cat(string1)
}

# Other classes -----------------------------------------------------------

#' @export
#' @keywords internal
print.anchors <- function(x, ...) {
  str1 <- paste0("An 'anchors' object of type '", attr(x, "type"),
                 "' of length ", nrow(x),":\n")
  class(x) <- "matrix"
  x <- x[T,]
  cat(str1)
  p <- print(head(x))
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

