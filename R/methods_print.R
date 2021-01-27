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
  centros <- if (is.null(x$CENTROMERES) || all(x$CENTROMERES$start == -1)) {
    paste0("- CENTROMERES:\tNo centromere information.\n")
  } else {
    paste0("- CENTROMERES:\tLocations of ", nrow(x$CENTROMERES), " centromeres.\n")
  }
  bal <- paste0("The original data were loaded in as ", if (attr(x, "balanced")) {
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

  myclass <- class(x)[[1]]
  miniclass <- strsplit(myclass, "_")[[1]][[1]]

  opening <- describe(x)
  
  if (myclass == "CSCAn_discovery" && length(dim(x$signal) == 4L)) {
    opening <- substr(opening, 1, nchar(opening) - 5) # strip newlines
    opening <- paste0(opening, '\nThe following groupings apply: "', 
                      paste0(dimnames(x$signal)[[3]], collapse = '", "'), 
                      '".\n\n')
  }
  
  descriptions <- c(
    signal = paste0("with summarised ", miniclass, " results at anchor", 
                    " positions"),
    obsexp = paste0("with summarised and normalised ", miniclass, " results"),
    shifted = paste0("with summarised ", miniclass, " results at shifted",
                     " anchor positions"),
    signal_raw = "containing raw results for each sample at anchor positions",
    shifted_raw = paste0("containing raw results for each sample at shifted",
                         " anchor positions")
  )

  slots0 <- paste0("Contains the following slots:\n")
  slots1 <- describe(x$signal, descriptions["signal"])
  slots2 <- if ("obsexp" %in% names(x)) {
    describe(x$obsexp, descriptions["obsexp"])
  } else ""
  slots3 <- if ("shifted" %in% names(x)) {
    describe(x$shifted, descriptions["shifted"])
  } else ""
  slots4 <- if ("signal_raw" %in% names(x)) {
    describe(x$signal_raw, descriptions["signal_raw"])
  } else ""
  slots5 <- if("shifted_raw" %in% names(x)) {
    describe(x$shifted_raw, descriptions["shifted_raw"])
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

  string <- describe(x)
  
  
  slot0 <- "Contains the following slots:\n"
  
  slot1 <- describe(x$raw, "with raw relative contact probabilities")
  slot2 <- describe(x$smooth, "containing smoothened RCPs")

  regs <- toString(unique(x$raw$region))

  print_regions = paste0("Computed on the following regions: ", regs, ". ")
  
  print_norm = paste0("Normalised in a ", attr(x, 'norm')," manner.\n")
  
  cat(string)
  cat(slot0)
  cat(slot1)
  cat(slot2)
  cat("\n")
  cat(print_regions)
  cat(print_norm)
  
}

#' @export
#' @keywords internal
print.CS_discovery <- function(x, ...) {
  
  string <- describe(x)
  
  slot0 <- "Contains the following slots:\n"
  
  slot1 <- if (attr(x, "signed")) {
    describe(x$compart_scores, "containing signed compartment scores")
  } else {
    describe(x$compart_scores, "containing unsigned compartment scores")
  }
  
  part <- attr(x, "partitioning")
  ncentro <- sum(grepl("centro$", part$values))
  centrobins <- sum(part$lengths[grepl("centro$", part$values)])
  
  details1 <- paste0("\n", ncentro, " centromeres spanning a total of ", 
                     centrobins, " bins have been ignored.\n")
  
  cat(string)
  cat(slot0)
  cat(slot1)
  cat(details1)
}

#' @export
#' @keywords internal
print.IS_discovery <- function(x, ...) {
  
  string <- describe(x)
  slot0 <- "Contains the following slots:\n"
  
  slot1 <- describe(x$insula_score, "containing insulation scores")

  details1 <- paste0("\nThe scores have been called using a ", attr(x, "window"),
                     " x ", attr(x, "window"), " sliding square.")

  cat(string)
  cat(slot0)
  cat(slot1)
  cat(details1)
}

#' @export
#' @keywords internal
print.saddle_discovery <- function(x, ...) {
  n_arms <- length(unique(x$saddle$chr))
  
  string <- describe(x)
  slot0 <- "Contains the following slots:\n"
  slot1 <- describe(x$saddle, paste0("containing quantile-quantile scores for ",
                                     n_arms, " chromosome arms"))
  
  cat(string)
  cat(slot0)
  cat(slot1)
}

#' @export
#' @keywords internal
print.domainogram_discovery <- function(x, ...) {

  pos <- format(range(x$scores$position), scientific = FALSE, trim = TRUE)
  chrom <- attr(x, "chr")
  locus <- paste0(chrom, ":", paste0(pos, collapse = "-"))
  wrange <- range(x$scores$window)
  locus <- paste0("Spans the locus ", locus, " at a window range from ",
                  wrange[1], " to ", wrange[2], ".\n")
  
  slot0 <- "Contains the following slot:\n"
  slot1 <- describe(x$scores, "with insulation scores at various window sizes")

  string <- describe(x)
  cat(string)
  cat(slot0)
  cat(slot1)
  cat("\n")
  cat(locus)
}

#' @export
#' @keywords internal
print.DI_discovery <- function(x, ...) {
  string <- describe(x)
  cat(string)
  print(as.data.table(x$DI))
}

#' @export
#' @keywords internal
print.virtual4C_discovery <- function(x, ...) {
  vp <- attr(x, "viewpoint")

  string <- describe(x)
  string1 <- paste0("The viewpoint of this virtual 4C is located at ",
                    vp[1, 1], ":", format(vp[1, 2], scientific = FALSE), "-",
                    format(vp[1, 3], scientific = FALSE), ".")
  
  cat(string)
  cat(string1)
}

#' @export
#' @keywords internal
print.IIT_discovery <- function(x, ...) {
  string <- describe(x)
  cat(string)
  
  slots0 <- paste0("Contains the following slots:\n")
  slots1 <- describe(x$results, "with TAD IDs and scores")
  slots2 <- describe(x$tads, "with TAD positions and IDs")
  
  cat(slots0)
  cat(slots1)
  cat(slots2)
}

#' @export
#' @keywords internal
print.chrommat_discovery <- function(x, ...) {
  string <- describe(x)
  cat(string)
  slots0 <- paste0("Contains the following slots:\n")
  slots1 <- describe(x$obs, "containing summed contacts")
  slots2 <- describe(x$exp, paste0(
    "containing expected proportions,\n\tcalculated with the '", 
    attr(x, "mode"), "' mode"
  ))
  
  cat(slots0)
  cat(slots1)
  cat(slots2)
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


# Helpers -----------------------------------------------------------------

describe <- function(x, text = "", name = NULL) {
  if (is.null(name)) {
    name <- deparse(substitute(x))
    name <- strsplit(name, "\\$")[[1]]
    name <- tail(name, 1)
  }
  describe_(x, text = text, name = name)
}

describe_ <- function(x, text = "", name = NULL) {
  UseMethod("describe_")
}

describe_.array <- function(x, text = "", name = NULL) {
 paste0(" - ", name, ":\tAn ", paste0(dim(x), collapse = " x "), 
                 " array ", text, ".\n")
}

describe_.matrix <- function(x, text = "", name = NULL) {
  paste0(" - ", name, ":\tA ", paste0(dim(x), collapse = " x "),
         " matrix ", text, ".\n")
}

describe_.list <- function(x, text = "", name = NULL) {
  paste0(" - ", name, ":\tA list of length ", length(x), " ", text, ".\n")
}

describe_.data.frame <- function(x, text = "", name = NULL) {
  paste0(" - ", name, ":\tA ", paste0(dim(x), collapse = " x "),
         " data.frame ", text, ".\n")
}

describe_.data.table <- function(x, text = "", name = NULL) {
  paste0(" - ", name, ":\tA ", paste0(dim(x), collapse = " x "),
         " data.table ", text, ".\n")
}

describe_.discovery <- function(x, text = "", name = NULL) {
  myclass <- class(x)[[1]]
  expnames <- expnames(x)
  res <- resolution(x)
  res <- if (res %% 1e6 == 0) {
    paste0(res / 1e6, " Mb")
  } else if (res %% 1e3 == 0) {
    paste0(res / 1e3, " kb")
  } else {
    paste0(res, " bp")
  }
  pkg <- attr(x, "package")
  pkg <- if (is.null(pkg)) attr(x, "PACKAGE") else pkg
  pkg <- paste0(pkg, " ")
  string <- paste0("A ", pkg, "'", myclass,
                   "' object involving the following ",
                   length(expnames), " experiment(s):\n'",
                   paste0(expnames, collapse = "', '"), "' at a ",
                   "resolution of ", res, ".\n\n")
}
