#' Construct a HiC-experiment.
#'
#' Make an object which holds the most needed information of a HiC-experiment at a given resolution.
#'
#' @param signalPath Full path to a HiC-pro-like matrix-file or .sig from juicerToGenova.py
#' @param indicesPath Full path the HiC-pro-like index-file or .bed from juicerToGenova.py
#' @param name The name of the sample.
#' @param Znorm Normalise the matrices per-chromosome with a Z-score.
#' @param ignore.checks Skip the checks for empty matrices: EXPERT-ONLY!
#' @param color Color associated with sample.
#' @param comments A place to store some comments.
#' @param centromeres A data.frame with three columns per chromosome: chromosome name, start-position and end-position of the centromeric region.
#' @param BPscaling Scale contacts to have agenome-wide sum of [BPscaling] reads (default: 1000M). Set to NULL to skip this.
#' @note
#' Some reference genomes have very small "random" or "patch" chromosomes, which can have zero contacts mapped to it (at certain resolutions). Construct.experiment checks this and omits these chromosomes in the resulting experiment-object. The RMCHROM-flag will also be set to TRUE: this will help other GENOVA-functions to deal better with this problem. There is a slight performance-cost during the construction of the experiment, however. Therefore, experienced users can set ignore.checks to TRUE to skip all of this, keeping in mind that some functions will not work properly/at all is there are these zero-coverage chromosomes in their data.
#' @examples
#' WT_10kb <- construct.experiment(ignore.checks = T, signalPath = "WT_10kb_iced.matrix", indicesPath = "WT_10kb_abs.bed", name = "WT", color = "black")
#' @return An experiment-object, which is a named list of contacts, indices and metadata for a Hi-C matrix of a given sample at a given resolution.
#' @export
construct.experiment <- function(signalPath, indicesPath, name, Znorm = F, ignore.checks = F, centromeres = NULL, color = 1, comments = NULL, BPscaling = 1e9) {
  # Check if files exist
  if (!file.exists(signalPath)) {
    stop("Signal file not found.")
  }
  if (!file.exists(indicesPath)) {
    stop("Index file not found.")
  }

  ICE <- read.hicpro.matrix(signalPath, norm = BPscaling)
  ABS <- data.table::fread(indicesPath, header = F, data.table = F)

  # check for similar binsizes BED and ICE
  ICEmaxID <- max(max(ICE$V1), max(ICE$V2))
  ABSmaxID <- max(ABS[, 4])
  if (ICEmaxID > ABSmaxID) {
    stop("Undefined HiC-bins.\nWe checked if the highest Hi-C bin in the signal-file could be found in the indices-file: it couldn't.\nPlease check if indices and signal are from the same reference genome and resolution!")
  }
  if (abs(ICEmaxID - ABSmaxID) > ICEmaxID * 0.1) {
    warning("Undefined BED-bins.\nWe checked if the highest Hi-C bin in the indices-file could be found in the signal-file: it couldn't.\nPlease check if indices and signal are from the same reference genome and resolution!")
  }
  chromVector <- as.character(unique(ABS$V1))
  res <- as.numeric(median(ABS$V3 - ABS$V2))

  if (!is.null(centromeres)) {
    centromeres <- clean_centromeres(centromeres, res)
  }

  # check is all chromosomes have actual data
  RMCHROM <- F
  if (!ignore.checks) {
    IIDX <- unique(c(ICE$V1, ICE$V2))
    for (C in chromVector) {
      CIDX <- ABS[ABS[, 1] == C, 4]

      if (!any(IIDX %in% CIDX)) {
        warning(paste0("No contacts of ", C))
        ABS <- ABS[ !ABS[, 1] == C, ]
      }
    }
  }
  chromVector.bk <- chromVector
  chromVector <- as.character(unique(ABS$V1))

  if (!length(chromVector) == length(chromVector.bk)) {
    RMCHROM <- T
    warning("Some chromosomes have been removed, since these had no data in the signal-matrix.\nTo turn this off, set ignore.checks to TRUE.")
  }

  # Z-score normalisation
  if (Znorm) {
    tmpABS <- data.table::as.data.table(ABS)
    colnames(tmpABS) <- c("C1", "S", "E", "V1")
    data.table::setkey(tmpABS, "V1")
    tmp <- merge(ICE, tmpABS[, c(1, 4)])
    data.table::setkey(tmp, "V2")

    colnames(tmpABS) <- c("C2", "S", "E", "V2")
    data.table::setkey(tmpABS, "V2")
    tmp <- merge(tmp, tmpABS[, c(1, 4)])

    ICEwithChrom <- tmp[tmp$C1 == tmp$C2, c(2, 1, 3, 4)]
    colnames(ICEwithChrom) <- c("V1", "V2", "V3", "C")
    data.table::setkey(ICEwithChrom, "C")

    ICEwithChrom$D <- abs(ICEwithChrom$V1 - ICEwithChrom$V2)
    thisZ <- dplyr::mutate(dplyr::group_by(ICEwithChrom, C, D), Z = scale(V3))

    thisZ <- data.table::as.data.table(thisZ[, c(1, 2, 6)])
    colnames(thisZ) <- c("V1", "V2", "V3")
    data.table::setkey(thisZ, "V1", "V2")
    ICE <- thisZ
  }

  # Contruct list
  structure(list(
    # Iced HiC-matrix in three-column format (i.e. from HiC-pro)
    ICE = ICE,

    # HiC-index in four-column format (i.e. from HiC-pro)
    ABS = ABS,

    # Name of sample
    NAME = name,

    # Resolution of sample
    RES = res,

    # Available chromosomes
    CHRS = chromVector,

    # color of sample (optional, but recommended for running RCP)
    COL = color,

    # comments
    COMM = comments,

    # Vector of masked bins
    MASK = vector(),

    # DF with chr, start, stop
    CENTROMERES = centromeres,

    # Are there now missing bins?
    RMCHROM = RMCHROM,

    # Did we normalise with Z?
    ZSCORE = Znorm
  ), class = "contacts", package = "GENOVA")
}

clean_centromeres <- function(centros, resolution) {
  # Essentially does the same as `reduce(a_granges_object, min.gapwidth = resolution)

  centros <- as.data.frame(centros)

  # Drop unused factor levels
  if (is.factor(centros[, 1])) {
    centros[, 1] <- droplevels(centros[, 1])
  }

  # Set column names and order on chromosome and start
  centros <- setNames(centros, paste0("V", 1:3))
  centros <- centros[order(centros[, 1], centros[, 2]), ]

  # Test wether adjacent entries should be merged
  centros$merge <- c(FALSE, vapply(seq_len(nrow(centros) - 1), function(i) {

    # Is the difference between end and next start sub-resolution?
    test <- abs(centros[i, 3] - centros[i + 1, 2]) < resolution

    # Are the entries on the same chromosome?
    test && centros[i, 1] == centros[i + 1, 1]
  }, logical(1)))

  # Merge the TRUE entries
  while (any(centros$merge)) {
    # What is the next entry to be merged?
    i <- which(centros$merge)[1]

    # Replace end of previous entry with end of current entry
    centros[i - 1, 3] <- centros[i, 3]

    # Remove current entry
    centros <- centros[-i, ]
  }

  # Return centromeres
  centros[, 1:3]
}
