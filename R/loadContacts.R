loadContacts = function(signalPath = NULL, indicesPath = NULL, juicerPath = NULL, resolution = 10e3, name){

  doJuicer = F
  doHiCpro = F

  # what are we dealing with?
  if(!is.null(juicerPath)){

    # a juicerpath is given. Does it exist?
    if(!file.exists(juicerPath)){
      stop("juicerPath doesn't point to an existing .hic-file.")
    }

    # is strawr installed?
    try_require('strawr', "loadContacts")

    doJuicer = T

    # did you also load in hicpro?
    if(!is.null(signalPath)){
      warning('Both hicpro- and juicebox-files provided. Will use juicebox.')
    }
  } else if(!is.null(signalPath)){

    # a signalPath is given. Does it exist?
    if(!file.exists(signalPath)){
      stop("signalPath doesn't point to an existing file.")
    }

    if(!file.exists(indicesPath)){
      stop("indicesPath doesn't point to an existing file.")
    }

    doHiCpro = T
  } else {
    stop('please provide either signalPath and indicesPath or a juicerPath.')
  }

  if(doJuicer){
    sig_ind = loadJuicer(juicerPath, resolution)
  }

  if(doHiCpro){
    sig_ind = loadHiCpro(signalPath, indicesPath)
  }

}



loadJuicer = function(juicerPath, resolution){

  require(strawr)

  # get metadata of juicer-file
  juicer_metadata = get_juicer_metadata(juicerPath)

  # check if resolution is found!
  if(!any(resolution == juicer_metadata[[2]])){
    stop('resolution is not found.\nAvailable resolutions are ',
         paste(rev(juicer_metadata[[2]]), collapse = ', '))
  }
  juicer_metadata[[1]] = juicer_metadata[[1]][juicer_metadata[[1]]$chrom != 'All',]

  expandedChromosomes = as.data.frame(t(combn(juicer_metadata[[1]]$chrom, m = 2)), stringsAsFactors = F)
  expandedChromosomes = rbind(as.data.frame(cbind(juicer_metadata[[1]]$chrom,juicer_metadata[[1]]$chrom)),
                              expandedChromosomes)

  juicerList = apply(expandedChromosomes, 1, function(ec){

    juicer_in <- tryCatch(
      {strawr::straw(norm = "KR",
                      fname =juicerPath,
                      chr1loc =  ec[2],
                      chr2loc =  ec[1],
                      unit = 'BP',
                      binsize = resolution)
      },error=function(cond) {

        return(NULL)
      })

    # remove NaN
    if(!is.null(juicer_in)){
      juicer_in = juicer_in[!is.nan(juicer_in[,3]),]
      if(nrow(juicer_in) == 0){
        juicer_in = NULL
      }
    }

    # decorate
    if(!is.null(juicer_in)){

      juicer_in$chrom_x = ec[1]
      juicer_in$chrom_y = ec[2]

      juicer_in = juicer_in[, c(4,1,5,2,3)]
    }
    juicer_in

  })

  juicer_data = data.table::rbindlist(juicerList)

  # split into index and signal-files




}




# taken from ggplot
try_require <- function(package, fun) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }

  stop("Package `", package, "` required for `", fun , "`.\n",
       "Please install and try again.", call. = FALSE)
}

get_juicer_metadata = function(juicerPath){
  file2read = file(juicerPath, "rb")
  number_of_integers_in_file = 125

  # check if magicstring is there
  magic_string = readBin(file2read, character(), n = 1, size = 3)
  if(magic_string != 'HIC'){
    stop('magic sting comparison failed. Check the .hic file!')
  }

  hic_version = readBin(file2read, integer(), n = 1, size = 4)

  master_index = readBin(file2read, integer(), n = 1, size = 8)

  genome = readBin(file2read, character())

  n_attributes = readBin(file2read, integer(), n = 1, size = 4)

  n_chrs = readBin(file2read, integer(), n = 1, size = 4)


  chrom_name_length = lapply(1:n_chrs, function(x){
    chrom = readBin(file2read, character(), size = 1)
    len = readBin(file2read, integer(), n = 1)
    data.table::data.table(chrom, len)
  })
  chrom_name_length = data.table::rbindlist(chrom_name_length)

  n_BP_Res  = readBin(file2read, integer(), n = 1, size = 4)


  resolutions = sapply(1:n_BP_Res, function(x){
    readBin(file2read, integer(), n = 1)
  })

  close(file2read)

  return(list(chrom_name_length, resolutions))
}
