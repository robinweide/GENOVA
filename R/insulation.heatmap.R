# taken from https://github.com/robinweide/RHWlib

#' Plot insulation-heatmaps.
#'
#' Takes a list of insulation-scores and plots a sorted heatmap and/or
#' average profile.
#'
#' @author Robin H. van der Weide, \email{r.vd.weide@nki.nl}
#' @param insulationList A named list of results from genome.wide.insulation().
#' @param bed A BED-df with TAD-calls. Only the first three columns
#' (chrom, start & stop) are important. Will use unique borders.
#' @param borders A df with borders-calls. Only the first three columns
#' (chrom, start & stop) are important: the middle of the region will be used
#' for the alignment. If given, this function will ignore the bed-agument. Will
#' use all borders, regardless of uniqueness.
#' @param focus Sort on which sample? (give the index-number of the list)
#' @param sortWidth Percentage of columns to sort on (aligned on center)
#' @param whatToPlot Do you want to plot a profile, heatmap or both?
#' @param profileFunct Function to make profile-plots (mean, median, sum)
#' @param flank Amount of flanking bp.
#' @param zlim The heatmap-zlims c(min,max).
#' @param profileCols Vector of line-colors for profiles
#' @param profileZlim The profile-zlims c(min,max).
#' @param heatmapCols A vector of colors for the heatmap-gradient.
#' @param leftNorm Normalise data on median on most upstream bin.
#' @param verbose Should this function be chatty?
#' @return A plot plus an (invisible) dataframe of the underlying matrix.
#' @examples
#' # Get the insulation score with window-size 25 of two experiments.
#' WT_10kb_ins = genome.wide.insulation(hic = WT_10kb, window.size = 25)
#' SCC4_10kb_ins = genome.wide.insulation(hic = SCC4_10kb, window.size = 25)
#'
#' # Store insulation-scores in a named list
#' inList = list(WT = WT_10kb_ins, SCC4 = SCC4_10kb_ins )
#'
#' # Plot the heatmap and profile and store the underlying df in df.out.
#' df.out = insulation.heatmap(inList, bed = WT_TADs,
#'                             zlim = c(-.5,0.25), profileZlim = c(-.75,-.1))
#' @export
insulation.heatmap <- function(insulationList, bed = NULL, borders = NULL, focus = 1, sortWidth = 10,
                               whatToPlot = 'both',  profileFunct = mean,
                               flank=500e3, zlim = c(-1,1), profileCols = NULL,
                               profileZlim = NULL, heatmapCols = rev(c('white',
                                                                       'white',
                                                                       'white',
                                                                       'red',
                                                                       'black')),
                               leftNorm =F, verbose = F ){

  require(ggplot2)
  require(reshape2)

  # What are we plotting? ------------------------------------------------------
  if(!whatToPlot %in% c('profile', 'heatmap', 'both')){
    stop('whatToPlot should be "profile", "heatmap" or "both".')
    }

  # Get borders ----------------------------------------------------------------
  BORDERS = NULL

  if(!is.null(borders)){
    # use borders
    BORDERS = data.frame(chrom = borders[,1], pos = apply(borders[,2:3],1, mean))
  } else if(!is.null(bed)){
    # use bed
    BORDERS = unique(rbind(setNames(bed[,c(1,2)], c('chrom', 'pos')),
                           setNames(bed[,c(1,3)], c('chrom', 'pos'))))
  } else {
    stop('Supply either bed or borders!')
  }

  # Get matrixes ---------------------------------------------------------------
  sampleList = lapply(insulationList, function(insdDat){

    align.insulation( ins.data = insdDat,
                      locs = BORDERS,
                      flank.length=flank,
                      verbose = verbose)

  })

  # Get order from focus -------------------------------------------------------

  if(focus > length(insulationList)){
    stop("Focus is outside of the number of samples!")
  } else {
    focusYoungGrasshopper <- sampleList[[focus]]

    focusYoungGrasshopper <- as.matrix(focusYoungGrasshopper[,-c(1,2)])
    # which of columns to sort on?
    perSample = ncol(focusYoungGrasshopper)
    nSortCol <-  ((perSample)/100)*sortWidth
    halfnSortCol <- round(nSortCol/2)
    midpoint <- ((ncol(focusYoungGrasshopper)-1)/2)+1
    sortCols <- (midpoint-halfnSortCol):(midpoint+halfnSortCol)

    rowSums <- apply(as.matrix(focusYoungGrasshopper)[, sortCols],1, mean)
    foundOrder <- base::order(rowSums,decreasing = T) # in/de flipped for ggplot
  }

  # order individual samples
  for(i in 1:length(sampleList)){
    sampleList[[i]] <- sampleList[[i]][foundOrder,]
    sampleList[[i]]$rowIDX = 1:nrow(sampleList[[i]])

  }

  # Melting the data -----------------------------------------------------------

  df = dplyr::bind_rows(sampleList, .id = 'sample'); outDF = df

  df = reshape2::melt(df[,-c(2:3)], id.vars = c('sample','rowIDX'))
  df$variable = as.numeric(gsub(df$variable, pattern = 'V', replacement = ''))

  df$sample <- factor(df$sample, levels = names(insulationList))

  # Normalise to first bin -----------------------------------------------------
  if(leftNorm){
    df = normHM2left(HM = df)
  }

  # Make the profile -----------------------------------------------------------
  groupedDF <- dplyr::group_by(df, sample, variable)
  groupedDF[is.na(groupedDF$value), "value"] = 0
  minVal = min(groupedDF[groupedDF$value != -Inf, "value"])
  groupedDF[groupedDF$value == -Inf, "value"] = minVal
  profDF <- dplyr::summarize(groupedDF, val = profileFunct(value))

  # GGPlot features- -----------------------------------------------------------
  if(is.null(profileCols)){
    profileCols <- rep('black', length(insulationList))
  }

  flankingSize = (perSample -1 )/2
  upTickPos = flankingSize/2
  downTickPos = perSample - upTickPos
  upTickLab = -1*((flank/2)/1e3)
  downTickLab = (flank/2)/1e3

  tickLabs = c(paste0(upTickLab, 'kb'),  "", paste0(downTickLab, 'kb'))
  tickPoss = c(upTickPos+1,  perSample - flankingSize ,  downTickPos )

  # Plot profile ---------------------------------------------------------------
  PRFLS <- NULL

  if(whatToPlot == 'profile' | whatToPlot == 'both'){

    if(is.null(profileZlim)){ profileZlim = zlim  }

    PRFLS <- ggplot2::ggplot(profDF, ggplot2::aes(x = variable, y = val, col = sample)) +
              ggplot2::geom_line()  +
              ggplot2::scale_x_continuous(expand=c(0,0),
                                          breaks = tickPoss,
                                          labels = tickLabs) +
              ggplot2::scale_y_continuous(expand=c(0,0), limits = profileZlim ) +
              ggplot2::labs( x = '', y = '') + ggplot2::guides(col = F) +
              ggplot2::scale_colour_manual(values = profileCols)+
              ggplot2::theme(aspect.ratio = 1) +
              ggplot2::facet_grid(. ~ sample ) +
              GENOVA:::GENOVA_THEME() +
              ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))+
              ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                             axis.title.x=ggplot2::element_blank())
  }

  # Plot HM --------------------------------------------------------------------
  HTMPS <- NULL

  if(whatToPlot == 'heatmap' | whatToPlot == 'both'){

    df$value[df$value < zlim[1]] = zlim[1]
    df$value[df$value > zlim[2]] = zlim[2]

    HTMPS <- ggplot2::ggplot(df, ggplot2::aes(variable, y = rowIDX, fill = value)) +
      ggplot2::geom_raster(interpolate = T) +
      ggplot2::facet_grid(. ~ sample) +
      ggplot2::scale_x_continuous(expand=c(0,0), breaks =tickPoss,
                                  labels = tickLabs) +
      ggplot2::scale_y_continuous(expand=c(0,0), breaks = NULL) +
      ggplot2::labs( x = '', y = '') +
      GENOVA:::GENOVA_THEME() +
      ggplot2::guides(fill = F) + ggplot2::theme_linedraw()+
      ggplot2::scale_fill_gradientn(colours = heatmapCols)+
      ggplot2::theme(strip.text = ggplot2::element_blank() ,
                     strip.background = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.margin = ggplot2::unit( c(0,0,0,0) , units = "lines" ))
  }


  if(whatToPlot == 'both'){
    df$dummy = factor('heatmap', levels = c('profile','heatmap'))
    profDF$dummy = factor('profile', levels = c('profile','heatmap'))


    BOT = ggplot2::ggplot(df, ggplot2::aes(variable, y = rowIDX, fill = value)) +
      ggplot2::geom_vline(xintercept = flankingSize+1, lty = 3)+
      ggplot2::geom_raster(interpolate = T) +
      ggplot2::facet_grid(dummy ~ sample, scales = 'free_y') +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::scale_x_continuous(expand=c(0,0),
                                  breaks = tickPoss,
                                  labels = tickLabs) +
      geom_line(data = profDF, mapping = aes(x = variable, y = val, col = sample), inherit.aes = F) +
      geom_blank(data    = data.frame(dummy = 'profile', x = 1, y = profileZlim),
                 mapping = aes(x = x, y = y),
                 inherit.aes = F)+
      geom_blank(data    = data.frame(dummy = 'heatmap', x = 1, y = zlim),
                 mapping = aes(x = x, y = y),
                 inherit.aes = F) +
      coord_cartesian(expand = F) +
      GENOVA:::GENOVA_THEME() +
      ggplot2::labs( x = '', y = '') +
      ggplot2::guides(col = F, fill = F) +
      ggplot2::scale_colour_manual(values = profileCols)+
      ggplot2::scale_fill_gradientn(colours = heatmapCols)+
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     axis.title.x     = ggplot2::element_blank(),
                     strip.text.y = element_blank())
    plot(BOT)
    invisible(outDF)
  } else if(whatToPlot == 'heatmap'){
    plot(HTMPS)
    invisible(outDF)
  } else {
    plot(PRFLS)
    invisible(outDF)
  }


}

normHM2left = function(HM, profileFunct = mean){

  tmp = lapply(split(HM, f = HM$sample), function(x){

    x$value[is.na(x$value)] = 0
    minvval = min(x$value[is.finite(x$value)])
    x$value[x$value == -Inf] = minvval

    left = profileFunct(unlist(x[x$variable ==1, 'value']))
    x$norm = x$value - (left-0)
    x$value = x$norm
    x$norm = NULL

    x

  })
  return(dplyr::bind_rows(tmp))
}







align.insulation.chrom <- function( ins.data, locs, flank = 10 ){
  n <- findInterval(locs[,2], ins.data[,2])

  #create a vector for the flanking sequence
  add <- -flank:flank
  #repeat it for every element in n
  add.long <- rep(add, length(n))
  #repeat every element in n for the length of add
  n.long <- rep(n, each=length(add))

  #add the selected value and the flanking sequences
  i.sel <- n.long + add.long

  # add a dummy out-of-range insu-entry
  ins.data = rbind(ins.data, data.frame('chrom' =  ins.data[1,1],
                             'start' = 007,
                             'end'   =  1337,
                             'insulation' = NA))
  i.sel[i.sel < 1] = nrow(ins.data)

  #put the selected elements in matrix
  align.mat <- matrix(ins.data[i.sel,4], nrow=length(n), ncol=length(add),
                      byrow=T)

  # add zeroes for NA's
  align.mat[is.na(align.mat)] = 0

  return(align.mat)
}

align.insulation <- function(ins.data, locs, flank.length = 200e3, verbose = F){
  #elaborate way of selecting the resolution of the Hi-C matrix
  res <- as.numeric(names(tail(sort(table(ins.data[,3]-ins.data[,2])),1)))
  flank = flank.length/res

  chrom.vec <- unique(ins.data[,1])
  align.mat <- NULL

  #for( chr in chrom.vec){

  AM = lapply( as.character(chrom.vec), function(chr){

    chromLoc = locs[locs[,1]==chr,]

    sub.mat <- align.insulation.chrom( ins.data = ins.data[ins.data[,1]==chr,],locs = chromLoc , flank = flank )

    DF = as.data.frame(sub.mat)
    DF = cbind(chromLoc, DF)

    DF
  }  )

  return(do.call('rbind', AM))
}



