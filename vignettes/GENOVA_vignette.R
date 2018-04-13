## ----global_options, include=FALSE-----------------------------------------
knitr::opts_chunk$set(fig.pos = 'h')

## ----echo = F--------------------------------------------------------------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(4, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(1.5,y,10,y+1/scale, col=lut[i], border=NA)
    }
}

## ----echo=FALSE, out.width='100%', fig.align='center'----------------------
knitr::include_graphics('/DATA/users/r.vd.weide/github/GENOVA/t1logo')

## ---- echo=T, warning=FALSE, error=F, results='hide'-----------------------
# devtools::install_github("robinweide/GENOVA", ref = 'dev')
library(GENOVA)

## ----centromere0, cache=T, echo = F----------------------------------------
centromeres = read.delim('data/hg19_cytobandAcen.bed', 
                         sep = '\t', 
                         h = F, 
                         stringsAsFactors = F)

## ----peakEXP3, collapse=F, results='markup', echo = F----------------------
str(Hap1_WT_40kb, width = 60,   vec.len=1, strict.width = 'wrap')

## ----J2G, eval=FALSE, highlight=FALSE--------------------------------------
#  # Convert data from Sanborn et al. normalised at 10kb resoltion:
#  juicerToGenova.py -C ucsc.hg19_onlyRealChromosomes.noChr.chromSizes \
#  -JT ~/bin/juicer/AWS/scripts/juicebox_tools.7.0.jar \
#  -H ~/Downloads/Sanborn_Hap1_combined_30.hic \
#  -R 10000 \
#  -force TRUE \
#  -norm KR \
#  -O Sanborn_Hap1_combined_30.hic_10kb_KR

## ---- echo=F---------------------------------------------------------------
options(scipen = 1)

## ----saddle, message=FALSE,  cache=T---------------------------------------
H3K27acPeaks = read.delim('data/H3K27ac_WT.narrowPeak', h = F)
CS = read.delim('/DATA/oidBackup/WAPL_Project/Hi-C/analysis/PCA/PCA/100kb/WT.PC1.bedGraph', h = F, skip = 1)
saddle_WT = saddleBins(exp = Hap1_WT_40kb, 
                       ChIP = H3K27acPeaks,CS = CS,
                       chromsToUse = paste0('chr', 1:15),
                       nBins = 50, 
                       verbose = T)

saddle_WAPL = saddleBins(exp = Hap1_WAPL_40kb, 
                       ChIP = H3K27acPeaks, CS = CS,
                       chromsToUse = paste0('chr', 1:15),
                       nBins = 50,  
                       verbose = F)

## ----saddlePlot, message=FALSE,  fig.asp=.65, cache=T, warning=FALSE, fig.cap= "A saddle-plot", fig.width=8----
visualise.saddle(SBoutList = list(saddle_WT, saddle_WAPL),
                 crossLines = T, 
                 addText = T,
                 zlim = c(-0.5,0.5), 
                 EVlim = c(-50,50))

## ----saddleStrength, message=FALSE,  cache=T, warning=FALSE, fig.cap= "The per-arm compartment strength", fig.small = T----
visualise.compartmentStrength(list(saddle_WT,
                                   saddle_WAPL))

visualise.compartmentStrength(list(saddle_WT,
                                   saddle_WAPL), showInteractions = T)

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(CTCF, 3), caption = 'A data.frame holding a standard BED6 format.'
)

## ----bigwrig, eval=F, echo = F---------------------------------------------
#  library(devtools)
#  install_github(repo ='bigwrig', username =  'jayhesselberth')

## ----HMPchip, message=T , fig.cap= "Hi-C matrixplot: ChIPseq. A BED-file of CTCF-sites is plotted at the top and a coverage-track of SMC1 ChIP-seq is plotted beneath this. The symmAnn-option leads to the same tracks being plotted on the left.",cache=T,fig.asp=1, dev = 'png', dpi=300,fig.small = F----
hic.matrixplot(exp1 = Hap1_WT_10kb,
               chrom = 'chr7',  start = 26.75e6,  end=28.5e6, 
               loops = WT_Loops, # see APA
               loops.color = '#998ec3', # purple loops
               loops.type = 'upper', # only plot in upper triangle
               loops.resize = 20e3, # expand for visibility
               type = 'triangle',
               chip = list('data/SMC1_WT.bw', # inner top
                           CTCF),# outer-top
               symmAnn = F, # place annotations also on left side
               cut.off = 65) # upper limit of contacts

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(martExport[,-c(1,2)], 5), caption = 'A data.frame holding the needed columns for plotting genes.'
)

## ----iiDIFF2, message=FALSE , fig.cap= "Differential TAD-analysis: scatterplot. Experiment 2 (WAPL) has more interactions between neighbouring TADs compared to wild type.",dev = 'png', dpi=300, cache=T, fig.retina=T----
par(mfrow = c(1,2), pty = 's')
differential.TAD.scatterplot(exp1 = TAD_N_WT, # x
                            exp2 = TAD_N_WAPL, 
                            allData = T, 
                            main = 'allData == T') # y
differential.TAD.scatterplot(exp1 = TAD_N_WT, # x
                            exp2 = TAD_N_WAPL, 
                            allData = F, 
                            main = 'allData == F') # y

## ---- echo=F---------------------------------------------------------------
options(scipen = 1e9)

## ---- echo=F---------------------------------------------------------------
options(scipen = 1)

## ----sesh, echo = F--------------------------------------------------------
sessionInfo()

