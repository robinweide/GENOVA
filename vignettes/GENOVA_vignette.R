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

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(CTCF, 3), caption = 'A data.frame holding a standard BED6 format.'
)

## ----bigwrig, eval=F, echo = F---------------------------------------------
#  library(devtools)
#  install_github(repo ='bigwrig', username =  'jayhesselberth')

## ---- echo =F--------------------------------------------------------------
knitr::kable(
  head(martExport[,-c(1,2)], 5), caption = 'A data.frame holding the needed columns for plotting genes.'
)

## ---- echo=F---------------------------------------------------------------
options(scipen = 1e9)

## ---- echo=F---------------------------------------------------------------
options(scipen = 1)

## ----quantAPA, message=FALSE , fig.cap= "With quantifyAPA In the WAPL-knockout, we see an increase of contacts at the loop.",cache=T, fig.retina=T----
quantifyAPA_out <- quantifyAPA(APAlist = list('WT' = APA_Hap1_WT_extended,
                                              'WAPL' = APA_Hap1_WAPL_extended), 
                               pixWidth = 3)
print(quantifyAPA_out$stats)

# pot boxplot with base-R (ggplot2 would be also easy)
boxplot(split(quantifyAPA_out$data$value, f = quantifyAPA_out$data$sample), 
        col = c('darkgrey', 'red'), outline = F)

## ----sesh, echo = F--------------------------------------------------------
sessionInfo()

