## ----global_options, include=FALSE-----------------------------------------
knitr::opts_chunk$set(fig.pos = '!h')

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

