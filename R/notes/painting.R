#painting.r

https://i.pinimg.com/564x/c0/1f/f1/c01ff14317d7a70d303a7bf37ef29ab4.jpg
rotate <- function(x) t(apply(x, 2, rev))

library(jpeg)
JPG = jpeg::readJPEG("data/nintchdbpict000277131838.jpg?w=811")
tmp = apply(JPG, 1:2, mean)
IMONABOAT = rotate(tmp)

save(antoni, file="../R/antoni.Rdat")

IMONABOATrect = IMONABOAT[,100:911]

load("../R/IMONABOAT.Rdat")
CR = colorRampPalette(c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000'))
image(  IMONABOATrect, col = rev(CR(10)))


IMONABOATrect = ResizeMat(IMONABOATrect,  c(350,350))
image(  tmp , col = rev(CR(10)))


antoni = IMONABOATrect
