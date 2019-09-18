select.subset.base <- function(exp, chrom, start, end){
  sel <- exp$ABS[exp$ABS[,1]==chrom &
                   exp$ABS[,2] >= start &
                   exp$ABS[,2] <= end,4]
  start.i <- min(sel); end.i <- max(sel)
  #create matrix

  tmp = exp$ICE[exp$ICE$V1 >= start.i &
                  exp$ICE$V1 <= end.i &
                  exp$ICE$V2 >= start.i &
                  exp$ICE$V2 <= end.i,]

  r1.s <- tmp[!is.na(tmp$V3),]

  r2.s <- r1.s
  sub.mat <- matrix(0, ncol=length(start.i:end.i), nrow=length(start.i:end.i))
  #fill the top triangle of the matrix
  index <- cbind(r1.s$V1-start.i+1,r1.s$V2-start.i+1)
  sub.mat[index] <- r1.s$V3
  #fill the bottom triangle of the matrix
  index <- cbind(r2.s$V2-start.i+1,r2.s$V1-start.i+1)
  sub.mat[index] <- r2.s$V3
  #take the middle of the window as the position int the genome
  pos <- (exp$ABS[exp$ABS[,4]%in%(start.i:end.i),2]+exp$ABS[exp$ABS[,4]%in%(start.i:end.i),3])/2
  list(x=pos, y=pos, z=sub.mat)
}

out = data.frame()

exp_nodt_10 = Hap1_WT_10kb
exp_nodt_10$ICE <- as.data.frame(exp_nodt_10$ICE)

exp_nodt_40 = Hap1_WT_40kb
exp_nodt_40$ICE <- as.data.frame(exp_nodt_40$ICE)
#exp = exp_nodt
for(d in seq(15e4,100e6, by = 5e6)){
  C <- paste0('chr', sample(1:13))[1]
  S <- 10
  E <- S + d
  # genova
  ptm     <- proc.time()
  tmp     <- select.subset(Hap1_WT_40kb, C, S, E)
  elapsed <- proc.time() - ptm
  outRow = data.frame('method'     = 'GENOVA',
                      'resolution' = Hap1_WT_40kb$RES,
                      'size'       = d,
                      'seconds'    = unname(elapsed[3]))
  out = rbind(out, outRow)
  # base R
  ptm     <- proc.time()
  tmp     <- select.subset.base(exp_nodt_40, C, S, E)
  elapsed <- proc.time() - ptm
  outRow = data.frame('method'     = 'base',
                      'resolution' = exp_nodt_40$RES,
                      'size'       = d,
                      'seconds'    = unname(elapsed[3]))
  out = rbind(out, outRow)


  #### 10kb

  # genova
  ptm     <- proc.time()
  tmp     <- select.subset(Hap1_WT_10kb, C, S, E)
  elapsed <- proc.time() - ptm
  outRow = data.frame('method'     = 'GENOVA',
                      'resolution' = Hap1_WT_10kb$RES,
                      'size'       = d,
                      'seconds'    = unname(elapsed[3]))
  out = rbind(out, outRow)
  # base R
  ptm     <- proc.time()
  tmp     <- select.subset.base(exp_nodt_10, C, S, E)
  elapsed <- proc.time() - ptm
  outRow = data.frame('method'     = 'base',
                      'resolution' = Hap1_WT_10kb$RES,
                      'size'       = d,
                      'seconds'    = unname(elapsed[3]))
  out = rbind(out, outRow)


  print(out)
}

library(reshape2)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
# m = melt(out)
# m = m[m$variable == 'seconds',c(1,3)]
# m =m %>% group_by(method) %>%
#   summarise(mean = mean(value),
#             sd = sd(value))
out.bk = out
out$resolution = as.factor(paste0(out.bk$resolution/1e3, 'kb'))
ggplot(out, aes( x = size, y = seconds, color = method)) +
  facet_grid(. ~ resolution) +
  #geom_boxplot() +
  geom_smooth(se = F) +
  scale_x_continuous(breaks = seq(1,100e6, length.out = 5),
                     labels = paste0(round(seq(1,100e6, length.out = 5)/1e6), '')) +
  geom_beeswarm() +
  scale_color_manual(values = c('skyblue','tomato')) +
  RHWlib::RHWtheme() +
  labs(y = 'seconds', x = 'distance (Mb)') +
  ggtitle('Region-lookup', subtitle = 'one region between 150kb and 100Mb')  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

