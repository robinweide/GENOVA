insulation.callTAD.plot <- function(exp , calls, chrom, start, stop, cut.off){
  hic.matrixplot(exp, chrom = chrom, start = start, end = stop, cut.off = cut.off, tads = calls$bedgraph[,1:3], tad.type = "lower")
  yMin <- calls$deltas %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(V4) %>% summarise(mi = min(V4)*1.5) %>% unlist()
  yMax <- calls$deltas %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(V4) %>% summarise(mi = max(V4)*1.5) %>% unlist()
  yMind <- calls$deltas %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(delta) %>% summarise(mi = min(delta)*1) %>% unlist()
  yMaxd <- calls$deltas %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(delta) %>% summarise(mi = max(delta)*1) %>% unlist()
  yMin <- min(yMin,yMind)
  yMax <- max(yMax,yMaxd)
  plot(calls0.05$deltas %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(V2, V4), type = 'l',axes = F, ylim = c(yMin,yMax))
  axis(2:4)
  lines(calls0.05$deltas %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(V2, delta), type = 'l', col = 'red')
  points(calls0.05$borders %>% filter(V1 == chrom, V2 >= start, V3 <= stop) %>% select(V2,V4), col = 'blue')
  abline(h = 0, col = 'red')
}
