stepped = do.call('rbind',lapply(1:length(tmpM), function(x){setNames(data.frame("window" = names(tmpM)[x],
                                                                       'border' = unname(tmpM[x])),
                                                            nm = c('window', 'border') ) }))
stepped$window = paste0(as.numeric(as.character(stepped$window))/1e6, 'Mb')
stepped$window = factor(stepped$window, levels = paste0(steps/1e6, 'Mb'))
arm =newTADs_chr19[grepl(rownames(newTADs_chr19), pattern = 'P'), 2]

stepped = rbind(stepped, data.frame('window' = 'arm', 'border' = arm))


intwo = as.numeric(names(table(unlist(tmpM[]))[table(unlist(tmpM[])) > 1]))
intwo = c(intwo, tmpM[[1]][tmpM[[1]] < steps[2]])
stepped = rbind(stepped, data.frame('window' = 'in2', 'border' = intwo))

stepped$inarm = stepped$border %in% arm


library(ggplot2)
ggplot(stepped, aes(x = border, y = window, col = inarm)) +
  geom_point(pch = "|")+
  scale_color_manual(values = c('red','black'))+
  RHWlib::RHWtheme('long')+
  coord_cartesian(xlim = c(0, E)) +
  geom_hline(yintercept = 10.5) +
  theme(aspect.ratio = 0.2)




