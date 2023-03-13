

### Plot
map_phys_position <- function(chr = 18){

 crossNBH <- subset(cross_NBH, chr = chr) 
 sm_NBH <- scanone(crossNBH, pheno.col=4, model="binary",method="mr")
 smy_NBH <- as.numeric(gsub(".*:","",markernames(crossNBH)))

 crossELR <- subset(cross_ELR, chr = chr) 
 sm_ELR <- scanone(crossELR, pheno.col=4, model="binary",method="mr")
 smy_ELR <- as.numeric(gsub(".*:","",markernames(crossELR)))

 crossNEW <- subset(cross_NEW, chr = chr) 
 sm_NEW <- scanone(crossNEW, pheno.col=4, model="binary",method="mr")
 smy_NEW <- as.numeric(gsub(".*:","",markernames(crossNEW)))

 crossBRP <- subset(cross_BRP, chr = chr) 
 sm_BRP <- scanone(crossBRP, pheno.col=4, model="binary",method="mr")
 smy_BRP <- as.numeric(gsub(".*:","",markernames(crossBRP)))
  
 xlims <- c(0,max(smy_NBH,smy_ELR,smy_NEW,smy_BRP, na.rm = T))

#pdf(paste0('mr',chr,'.pdf'), height = 7, width = 21)
#plot_test_svg(paste0(chr,'mr'), width = 14, height = 7)
#plot_test(paste0(chr,'mr'), width = 1000)
 plot(xlims, c(0,18),type="n", bty="n", axes=F, ylab = NA, xlab = NA)
  axis(1,pos=0, lwd = 4, labels = F)
  axis(2,pos=0, lwd = 2, at = c(0,5,10), labels = F)
  points(smy_NBH, sm_NBH$lod, col = pop_colors['NBH_tol'], pch = 16)
  points(smy_BRP, sm_BRP$lod, col = pop_colors['BRP_tol'], pch = 16)
  points(smy_NEW, sm_NEW$lod, col = pop_colors['NEW_tol'], pch = 16)
  points(smy_ELR, sm_ELR$lod, col = pop_colors['ELR_tol'], pch = 16)
 dev.off()

}





## ## cover part of scale single fig
## trans <- function(x){pmin(x,6) + 0.5*pmax(x-6,0)}
## yticks <- c(0, 3, 6, 12.5, 15, 17.5, 20)
## longpops$transy <- trans(longpops$lod)
## ggplot(longpops, aes(x=pos, y=transy, group=popl)) +
##   geom_point(aes(color=popl)) + 
##   #geom_smooth(aes(color=popl), method = 'loess', span = 0.3)
##   stat_smooth(aes(color=popl),method = 'lm', formula = y ~ ns(x, 10)) +
##   geom_rect(aes(xmin=0, xmax=max(pos), ymin=5.9, ymax=8), fill="white") +
##   scale_y_continuous(limits=c(0,NA), breaks = trans(yticks), labels = yticks) +
##   theme_classic()

