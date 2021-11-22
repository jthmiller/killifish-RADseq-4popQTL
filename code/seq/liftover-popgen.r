#!/bin/R

## Rscript to write popgen files ##
load("/home/nreid/popgen/variants/bowfree/angsd/5kb1kboutliers.Robj")
stats <- c('pbstat', 'piper','tajstat','pfst','picou')
dir <- '/home/jmiller1/QTL_Map_Raw/popgen/outlier_table/'

for (i in stats){
        write.table(get(i),file=paste(dir,i,'.txt',sep=''),sep='\t',row.names =F, quote=F)

} 
