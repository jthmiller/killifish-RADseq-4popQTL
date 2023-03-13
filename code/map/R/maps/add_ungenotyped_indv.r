
#!/bin/R
### Map QTLs 1 of 3
pop <- 'NBH'
debug.cross <- T
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
setwd('/home/jmiller1/QTL_agri')
#source("/home/jmiller1/QTL_agri/MAP/control_file.R")
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")

library('qtl')
mpath <- '/home/jmiller1/QTL_agri/data'


### Add in the genotypes performed by whoi
add.genotypes(cross, pheno.csv, phenos, genos){
    cross.genos <- cross$geno


}


#Pheno,sex,ID,bin,pheno_norm,1:210602,1:317181,1:363497,1:515871,1:857165,1:1066755,1:1317831,1:1582448,1:1912131,1:1926606,1:1959014,1:2040515,1:2040694,1:2040701,1:2040717,1:2040771,1:2040898,1:2040912,1:2049642,1
#,,,,,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
#,,,,,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
#4,0,ELR_10871,NA,NA,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,AB,AB,AA,AB,AB
#0,0,ELR_10968,NA,NA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA
#0,0,ELR_10978,NA,NA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,AB,AB,AB,AB,AB
#0,0,ELR_10991,NA,NA,AB,AB,AB,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,BB,BB,BB,AB,AB,AB,AA,AB,AB,AB,AB,AB,AA,AA,AB,BB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,BB,BB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB
#4,0,ELR_10884,NA,NA,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,AA,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AA,BB,AB,AB,AB,AB,AB,AB,AB,AB,AB,BB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB
#5,0,ELR_10928,NA,NA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA,AA
#1,0,ELR_11578,NA,NA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB
#0,0,ELR_11114,NA,NA,AB,AB,AB,AB,BB,AB,AB,AB,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,BB,AB,AB,AB,AA,AA,AA,AB,AB,AB,AA,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB,AB
#0,0,ELR_10977,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#0,0,ELR_10988,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#0,0,ELR_10983,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#5,0,ELR_10980,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#5,0,ELR_10974,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#1,0,ELR_10987,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#4,0,ELR_10971,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,
#1,0,ELR_10972,NA,NA,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,


add.ungenotyped.ind <- function(cross, pheno.csv, phenos, filename, pop){

    cross.phenos <- cross$pheno
    rownames(cross.phenos) <- cross.phenos$ID 

    pheno.csv <- pheno.csv[grep(pop,rownames(pheno.csv)),phenos]
    
    ## names of phenotypes to add
    col.add.phenos <- colnames(pheno.csv)[which(!colnames(pheno.csv)%in%colnames(cross.phenos))]
    
    ## values of new phenotypes
    val.add.phenos <- setNames(pheno.csv[,col.add.phenos], rownames(pheno.csv))

    ## which ids are not in current cross object
    ids.to.add <- as.character(pheno.csv$ID[which(!pheno.csv$ID %in% cross$pheno$ID)])

    ## add new phenotype to previous cross
    cross$pheno[,col.add.phenos] <- pheno.csv[cross$pheno$ID,col.add.phenos]

    ## added individual phenos
    ## common.ids <- colnames(pheno.csv)[colnames(pheno.csv) %in% colnames(cross$pheno)]

    uniq.pheno <- colnames(cross$pheno)[ which(!colnames(cross$pheno) %in% colnames(pheno.csv))]
    pheno.csv[,uniq.pheno] <- NA
    new.pheno <- pheno.csv[ids.to.add,colnames(cross$pheno)]
    new.pheno$ID <- as.character(new.pheno$ID)
    new.pheno$GT_NG_ALT <- as.character(new.pheno$GT_NG_ALT)

    empty.genos <- rep('-',sum(nmar(cross)))
    names(empty.genos) <- markernames(cross)

    ## write cross object with new pheno, and then append new individuals
    write.cross(cross, format="csv", filestem = filename )


    # bind together phenos and empty genos, append to file
    lapply(ids.to.add,function(id){ 
    
        text <- paste(unlist(c(new.pheno[id,],empty.genos)), collapse = ',')

        cat(text , file = paste0(filename,'.csv'), sep = "\n", fill = F, labels = NULL, append = T)

        }
    )

    cross <- read.cross(
        file = paste0(filename,'.csv'),
        format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
        estimate.map = FALSE
    )

    return(cross)

}



pheno.csv <- read.csv('/home/jmiller1/QTL_agri/data/allpop_pheno_whoiGT.csv')
rownames(pheno.csv) <- pheno.csv[,1]
phenos <- colnames(pheno.csv)[1:3]

## cross
pop <- 'NBH'
mapfile <- paste0(pop,'_reorder_noimp_nopar')
filename <- file.path(mpath,mapfile)
cross <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

## ug cross
pop <- 'NBH' 
mapfile <- paste0(pop,'_ungenotyped')
filename <- file.path(mpath,mapfile)
cross_ug <- add.ungenotyped.ind(cross = cross, pheno.csv = pheno.csv, phenos = phenos, filename = filename, pop = 'NBH')
