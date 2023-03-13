
#!/bin/R

### Compare model fit that includes ungenotyped individuals phenotype 
### 

## Only options for ungenotyped is imp and hk. hk should be avoided, imp needs stratfied permutations

## Binary trait mapping

### Map QTLs 1 of 3
pop <- 'NBH'
debug.cross <- T
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
setwd('/home/jmiller1/QTL_agri')
source("/home/jmiller1/QTL_agri/MAP/R/control_file.R")
cores <- 16

library('qtl')
mpath <- '/home/jmiller1/QTL_agri/data'

## noperms cross ###############################################################
mapfile <- paste0(pop,'_reorder_noimp_nopar')
filename <- file.path(mpath,mapfile)
noperms <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

## noperms cross ###############################################################
mapfile <- paste0(pop,'_ungenotyped')
filename <- file.path(mpath,mapfile)
cross_ug <- read.cross(file = paste0(mapfile,'.csv'), format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
################################################################################

#### PHENO #####################################################################
noperms$pheno$bin <- ifelse(noperms$pheno$Pheno > 2, 1 , 0)
noperms$pheno$pheno_norm <- round(nqrank(noperms$pheno$Pheno))
################################################################################

#### PHENO #####################################################################
cross_ug$pheno$bin <- ifelse(cross_ug$pheno$Pheno > 2, 1 , 0)
cross_ug$pheno$pheno_norm <- round(nqrank(cross_ug$pheno$Pheno))
################################################################################

i <- 1 ; plotit(cross_ug,'_cross_ug')
i <- 2 ; plotit(noperms,'_noperms')

#spuriously mapped markers
erprob_noperm <- 0.05 #determined from mapping step
noperms <- calc.genoprob(noperms, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
sone_bin_noperm <- scanone(noperms, pheno.col=4, method="em", model="bin", n.cluster = cores)

# cleanup
drops <- rownames(summary(sone_bin_noperm))[c(1,3,4,5,6,7,9,10,11,15,16,17,19,21,22)]
cross_ug <- drop.markers(cross_ug,drops)
noperms <- drop.markers(noperms,drops)

################################################################################
## Flip if the chr is backwards
direc <- sapply(c(1:24),function(i) {
 pos <- as.numeric(gsub(".*:","",markernames(noperms,i)))
 map <- as.numeric(pull.map(noperms)[[i]])
 cor(pos,map, use="complete.obs")
})

if(any(direc < 0)) noperms <- flip.order(noperms,which(direc < 0))
if(any(direc < 0)) cross_ug <- flip.order(cross_ug,which(direc < 0))
################################################################################

i <- 1 ; plotit(cross_ug,'_cross_ug')
i <- 2 ; plotit(noperms,'_noperms')

save.image('NBH_compare_ungenotyped.rsave')

################################################################################



## GENOTYPED DATASET ONLY (no stratified permutations needed)
erprob_noperm <- 0.05 #determined from mapping step
noperms <- sim.geno(noperms, n.draws = 256, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
noperms <- calc.genoprob(noperms, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")

## SCAN
## sone_bin_noperm <- scanone(noperms, pheno.col=4, method="em", model="bin", n.cluster = cores)
## sone_bin_noperm <- scanone(noperms, pheno.col=4, method="imp", model="bin", n.cluster = cores)

# SCAN
### no binary model for imp, uses em and transformed normal
sone_noperm_bin <- scanone(noperms, pheno.col=4, method="em",  model="bin")
sone_noperm_bin_np <- scanone(noperms, pheno.col=1, model="np", n.cluster = cores)

## Compare with HK (violate Lander and Botstein)
sone_noperm_bin_hk <- scanone(noperms, pheno.col=4, method="hk",  model="bin", n.cluster = cores)

##  Method imp not available for binary model; using nqrank transformed normal
sone_noperm_normal <- scanone(noperms, pheno.col=5, method="imp",  model="normal", n.cluster = cores)

# PERMUTE LOD CUTOFF
## Drop chromosome 2 and X for permutations
perm_noperms <- subset(noperms, chr = c(1,3:4,6:24))
#sone.perms_noperm <- scanone(perm_noperms, pheno.col=4, method="em",  model="bin", n.perm=1000, n.cluster = cores)
#sone.perms_noperm <- scanone(perm_noperms, pheno.col=4, method="imp",  model="bin", n.perm=1000, n.cluster = cores)

### no binary model for imp, uses em and transformed normal
sone.perms_noperm_bin <- scanone(perm_noperms, pheno.col=4, method="em",  model="bin", n.perm=10000)
lod_noperm_bin <- summary(sone.perms_noperm_bin)[[1]]
qtl.add_noperm_bin <- summary(sone_noperm_bin, lod_noperm_bin)

### use binary model for hk (lander and botstein violation)
sone.perms_noperm_bin_hk <- scanone(perm_noperms, pheno.col=4, method="hk",  model="bin", n.perm=10000)
lod_noperm_bin_hk <- summary(sone.perms_noperm_bin_hk)[[1]]
qtl.add_noperm_bin_hk <- summary(sone_noperm_bin_hk, lod_noperm_bin_hk)

##  Method imp not available for binary model; using nqrank transformed normal
sone.perms_noperm_normal <- scanone(perm_noperms, pheno.col=5, method="imp",  model="normal", n.perm=10000, n.cluster = cores)
lod_noperm_normal <- summary(sone.perms_noperm_normal)[[1]]
qtl.add_noperm_normal <- summary(sone_noperm_normal, lod_noperm_normal)

# FIT AND COMPARE MODELS
## make QTL object for the model (chr 2,8,18)
em.qtl.4_noperms_bin <- makeqtl(noperms, chr = qtl.add_noperm_bin[['chr']], pos = qtl.add_noperm_bin[['pos']], what="prob")
hk.qtl.4_noperms_bin <- makeqtl(noperms, chr = qtl.add_noperm_bin_hk[['chr']], pos = qtl.add_noperm_bin_hk[['pos']], what="prob")
imp.qtl.4_noperms_normal <- makeqtl(noperms, chr = qtl.add_noperm_normal[['chr']], pos = qtl.add_noperm_normal[['pos']], what="draws")

## refine the position of QTL considering all 4 QTL (only works with imp)
imp.qtl.4_noperms_normal <- refineqtl(noperms, pheno.col = 5, qtl = imp.qtl.4_noperms_normal, method = "imp", model='normal', incl.markers=T)
hk.qtl.4_noperms_bin <- refineqtl(noperms, pheno.col = 4, qtl = hk.qtl.4_noperms_bin, method = "hk", model='bin', incl.markers=T)

## Fit model (only works with imp or hk, hk not appropriate)
fit_imp_normal_noperm <- fitqtl(noperms, pheno.col=5, method="imp", model="normal", qtl = imp.qtl.4_noperms_normal,
                covar=NULL, formula = y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)
summary(fit_imp_normal_noperm)


fit_hk_bin_noperm <- fitqtl(noperms, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.4_noperms_bin,
                covar=NULL, formula = y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=100000, forceXcovar=FALSE)
summary(fit_hk_bin_noperm)

save.image('NBH_compare_ungenotyped.rsave')



## UNGENOTYPED DATASET ONLY
erprob_noperm <- 0.05 #determined from mapping step
cross_ug <- sim.geno(cross_ug, n.draws = 256, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
cross_ug <- calc.genoprob(cross_ug, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")

# SCAN
### no binary model for imp, uses em and transformed normal
sone_ug_bin <- scanone(cross_ug, pheno.col=4, method="em",  model="bin")

## Compare with HK (violate Lander and Botstein)
sone_ug_bin_hk <- scanone(cross_ug, pheno.col=4, method="hk",  model="bin", n.cluster = cores)

##  Method imp not available for binary model; using nqrank transformed normal
sone_ug_normal <- scanone(cross_ug, pheno.col=5, method="imp",  model="normal", n.cluster = cores)
sone_ug_bin_np <- scanone(cross_ug, pheno.col=1, model="np", n.cluster = cores)

# PERMUTE LOD CUTOFF
## Drop chromosome 2,X, and 18 for permutations
perm_cross_ug <- subset(cross_ug, chr = c(1,3:4,6:24))

### no binary model for imp, uses em and transformed normal
sone.perms_ug_bin <- scanone(perm_cross_ug, pheno.col=4, method="em",  model="bin", n.perm=10000, perm.strata = perm_cross_ug$pheno$GT_NG_ALT)
lod_ug_bin <- summary(sone.perms_ug_bin)[[1]]
qtl.add_ug_bin <- summary(sone_ug_bin, lod_ug_bin)

### us binary model for hk
sone.perms_ug_bin_hk <- scanone(perm_cross_ug , pheno.col=4, method="hk",  model="bin", n.perm=10000, n.cluster = cores, perm.strata = perm_cross_ug$pheno$GT_NG_ALT)
lod_ug_bin_hk <- summary(sone.perms_ug_bin_hk)[[1]]
qtl.add_ug_bin_hk <- summary(sone_ug_bin_hk, lod_ug_bin_hk)

##  Method imp not available for binary model; using nqrank transformed normal
sone.perms_ug_normal <- scanone(perm_cross_ug, pheno.col=5, method="imp",  model="normal", n.perm=10000, n.cluster = cores, perm.strata = perm_cross_ug$pheno$GT_NG_ALT)
lod_ug_normal <- summary(sone.perms_ug_normal)[[1]]
qtl.add_ug_normal <- summary(sone_ug_normal, lod_ug_normal)

# FIT AND COMPARE MODELS
## make QTL object for the model (chr 2,8,18)
em.qtl.4_cross_ug <- makeqtl(cross_ug, chr = qtl.add_ug_bin[['chr']], pos = qtl.add_ug_bin[['pos']], what="prob")
hk.qtl.4_cross_ug_bin <- makeqtl(cross_ug, chr = qtl.add_ug_bin_hk[['chr']], pos = qtl.add_ug_bin_hk[['pos']], what="prob")
imp.qtl.4_cross_ug <- makeqtl(cross_ug, chr = qtl.add_ug_normal[['chr']], pos = qtl.add_ug_normal[['pos']], what="draws")


## refine the position of QTL considering all 4 QTL
imp.qtl.4_cross_ug <- refineqtl(cross_ug, pheno.col = 5, qtl = imp.qtl.4_cross_ug, method = "imp", model='normal', incl.markers=T)
hk.qtl.4_cross_ug_bin <- refineqtl(cross_ug, pheno.col = 4, qtl = hk.qtl.4_cross_ug_bin, method = "hk", model='bin', incl.markers=T)

## Fit model
fit_imp_normal_ug <- fitqtl(cross_ug, pheno.col=5, method="imp", model="normal", qtl = imp.qtl.4_cross_ug,
                covar=NULL, formula = y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

fit_hk_bin_ug <- fitqtl(cross_ug, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.4_cross_ug_bin,
                covar=NULL, formula = y~Q1+Q2+Q3, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=100000, forceXcovar=FALSE)

summary(fit_imp_normal_ug)
summary(fit_hk_bin_ug)

save.image('NBH_compare_ungenotyped.rsave')

png('~/public_html/dist.png')
hist(cross_ug$pheno$pheno_norm)
dev.off()

#######################################################################################
#######################################################################################

  cim_noperms_notransf <-   cim(noperms, pheno.col=1, n.marcovar=3, window=500,
         method= "em",
         imp.method="argmax", error.prob=0.05,
         map.function= "kosambi")

  cim_noperms <- cim(noperms, pheno.col=4, n.marcovar=3, window=500,
         method= "em",
         imp.method="argmax", error.prob=0.05,
         map.function= "kosambi")

  cim_ug <-   cim(noperms, pheno.col=5, n.marcovar=3, window=500,
         method= "em",
         imp.method="argmax", error.prob=0.05,
         map.function= "kosambi")

png('~/public_html/cim.png', width = 2000)
plot(cim_noperms_notransf, cim_noperms, cim_ug, col=c("green","blue", "red"))
add.cim.covar(cim_noperms_notransf, col="green")
add.cim.covar(cim_noperms, col="blue")
add.cim.covar(cim_ug, col="red")
dev.off()


plotLodProfile(rqtl)


#### IMPUTATION: most similar between normal-binary

  cim_noperms_notransf <-   cim(noperms, pheno.col=1, n.marcovar=3, window=500,
         method= "imp",
         imp.method="imp", error.prob=0.05,
         map.function= "kosambi")

  cim_noperms <- cim(noperms, pheno.col=4, n.marcovar=3, window=500,
         method= "imp",
         imp.method="imp", error.prob=0.05,
         map.function= "kosambi")

  cim_ug <-   cim(noperms, pheno.col=5, n.marcovar=3, window=500,
         method= "imp",
         imp.method="imp", error.prob=0.05,
         map.function= "kosambi")

png('~/public_html/cim-imp.png', width = 2000)
plot(cim_noperms_notransf, cim_noperms, cim_ug, col=c("green","blue", "red"))
add.cim.covar(cim_noperms_notransf, col="green")
add.cim.covar(cim_noperms, col="blue")
add.cim.covar(cim_ug, col="red")
dev.off()



sone.perms_ug_bin <- scanone(perm_cross_ug, pheno.col=4, method="em",  model="bin", n.perm=10000, perm.strata = perm_cross_ug$pheno$GT_NG_ALT)



plt <- scanonevar(cross_ug)

png('~/public_html/scanvar.png', width = 2000)
plot(plt)
dev.off()

png('~/public_html/cim.png', width = 2000)
plot(cim_noperms_notransf, cim_noperms, cim_ug, col=c("green","blue", "red"))
add.cim.covar(cim_noperms_notransf, col="green")
add.cim.covar(cim_noperms, col="blue")
add.cim.covar(cim_ug, col="red")
dev.off()






  cim_noperms <-   cim(noperms, pheno.col=4, n.marcovar=3, window=10,
         method= "em",
         imp.method=c("imp", "argmax"), error.prob=0.05,
         map.function= "kosambi")


 can I make window high to select one QTL per chromosome. 

 add.cim.covar(out.cim, chr=c(1,4,6,15))



## Only options for ungenotyped is imp and hk. hk should be avoided, imp needs stratfied permutations

cross_ug

erprob_ug <- 0.05 #determined from mapping step
cross_ug <- sim.geno(cross_ug, n.draws = 150, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
cross_ug <- calc.genoprob(cross_ug, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")

## SCAN
sone_bin_ug <- scanone(cross_ug, pheno.col=4, method="hk", model="bin", n.cluster = cores)

## Drop chromosome 2,X, and 18 for permutations
perm_cross_ug <- subset(cross_ug, chr = c(1,3:4,6:17,19:24))
sone.perms_ug <- scanone(perm_cross_ug, pheno.col=4, method="hk",  model="bin", n.perm=1000, n.cluster = cores, perm.strata = perm_cross_ug$pheno$GT_NG_ALT)


## LOD 
lod_ug <- summary(sone.perms_ug)[[2]]
qtl.add_ug <- summary(sone_bin_ug, lod_noperm)

## make QTL object for the model
hk.qtl.4_cross_ug <- makeqtl(cross_ug, chr = qtl.add_noperm[['chr']], pos = qtl.add_noperm[['pos']], what = "prob")
## refine the position of QTL considering all 4 QTL
hk.qtl.4_cross_ug <- refineqtl(cross_ug, pheno.col = 4, qtl = hk.qtl.4_cross_ug, method = "hk", model='binary', incl.markers=T)

fit_hk_bin_noperm_ug <- fitqtl(cross_ug, pheno.col=4, method="hk", model="bin", qtl = hk.qtl.4_cross_ug,
                covar = NULL, formula = y~Q1+Q2+Q3+Q4, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)



save.image('NBH_ungenotyped.rsave')




sone_bin_noperm_em <- scanone(noperms, pheno.col=4, method="em", model="bin")
sone.perms_ug_em <- scanone(cross_ug, pheno.col=4, method="em",  model="bin")



sone_bin_noperm_em_np <- scanone(noperms, pheno.col=4, method="em", model="np")
sone.perms_ug_em_np <- scanone(cross_ug, pheno.col=4, method="em",  model="np")




i <- 100 ; plotit(cross_ug,'cross_ug')
i <- 101 ; plotit(noperms,'cross_noperms')



https://www.biostat.wisc.edu/~kbroman/teaching/misc/Jax/2003/qtlhandouts.pdf


Selective genotyping
- perform IM with all genotypes

Non-normal traits:
- non-parametric





cross_ug

erprob_ug <- 0.05 #determined from mapping step
cross_ug <- sim.geno(cross_ug, n.draws = 150, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")
cross_ug <- calc.genoprob(cross_ug, error.prob = erprob_noperm, map.function="kosambi", stepwidth="fixed")

## SCAN
#sone_bin_ug <- scanone(cross_ug, pheno.col=5, method="imp", model="bin", n.cluster = cores)
sone_bin_ug <- scanone(cross_ug, pheno.col=5, method="imp", model="normal", n.cluster = cores)
#sone_bin_ug <- scanone(cross_ug, pheno.col=1, model="np", n.cluster = cores)
#sone_bin_ug <- scanone(cross_ug, pheno.col=4, method="hk", model="bin", n.cluster = cores)

## Drop chromosome 2,X, and 18 for permutations
perm_cross_ug <- subset(cross_ug, chr = c(1,3:4,6:17,19:24))
sone.perms_ug <- scanone(perm_cross_ug, pheno.col=5, method="imp",  model="normal", n.perm=1000, n.cluster = cores, perm.strata = perm_cross_ug$pheno$GT_NG_ALT)


## LOD 
lod_ug <- summary(sone.perms_ug)[[2]]
qtl.add_ug <- summary(sone_bin_ug, lod_ug)

## make QTL object for the model
imp.qtl.4_cross_ug <- makeqtl(cross_ug, chr = qtl.add_ug [['chr']], pos = qtl.add_ug [['pos']], what = "prob")
## refine the position of QTL considering all 4 QTL
imp.qtl.4_cross_ug <- refineqtl(cross_ug, pheno.col = 5, qtl = imp.qtl.4_cross_ug, method = "imp", model='normal', incl.markers=T)

fit_imp_normal_noperm_ug <- fitqtl(cross_ug, pheno.col = 5, method="imp", model="normal", qtl = imp.qtl.4_cross_ug,
                covar = NULL, formula = y~Q1+Q2+Q3+Q4+Q5, dropone=TRUE, get.ests=T,
                run.checks=TRUE, tol=1e-4, maxit=10000, forceXcovar=FALSE)

