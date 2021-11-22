#!/bin/bash

module load vcftools/0.1.13
vcftools --version

set -o allexport
source code/seq/analysis.env
set +o allexport

### Output allele counts to explore data in R
for X in NBH BRP NEW ELR
do
vcftools --gzvcf $vcfdir/SOMM.vcf.gz --keep $infiles/$X.samples --remove-filtered-all --recode --stdout \
| vcftools --vcf - --max-meanDP 90 --maxDP 90 --stdout --minGQ 20 $maf --remove-filtered-all --recode \
| vcf-sort \
| vcftools --vcf - --counts > ~/QTL_agri/data/$X.counts
done

### Output site mean depth to explore data in R
for X in NBH BRP NEW ELR
do
vcftools --gzvcf $vcfdir/SOMM.vcf.gz --keep $infiles/$X.samples --remove-filtered-all --recode --stdout \
| vcftools --vcf - --max-meanDP 90 --maxDP 90 --stdout --minGQ 20 $maf --remove-filtered-all --recode --stdout \
| vcf-sort \
| vcftools --vcf - --site-mean-depth > ~/QTL_agri/data/$X.site-mean-depth
done

### Filter
for X in NBH BRP NEW ELR
do
vcftools --gzvcf $vcfdir/SOMM.vcf.gz --keep $infiles/$X.samples --remove-filtered-all --recode --stdout \
| vcftools --vcf - --max-meanDP 90 --maxDP 90 --stdout --minGQ 20 $maf --remove-filtered-all --recode --stdout \
| gzip -c > $vcfdir/$X.vcf.gz> ~/QTL_agri/data/$X.vcf
done

## Use PLINK to seperate familes for rQTL. Alignment quality drops off for loci that have mean depth over 90. 
for X in NBH BRP NEW
do
	$plink --vcf $vcfdir/$X.vcf.gz --out $indpops/$X  $init_flagset --pheno $pheno --all-pheno --update-ids $infiles/SOMM.txt
	$plink --bfile $indpops/$X --out $indpops/$X $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names $X --make-founders
	$plink --bfile $indpops/$X --out $indpops/$X $flagset --pheno $pheno --all-pheno --keep-cluster-names $X $geno $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders
done

### REFILTER ELR. It contains the parent 'BLI' code (in addition to ELR fam)
$plink --vcf $vcfdir/SOMM.vcf.gz --out $indpops/$X  $init_flagset --make-bed --pheno $pheno --all-pheno --update-ids $infiles/SOMM.txt
$plink --bfile $indpops/$X --out $indpops/$X $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names ELR BLI --make-founders
$plink --bfile $indpops/$X --out $indpops/$X $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names ELR BLI $geno $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders
#$plink --bfile $indpops/$X --out $indpops/$X.prune $flagset --pheno $pheno --all-pheno --make-bed --keep-cluster-names ELR BLI --geno 0.9 $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders --indep 50 5 2

#### ELR IBD Analysis
$plink --bfile $indpops/$X.prune --out $indpops/$X --extract $indpops/$X.prune.in $flagset --pheno $pheno --all-pheno --keep-cluster-names ELR BLI $geno $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders --distance square bin 1-ibs ibs

##### ELR NBH IBD Analysis
$plink --vcf $vcfdir/SOMM.vcf.gz --out $indpops/NBH_ELR_NEW  $init_flagset --pheno $pheno --all-pheno --update-ids $infiles/SOMM.txt
$plink --bfile $indpops/NBH_ELR_NEW --out $indpops/NBH_ELR_NEW $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names NEW NBH ELR BLI --make-founders --geno 0.9 --recode --biallelic-only strict --snps-only just-acgt --nonfounders
$plink --bfile $indpops/NBH_ELR_NEW --out $indpops/NBH_ELR_NEW $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names NEW NBH ELR BLI --geno 0.9 --recode --biallelic-only strict --snps-only just-acgt --nonfounders --indep 50 5 2
$plink --bfile $indpops/NBH_ELR_NEW --out $indpops/NBH_ELR_NEW --extract $indpops/NBH_ELR_NEW.prune.in $flagset --pheno $pheno --all-pheno --keep-cluster-names NEW NBH ELR BLI --geno 0.9 --mind .9 --recode --snps-only just-acgt --nonfounders --distance square bin 1-ibs ibs

#### IBD Analysis for all pops
#$plink --bfile $indpops/$X --out $indpops/$X $flagset --pheno $pheno --all-pheno --keep-cluster-names $X $geno $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders
$plink --vcf $vcfdir/SOMM.vcf.gz --out $indpops/ALL  $init_flagset --pheno $pheno --all-pheno --update-ids $infiles/SOMM.txt
$plink --bfile $indpops/ALL --out $indpops/ALL $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names NEW NBH ELR BLI BRP --make-founders --geno 0.75 --recode --biallelic-only strict --snps-only just-acgt --nonfounders
$plink --bfile $indpops/ALL --out $indpops/ALL $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names NEW NBH ELR BLI BRP --geno 0.75 --recode --biallelic-only strict --snps-only just-acgt --nonfounders --indep 50 5 2
$plink --bfile $indpops/ALL --out $indpops/ALL --extract $indpops/ALL.prune.in $flagset --pheno $pheno --all-pheno --keep-cluster-names NEW NBH ELR BLI BRP --geno 0.3 --mind .9 --recode --snps-only just-acgt --nonfounders --distance square bin 1-ibs ibs flat-missing

### Parents each chromosome
for i in {1..24}
do
	$plink --bfile $indpops/ALL --out $indpops/$i.parents --make-bed --chr $i --set-missing-var-ids @:# --allow-extra-chr --autosome-num 24 --allow-no-sex --pheno $pheno --all-pheno --keep $infiles/parents.txt --keep-cluster-names NEW NBH ELR BLI --recode --family --snps-only just-acgt --nonfounders --distance square bin 1-ibs ibs flat-missing
done

