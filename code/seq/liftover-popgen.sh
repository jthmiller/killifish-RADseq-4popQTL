#!/bin/bash

# Copy over files with R
Rscript liftover-popgen.r 

# Liftover as bedfiles

# useage:  liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed

liftover='/home/jmiller1/bin/liftOver'
chainfile='/home/jmiller1/bin/code/ALLMAPS_OUT/unsplit_merge.chain'
indir='/home/jmiller1/QTL_Map_Raw/popgen/outlier_table'
outdir='/home/jmiller1/QTL_Map_Raw/popgen/outfiles'

### 
file='/home/jmiller1/QTL_Map_Raw/popgen/tables/commas.txt'
tr '/r' '/n' < $file > commmas.fix.txt

file='/home/jmiller1/QTL_Map_Raw/popgen/tables/commas.fix.txt'
head -n1 $file > genemod.ncbi.header
sed -i '1d' $file

### liftover bed coords
$liftover -bedPlus=1 $file $chainfile $file.lifted $file.unmaped

### Add header to intermediat file
cat /home/jmiller1/QTL_Map_Raw/popgen/tables/genemod.ncbi.header $file > commas.txt 

### Copy outlier tables 
cp /home/nreid/popgen/outliers*.txt $indir

## awk 'FNR==NR                   #### Checking condition of FNR==NR which will be TRUE when first Input_file1 is being read.
## {A[$2]=$1;                     #### create an array named A with index of field  2 and have value as first field.
## next}                          #### using next keyword(built-in awk) so it will skip all next statements.
## ($1 in A)                      #### Now checking if first field of file2 is present in array A, this will be checked only when Input_file2 is being read.
## {print A[$1], $2               #### printing value of array A's value whose index is $1 and $2 of Input_file2.
## }' Input_file1  Input_file2    #### Mentioning the Input_file1 and Input_file2 here.

## Liftover files of AHR and candidate outliers

FILES=$indir/*.txt
scaf="$indir/scaf"
ncbi="$indir/ncbi"

for file in $FILES
do
    sed -i 's|["'\'']||g' $file
    echo "striping header to $file.ncbi.header"
    head -n1 $file > $file.ncbi.header
    sed -i '1d' $file
    awk -v OFS='\t' 'FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' convert.table $file > $file.ncbi
    \rm $file 
done

converted=$indir/*.txt.ncbi

for file in $converted ; do
    echo "lifting over $file"
    $liftover -bedPlus=1 $file $chainfile $file.lifted $file.unmaped
    cat $file.header $file.lifted > $file.temp
    mv $file.temp $file.lifted
    \rm $file.header
    \rm $file
done
