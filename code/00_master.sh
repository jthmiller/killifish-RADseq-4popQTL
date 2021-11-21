





### Demultiplex
### /home/jmiller1/QTL_Map_Raw/QTL_scripts
### Depends: BarcodeSplitListBestRadPairedEnd.pl  (author: Mike Miller)
sbatch code/01_BP_run_BestRadSplit14.sh  
sbatch code/01_ER_run_BestRadSplit14.sh  
sbatch code/01_NBH_run_BestRadSplit14.sh  
sbatch code/01_NEW_run_BestRadSplit14.sh

### Alignment 
### code/array.align.sh <directory_input_bams> <output_dir>
### Depends: analysis.env
sbatch code/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/NBH /home/jmiller1/QTL_Map_Raw/align
sbatch code/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/ELR /home/jmiller1/QTL_Map_Raw/align
sbatch code/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/BRP /home/jmiller1/QTL_Map_Raw/align
sbatch code/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/NEW /home/jmiller1/QTL_Map_Raw/align
### bams currently stored at /home/jmiller1/compressed.data/qtl.data/align

### Call Variants
sbatch code/03_freebayes.array.sh 

## Concat chromosomal vcfs
module load vcftools
dir=/home/jmiller1/QTL_Map_Raw/vcf/freebayes.array
vcf-concat $dir/chr*.vcf | gzip > $dir/SOMM.vcf.gz




### Convert PLINK format to rQTL format
### Plink2Rqtl.JM.R <inputdir> <population> <output.tag>
### Depends: PLINK2RQTL.f2.R (modified from - Brockmann group - HU Berlin, Danny Arends)
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops NBH '.unphased.f2.csvr'
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops ELR '.unphased.f2.csvr'
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops NEW '.unphased.f2.csvr'
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops BRP '.unphased.f2.csvr'
