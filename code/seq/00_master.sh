### Workfolw for sequence analysis 

### Path to data at farm.cse.ucdavis.edu: 
/home/jmiller1/compressed.data/qtl.data/titan.bch.msu.edu/20160211_RADSeq_PE/

### Md5sums
74e313d89665dec8f7655b56fe010828  SOMM075_NBH_AGTCAA_L005_R1_001.fastq.gz
511f15f1ad16ef0e8a1566d515dbd4cd  SOMM075_NBH_AGTCAA_L005_R2_001.fastq.gz
6aa66917e7a9d860b57e4ae2ff98d25f  SOMM086_New_CCGAGA_L006_R1_001.fastq.gz
85cb8dac0aba2b149c535d0a8c6ed9c5  SOMM086_New_CCGAGA_L006_R2_001.fastq.gz
f9d0eae392628364ae56190f4596fb9c  SOMM087_BP_AGAAGA_L007_R1_001.fastq.gz
a534c97f23a062d6cb4af6ab06571ade  SOMM087_BP_AGAAGA_L007_R2_001.fastq.gz
e16cffe94c9e65cbbcf087ed5b935627  SOMM088_ER_GGTTCA_L008_R1_001.fastq.gz
59f0482098525c5e6f8ad6a83f99e1ac  SOMM088_ER_GGTTCA_L008_R2_001.fastq.gz

### Reference genome: 
/home/jmiller1/genomes_jm/mapped/ALLMAPS_OUT/unsplit_merge.fasta
Chromosome assembled genome from Miller et al 2018. 

### Demultiplex
### Depends: BarcodeSplitListBestRadPairedEnd.pl  (author: Mike Miller)
sbatch code/seq/01_BP_run_BestRadSplit14.sh  
sbatch code/seq/01_ER_run_BestRadSplit14.sh  
sbatch code/seq/01_NBH_run_BestRadSplit14.sh  
sbatch code/seq/01_NEW_run_BestRadSplit14.sh

### Alignment 
### code/seq/array.align.sh <directory_input_bams> <output_dir>
### Depends: analysis.env
sbatch code/seq/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/NBH /home/jmiller1/QTL_Map_Raw/align
sbatch code/seq/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/ELR /home/jmiller1/QTL_Map_Raw/align
sbatch code/seq/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/BRP /home/jmiller1/QTL_Map_Raw/align
sbatch code/seq/02_array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/NEW /home/jmiller1/QTL_Map_Raw/align
### bams currently stored at /home/jmiller1/compressed.data/qtl.data/align

### Call joint variants
sbatch code/seq/03_freebayes.array.sh 

## Concat chromosomal vcfs
module load vcftools/0.1.13
vcf-concat $dir/chr*.vcf | gzip > $dir/SOMM.vcf.gz

### Convert PLINK format to rQTL format
### Plink2Rqtl.JM.R <inputdir> <population> <output.tag>
### Depends: PLINK2RQTL.f2.R (modified from - Brockmann group - HU Berlin, Danny Arends)
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops NBH '.unphased.f2.csvr'
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops ELR '.unphased.f2.csvr'
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops NEW '.unphased.f2.csvr'
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops BRP '.unphased.f2.csvr'
