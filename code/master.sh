





### Demultiplex
### /home/jmiller1/QTL_Map_Raw/QTL_scripts

BarcodeSplitListBestRadPairedEnd.pl  BP_run_BestRadSplit14.sh  demulti.sh  ER_run_BestRadSplit14.sh  NBH_run_BestRadSplit14.sh  NEW_run_BestRadSplit14.sh

### Alignment 
### code/array.align.sh <directory_input_bams> <output_dir>
code/array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/NBH /home/jmiller1/QTL_Map_Raw/align
code/array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/ELR /home/jmiller1/QTL_Map_Raw/align
code/array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/BRP /home/jmiller1/QTL_Map_Raw/align
code/array.align.sh /home/jmiller1/QTL_Map_Raw/demulti_out/NEW /home/jmiller1/QTL_Map_Raw/align

### bams stored at 
/home/jmiller1/compressed.data/qtl.data/align



#### Plink2Rqtl.JM.R <inputdir> <population> <output.tag>
Rscript --vanilla Plink2Rqtl.JM.R popgen/plinkfiles/ind.pops NBH '.unphased.f2.csvr'
