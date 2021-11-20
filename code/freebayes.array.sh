#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_err_%A_%a.txt
#SBATCH --array=1-8577
#SBATCH -p med
#SBATCH --mem=30000
#SBATCH -t 48:00:00
###### number of nodes
###### number of processors
#SBATCH --cpus-per-task=6

## Set environmental variables from file
export $(grep -v '^#' code/analysis.env | xargs)

$my_freebayes --version
$my_bedtools --version


scafnum=$(expr $SLURM_ARRAY_TASK_ID + -1)
scaf=Scaffold$scafnum
endpos=$(expr $(grep -P "$scaf\t" ${bwagenind}.fai | cut -f 2) - 1) || echo "ref index missing"
region=$scaf:1..$endpos

echo $region

outfile=$scaf.vcf

$my_bamtools merge -list $bam_list -region $region| \
$my_bamtools filter -in stdin -mapQuality '>30' -isProperPair true | \
$my_bedtools intersect -sorted -a stdin -b $bed_regions | \
$my_freebayes -f $bwagenind --use-best-n-alleles 4 --pooled-discrete --stdin > $vcf_out/$outfile && \
echo 'Done' || echo 'pipe failed'

date 

