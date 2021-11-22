#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o code/slurm_out/array_job_out_%A_%a.txt
#SBATCH -e code/slurm_out/array_job_err_%A_%a.txt
#SBATCH --array=1-8577
#SBATCH -p med
#SBATCH --mem=30g
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=6

## Set environmental variables from file
export $(grep -v '^#' code/seq/analysis.env | xargs)

set -o allexport
source code/seq/analysis.env
set +o allexport


## Echo version to code/seq/array_job_out
$my_freebayes --version
$my_bedtools --version
$my_bamtools --version

scaf=$(sed -n "$SLURM_ARRAY_TASK_ID p" $all_scaf | cut -f1)
endpos=$(sed -n "$SLURM_ARRAY_TASK_ID p" $all_scaf | cut -f2)

region=$scaf:1..$endpos
echo $er_out_dir/$region

outfile=$scaf.vcf

$my_bamtools merge -list $bam_list -region $region| \
$my_bamtools filter -in stdin -mapQuality '>30' -isProperPair true | \
$my_bedtools intersect -sorted -a stdin -b $bed_regions | \
$my_freebayes -f $bwagenind --use-best-n-alleles 4 --pooled-discrete --stdin > $vcf_out/$outfile && \
echo 'Done' || echo 'pipe failed'

date 

