#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_err_%A_%a.txt
#SBATCH --array=1-96
#SBATCH -p hi
#SBATCH --mem=30000
#SBATCH --cpus-per-task=6
###### number of nodes
###### number of processors
####SBATCH -n 1

### Set on command line
indir=/

outdir=/home/jmiller1/QTL_Map_Raw/align/


lib=$(echo $fq1 | sed 's/_/\t/'| cut -f 1)
fq1=$(find $indir -name "*RA*fastq.gz" | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)
fq2=$(echo $fq1 | sed 's/RA/RB/')
root=$(echo $fq1 | sed 's/.*\///' | sed 's/_R._/_/' | cut -c 1-8,11-18)
rg=$(echo \@RG\\tID:$root\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$root)
tempsort=$root.temp
outfile=$outdir/$root.bam

echo $root
echo $fq1
echo $fq2
echo $rg
echo $tempsort
echo $outfile

## Alignment 
$my_bwa mem $bwagenind -t 6 -R $rg $fq1 $fq2 | \
    $my_samblstr | \
    $my_samtools view -S -h -u - | \
    $my_samtools sort - -T /scratch/$tempsort -O bam -o $outfile

$my_samtools index $outfile
