#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o code/out_er/array_job_out_%A_%a.txt
#SBATCH -e code/out_er/array_job_err_%A_%a.txt
#SBATCH --array=1-96
#SBATCH -p hi
#SBATCH --mem=30000
#SBATCH --cpus-per-task=6
####SBATCH -n 1

### Set on command line
indir=${1}
outdir=${2}

fq1=$(find $indir -name "*RA*fastq.gz" | sed -n $(echo $SLURM_ARRAY_TASK_ID)p)
fq2=$(echo $fq1 | sed 's/RA/RB/')

lib=$(echo $fq1 | sed 's/_/\t/'| cut -f 1)
root=$(echo $fq1 | sed 's/.*\///' | sed 's/_R._/_/' | cut -c 1-8,11-18)

### Set the read group information
rg=$(echo \@RG\\tID:$root\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$root)

### Temp dir for sorting the bam
tempsort=$root.temp

### Write the bam file
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
    $my_samtools sort - -T /scratch/$tempsort -O bam -o $outfile \
&& $my_samtools index $outfile \
|| echo 'alignment failure'
