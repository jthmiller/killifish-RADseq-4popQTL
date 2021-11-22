#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o code/seq/slurm_out/array_job_out_%A_%a.txt
#SBATCH -e code/seq/slurm_out/array_job_err_%A_%a.txt
#SBATCH --array=1-8577
#SBATCH -p med
#SBATCH --mem=30000
#SBATCH -t 48:00:00
#SBATCH --cpus-per-task=2
