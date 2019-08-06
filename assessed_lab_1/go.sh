#!/bin/bash
#SBATCH --time=0:5:00
#SBATCH --job-name=openmptest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=coursework

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

DATE=$(date +"%Y%m%d%H%M")
echo "time started  "$DATE

time ./openmp-test.exe

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
