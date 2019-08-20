#!/bin/sh
#SBATCH --time=00:02:30
#SBATCH --partition=coursework
#SBATCH --job-name=ga_serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

time ./ga_serial.exe 6 5 1 -v
