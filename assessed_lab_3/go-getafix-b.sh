#!/bin/bash
#SBATCH --partition=cosc
#SBATCH --job-name=hello4x1x1
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=0-2:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

DATE=$(date +"%Y%m%d%H%M")
echo "time started  "$DATE
echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "running with OMP_NUM_THREADS     = $OMP_NUM_THREADS "
echo "running with SLURM_TASKS_PER_NODE= $SLURM_TASKS_PER_NODE "
echo "running with SLURM_NTASKS total  = $SLURM_NTASKS "
echo "Now we start the show:"
export TIMEFORMAT="%E sec"

module load mpi/openmpi-x86_64
time mpirun -n ${SLURM_NTASKS} ./hello-world.exe

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE

