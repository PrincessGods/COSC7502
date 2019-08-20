#!/bin/sh
#SBATCH --time=00:02:30
#SBATCH --partition=coursework
#SBATCH --job-name=ga_island
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "running with SLURM_TASKS_PER_NODE= $SLURM_TASKS_PER_NODE "
echo "running with SLURM_NPROCS= $SLURM_NPROCS "

module load mpi
time mpirun -n ${SLURM_NPROCS} ./ga_island.exe 16 100 1
