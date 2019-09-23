/*****************************************************************************
* FILE: hello-world-mpi-openmp-hybrid.c
* DESCRIPTION:
*   This simple program is inspired by the more complicated
*      - omp_dotprod_hybrid.c  - Hybrid MPI and OpenMP version
*  https://computing.llnl.gov/tutorials/openMP/samples/C/omp_dotprod_hybrid.c
*
* dotprod SOURCE: Blaise Barney
* hello world REVISED:  19/10/17 Michael Bromley
*
* Compile eg. with
* module load mpi/openmpi-x86_64
* mpicc -fopenmp -o hello-world.exe hello-world-mpi-openmp-hybrid.c
* sbatch goa.sh
* sbatch gob.sh
* sbatch goc.sh
* sbatch god.sh
* sbatch goe.sh
******************************************************************************/

#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main (int argc, char* argv[]){
  int myid, tid, numprocs, nthreads, ntotthreads;
  int i,mychecksum,myallsum;

/* MPI Initialization */
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);

  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);

  if (myid == 0){
    printf("Starting - %d MPI tasks - master host %s \n",numprocs,hostname);
  }

  #pragma omp parallel private(tid) shared(nthreads,myid,hostname)
  {
    tid = omp_get_thread_num();
    if (tid ==0){
      nthreads = omp_get_num_threads();
      printf("Task %d has %d OpenMP threads on host %s \n",myid,nthreads,hostname);
    }
    #pragma omp critical
    {
      printf("Hello World from MPI Task %d and OpenMP thread %d on host %s \n",myid,tid,hostname);
    }
  }

/* finally perform an MPI summation of total number of threads over all tasks */
  ntotthreads = 0;
  MPI_Reduce (&nthreads, &ntotthreads, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if (myid == 0){
    printf ("Done. Hybrid version: global number of threads = %d \n", ntotthreads);
  }
  MPI_Finalize();
}

