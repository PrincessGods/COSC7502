/* Copyleft  Bromley brom@physics.uq.edu.au August 2017
 * Code to test out OpenMP/gprof as per lecture notes
 *
 * basic compile with
 * gcc -Wall -Wextra openmp-test.c -fopenmp -o openmp-test.exe
 * run after setting, eg.
 * export OMP_NUM_THREADS=4
 *
 * gprof compile with
 * gcc -pg -O0 openmp-test.c -fopenmp  -o openmp-test.exe
 * after run
 * gprof openmp-test.exe gmon.out > gmon-analysis.txt
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>

int main(){
  int nt, tid;
  time_t epochsecs;

  nt = omp_get_max_threads(); // not omp_get_num_threads();
  /* nt is shared by all threads */
  #pragma omp parallel private(tid,epochsecs) shared(nt)
  {
     /* value of tid is set by each thread */
     tid = omp_get_thread_num();
     /* each thread will tend to call at the same time*/
     epochsecs = time(NULL);
     printf("Hello from thread %d of %d at time %ld \n",tid,nt,(long int)epochsecs);
  }
  /* now value of tid is undefined */
  exit(EXIT_SUCCESS);
}

