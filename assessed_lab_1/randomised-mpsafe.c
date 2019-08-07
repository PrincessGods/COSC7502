/*  Copyleft Bromley brom@physics.uq.edu.au August 2018
 *  Implementing a couple of the simplest random methods in c

 *  basic compile eg. with:
 *  gcc -Wall -Wextra randomised-simple.c -o randomised-simple.exe -lm
 *  then test usage eg.:
 *  ./randomised-simple.exe 12345 10
 *  
 *  older gnu compilers might need -std=gnu99 for drand48 support
 *
 *  ************80 char length page width best for a2ps ***********************
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(int argc, char *argv[]) {
  unsigned int iseed;  /* input 1 random number seed */
  int nloops;          /* input 2 number of loops */

  int jloops;

  //int irandcur;
  //double drandcur;

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <iseed> <nloops>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  iseed = atoi(argv[1]);
  nloops = atoi(argv[2]);

  printf("Running %s with iseed=%u and nloops=%d \n", argv[0],iseed,nloops);

/* method 1 rand() generation.... see ** man rand **
 * RAND_MAX is defined in stdlib.h
 * First generate the seed and then run off nloops*/
//   printf("method1 using rand() with RAND_MAX=%d \n", RAND_MAX);
//   srand(iseed);
//   for (jloops = 0; jloops < nloops; jloops++) {
//     irandcur =  rand();
//     drandcur = (double)rand() / ((double)RAND_MAX+1.0);
//     printf("method1 int=%d, dble=%f, \n", irandcur,drandcur);
//   }

//   printf("\n");

/* method 2 drand48() generation (uses long int seed)
 * ** man drand48 **
 * for more information about (lack of) thread safety */
  double x, drandcur, sum, result;
  long int iseedlong;
  int num;

  printf("method2 using drand48() \n");
  
  sum = 0;
  #pragma omp parallel private(num, x,iseedlong)
  {
    num = omp_get_thread_num();
    iseedlong = (long int) iseed + num;

    struct drand48_data* randBuffer = malloc(sizeof(struct drand48_data));
    srand48_r(iseedlong, randBuffer);
    
    #pragma omp for
    for (jloops = 0; jloops < nloops; jloops++) {
        drand48_r(randBuffer, &x); // random number between 0,1
        drandcur = pow(x, 2);
        //printf("method2 drandcur=%g, num=%d, jloops=%d\n", drandcur, num, jloops);
        #pragma omp critical
        {   
            sum += drandcur;
        }
    }
    free(randBuffer); 
  }

  result = sum / nloops;

  printf("method2 result=%g\n",result);

  exit(EXIT_SUCCESS);
}
