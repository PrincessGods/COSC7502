/* Original Code: darcy_mpi.c
 * Parallel Darcy flow simulation by Lutz Gross
 *
 * Moified imagtime-mpi-double.c by Michael Bromley 2017-19
 *
 * Solves the Schrodinger Equation in
 * negative imaginary time using basic MPI parallelisation.
 *
 * Uses explicit method with purely double calculations.
 *
 * Input parameters are hardcoded in main.
 *
 * NOTE that nranks is also hardcoded!
 *
 * To compile on clusters :  
 * module load mpi/openmpi-x86_64
 * mpicc -Wall -fopenmp -Wextra -o imagtime-mpi-double.exe imagtime-mpi-double.c -lm
 * mpicc -O3 -o imagtime-mpi-double.exe imagtime-mpi-double.c -lm
 *
 * to run in sequential (match the 1 with whatever nranks is)
 * time mpirun -n 1 ./imagtime-mpi-double.exe
 * To run on clusters, check setup/queue/modify/run:
 * sbatch go.sh
 *
 */
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

/* global variables */
const MPI_Comm comm=MPI_COMM_WORLD;
int myrank, mysize;

void updatePsi(double *psinew, const double *psiold, const double *vpot,
                    const int N0, const int N1,
                    const double a, const double b, const double c)
{
    int i,j;
    for (j=1; j<N1-1; ++j) {
        for (i=1; i<N0-1; ++i) {
            psinew[i+N0*j] =
                      a * psiold[i+N0*j]      +   /* center */
                      b *( psiold[i-1+N0*j]   +   /* west */
                           psiold[i+1+N0*j]   +   /* east */
                           psiold[i+N0*(j-1)] +   /* south */
                           psiold[i+N0*(j+1)] ) +  /* north */
                      c * vpot[i+N0*j]*psiold[i+N0*j]; /* potential */
        }
    }
}

void updateNorthPsi(double *psinew, 
                         const double *psiold, const double *psinorth, const double *vpot,
                         const int N0, const int N1,
                         const double a, const double b, const double c)
{
    int i; /* j = N1-1 */
    for (i=1; i<N0-1; ++i) {
        psinew[i+N0*(N1-1)] =
                  a * psiold[i+N0*(N1-1)]    +   /* center */
                  b *( psiold[i-1+N0*(N1-1)] +   /* west */
                       psiold[i+1+N0*(N1-1)] +   /* east */
                       psiold[i+N0*(N1-2)]   +   /* south */
                       psinorth[i] )         +   /* north */
                  c * vpot[i+N0*(N1-1)]*psiold[i+N0*(N1-1)]; /* potential */
    }
}

void updateSouthPsi(double *psinew, 
                    const double *psiold, const double *psisouth, const double *vpot,
                    const int N0, const int N1,
                    const double a, const double b, const double c)
{
    int i; /* j = 0 */
    for (i=1; i<N0-1; ++i) {
        psinew[i] = a * psiold[i]      +  /* center */
                  b *( psiold[i-1]     +  /* west */
                       psiold[i+1]     +  /* east */
                       psisouth[i]     +  /* south */
                       psiold[i+N0] )  + /* north */
                  c * vpot[i]*psiold[i]; /* potential */
    }
}

void copyPsi(double *p_dst, const double *p_src,
             const int N0,  const int N1)
{
    int i,j;
    for (j=0; j<N1; ++j) {
        for (i=0; i<N0; ++i) {
            p_dst[i+N0*j] = p_src[i+N0*j];
        }
    }
}

/* currently only initialises Psi at one point in middle of grid */
void initialisePsi(double *pinit, const int N0, const int N1,
                   const double xymin, const double dxy)
{
    int i, j, jactual;
    double xactual, yactual;
    
    printf("initialising Psi: %d %d \n",myrank,mysize);

    for (j=0; j<N1; ++j) {
        for (i=0; i<N0; ++i) {
           pinit[i+N0*j]   = 0.0e0;  // initalise including boundary
        }
    }

/*    for (j=1; j<(N1-1); ++j) { */
    for (j=0; j<N1; ++j) {
        jactual = j + N1*myrank;
        if( jactual>0 && jactual<(N0-1) ) {  // ie. do not include boundaries
            yactual = xymin + ((double)jactual)*dxy;
            for (i=1; i<(N0-1); ++i) {  // do not include boundaries
                xactual = xymin + ((double)i)*dxy;
/*                pinit[i+N0*j] = 1.0e0;    // being boring */
                pinit[i+N0*j] = exp(-0.5e0*(xactual*xactual+yactual*yactual));  // test
/*              if( i == 32 ) {
                if( jactual == 32 ) { 
                    printf("debug: initialising Psi: %d %d %d %e %e %e \n",
                       myrank,i,jactual,xactual,yactual,pinit[i+N0*j]);
                } */
            }
        }
    }
}

void rescalePsi(double *pinit, const int N0, const int N1, const double psiscale)
{
    int i, j;
    for (j=0; j<N1; ++j) {
        for (i=0; i<N0; ++i) {
           pinit[i+N0*j] = pinit[i+N0*j] * psiscale; // scale including boundary
        }
    }
}

/* output all psi to grid */
void printPsi(double *psiout, const int N0, const int N1,
                      const double xymin, const double dxy)
{
    int i,j,jactual;
    double xactual, yactual;

    for (j=0; j<N1; ++j) {
        jactual = j + N1*myrank;
        yactual = xymin + ((double)jactual)*dxy;
        for (i=0; i<N0; ++i) {
/*            if ( psiout[i+N0*j] > 0 ){ */
            xactual = xymin + ((double)i)*dxy;
            printf("%d %d %d %e %e %e \n",
                             myrank,i,jactual,xactual,yactual,psiout[i+N0*j]);
/*            } */
        }
    }
}

/* initialises potential across local grid */
void initialisePot(double *vpot, const int N0, const int N1,
                   const double xymin, const double dxy)
{
    int i, j, jactual;
    double xactual, yactual, ractual;
    
    printf("initialising Psi: %d %d \n",myrank,mysize);

    for (j=0; j<N1; ++j) {
        for (i=0; i<N0; ++i) {
           vpot[i+N0*j]   = 0.0e0;  // initalise including boundary
        }
    }

    for (j=1; j<(N1-1); ++j) {
        jactual = j + N1*myrank;
        yactual = xymin + ((double)jactual)*dxy;
        for (i=1; i<(N0-1); ++i) {
            xactual = xymin + ((double)i)*dxy;
            ractual = sqrt(xactual*xactual+yactual*yactual);
            vpot[i+N0*j] = 0.5e0*ractual*ractual;
/*            printf("debug: initialising Vpot: %d %d %d %e %e %e \n",
                       myrank,i,jactual,xactual,yactual,vpot[i+N0*j]); */
        }
    }
}

/* calculate the kinetic energy integral of psi over the domain */
/* note that this routine ignores the north/south contribution */
/* thus only approximates the total kinetic energy - ok for now */
double getTenergy(double *psi, const int N0, const int N1,
                  const double pkinetic, const double dxy){
    double psitint=0.0e0;
    int i, j;
    #pragma omp parallel private(i, j)
    {   
        #pragma omp for reduction (+: psitint)
        for (j=1; j<(N1-1); ++j) {
            for (i=1; i<(N0-1); ++i) {
                psitint = psitint + pkinetic*psi[i+N0*j]*
                        ( -4.0e0 * psi[i+N0*j] + /* center */
                            psi[(i-1)+N0*j] +   /* west */
                            psi[(i+1)+N0*j] +   /* east */
                            psi[i+N0*(j-1)] +   /* south */
                            psi[i+N0*(j+1)] );   /* north */
            }
        }
    }

    psitint = psitint*dxy*dxy;
        /* enter your MPI code here */
    double allsum;
    MPI_Allreduce (&psitint, &allsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return allsum;
}

/* calculate the potential energy integral of psi over the domain */
double getVenergy(double *psi, const int N0, const int N1,
                  double *vpot, const double dxy){
    double psivint=0.0e0;
    int i,j;
    #pragma omp parallel private(i, j)
    {   
        
        #pragma omp for reduction (+: psivint)
        for (j=1; j<(N1-1); ++j) {
            for (i=1; i<(N0-1); ++i) {
                psivint = psivint + vpot[i+N0*j]*psi[i+N0*j]*psi[i+N0*j];
            }
    /*        printf("debug int j: %d %d %d %e \n",myrank,i,j,psivint); */
        }
    }

    psivint = psivint*dxy*dxy;
        /* enter your MPI code here */
    double allsum;
    MPI_Allreduce (&psivint, &allsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return allsum;
}

/* calculate the integral of psi^2 over the domain */
double getPsi2Integral(double *psi, const int N0, const int N1,
                            const double dxy) {
    double psi2int=0.0e0;
    int i,j;
    #pragma omp parallel private(i, j)
    {   
        
        #pragma omp for reduction (+: psi2int)
        for (j=1; j<(N1-1); ++j) {
            for (i=1; i<(N0-1); ++i) {
                psi2int = psi2int + psi[i+N0*j]*psi[i+N0*j];
            }
        }
    }
    
    psi2int = psi2int*dxy*dxy;
    
    /* enter your MPI code here - code will not work without this! */
    double allsum;
    MPI_Allreduce (&psi2int, &allsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return allsum;
}

int main(int argc, char** argv)
{
    const int N0=atoi(argv[1]);          /* hardcoded grid size max N0<2048*/
    const int nranks=atoi(argv[2]);       /* number of splits to do = MPI procs */
    const int N1=N0/nranks;   /* the other dimension of the rectangle */
    const double mass=1.0e0;  /* mass in arbitary units */
    const double hbar=1.0e0;  /* plancks constant in arbitary units */
    const double dtau=0.01e0; /* the (imaginary) time step size */
    const double tend=0.20e0; /* end (imaginary) time */
    const double xymax=5.0e0; /* boundary max location */
    const double xymin=-5.0e0;/* boundary min location */
    const double dxy=(xymax-xymin)/((double)(N0-1)); /* grid spacing */
    const double pbeta=0.5e0*hbar*dtau/(mass*dxy*dxy);     /* defn */
    const double palpha=1.0e0-4.0e0*pbeta; /* defn */
    const double pgamma=-dtau/hbar;        /* defn */
    const double pkinetic=-0.5e0*hbar*hbar/(mass*dxy*dxy); /* defn */

    const int N1_all=N1*nranks;
    const int northTag=0, southTag=1;
    int ncur=0;
    double tcur=0.0e0;
    double psi2integ=0.0e0;
    double* psiold = malloc(sizeof(double) * N0*N1);
    double* psinew = malloc(sizeof(double) * N0*N1);
    double* psisouth = malloc(sizeof(double) * N0);
    double* psinorth = malloc(sizeof(double) * N0);
    double* vpot = malloc(sizeof(double) * N0*N1);
    double psiscale;
    double ekinetic,epotent;
    MPI_Request Rrequests[2], Srequests[2];
    MPI_Status status[2];

    // printf("argv1: %s\n", );
    // printf("argv2: %s\n", );

    int provided;
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &mysize);

    if (mysize != nranks) {
        printf("This program is meant to be run %d ranks.\n", nranks);
        printf("you need to set hardcoded nranks == MPI_Comm_size \n");
        MPI_Finalize();
        return 15;
    }

    if (myrank == 0) {
        printf("Global grid size: %d x %d.\n",N0, N1_all);
        printf("Local slice size: %d x %d.\n",N0, N1);
        printf("resolution dxy = %e.\n",dxy);  
    }

    initialisePsi(psinew, N0, N1, xymin, dxy);
    /* get psi2 integral and then rescale all Psi*/
    psi2integ=getPsi2Integral(psinew, N0, N1, dxy);
    printf("psi2integ before: %e \n",psi2integ);
    psiscale = 1.0e0/sqrt(psi2integ);
    rescalePsi(psinew, N0, N1, psiscale);
    psi2integ=getPsi2Integral(psinew, N0, N1, dxy);
    printf("psi2integ after : %e \n",psi2integ);

    /* calculate the initial kinetic energy */
    ekinetic=getTenergy(psinew, N0, N1, pkinetic, dxy);
    /* initialise/calculate the initial potential energy */
    initialisePot(vpot, N0, N1, xymin, dxy);
    epotent=getVenergy(psinew, N0, N1, vpot, dxy);
/* DEBUG CHECK: output the initial state wavefunction */
/*    printPsi(psinew,N0,N1,xymin,dxy);
    return 0;     */

    if (myrank == 0) {
        printf("imag time = %e Energy = %e + %e = %e \n",tcur,
                                 ekinetic,epotent,(ekinetic+epotent));  
    }

    copyPsi(psiold, psinew, N0, N1);

    while (tcur < tend) {

        int nsends=0;
        /* exchange to the north: */
        if (myrank < mysize-1) {
            MPI_Isend(&psiold[N0*(N1-1)], N0, MPI_DOUBLE, myrank+1, northTag, comm, &Srequests[nsends]);
            nsends++;
            MPI_Irecv(&psinorth[0], N0, MPI_DOUBLE, myrank+1, southTag, comm, &Rrequests[0]);
        }
        
        /* exchange to the south: */
        if (myrank > 0) {
            MPI_Isend(&psiold[0], N0, MPI_DOUBLE, myrank-1, southTag, comm, &Srequests[nsends]);
            nsends++;
            MPI_Irecv(&psisouth[0], N0, MPI_DOUBLE, myrank-1, northTag, comm, &Rrequests[1]);
        }

        /* update psi for the inner grid points */
        updatePsi(psinew, psiold, vpot, N0, N1, palpha, pbeta, pgamma);
        
        
        /* update psi for the grid points in the north
           first we need to wait for the data receive to complete. */
        if (myrank < mysize-1) {
          MPI_Wait(&Rrequests[0], &status[1]); 
          updateNorthPsi(psinew, psiold, psinorth, vpot, N0, N1, palpha, pbeta, pgamma);
        }
        
        /* update psi for the grid points in the south
           first we need to wait for the data receive to complete. */
        if (myrank > 0) {
            MPI_Wait(&Rrequests[1], &status[1]); 
            updateSouthPsi(psinew, psiold, psisouth, vpot, N0, N1, palpha, pbeta, pgamma);
        }
        
        ncur++;
        tcur+=dtau;
        
        /* get psi2 integral and then rescale all Psi*/
        psi2integ=getPsi2Integral(psinew, N0, N1, dxy);
        psiscale = 1.0e0/sqrt(psi2integ);
        rescalePsi(psinew, N0, N1, psiscale);
        psi2integ=getPsi2Integral(psinew, N0, N1, dxy);

        /* calculate the updated energies */
        ekinetic=getTenergy(psinew, N0, N1, pkinetic, dxy);
        printf("ekinetic: N0=%d, N1=%d, pkinetic=%e, dxy=%e\n",N0,N1,pkinetic,dxy);
        epotent=getVenergy(psinew, N0, N1, vpot, dxy);
        printf("epotent: N0=%d, N1=%d, dxy=%e\n",N0,N1,dxy);
        // if (myrank == 0) {
        //   printf("imag time = %e Energy = %e + %e = %e \n",tcur,
        //                          ekinetic,epotent,(ekinetic+epotent));  
        // } 

        if (myrank == 0) {
            printf("iteration number = %d, imag time = %e, psi2integ = %e, ekinetic = %e, epotent = %e, Energy = %e \n", 
                                ncur, tcur, psi2integ, ekinetic, epotent, (ekinetic+epotent));  
        } 

        /* make sure all sends are completed: */
        if (nsends > 0) { 
            MPI_Waitall(nsends, Srequests, status);
        }
        copyPsi(psiold, psinew, N0, N1);

        // printf("iteration number = %d, imag time = %e, nSolution = %e, ekinetuc = %e, epotential = %e, etotal = %e\n",
        //             ncur, tcur, psi2integ, ekinetic, epotent, (ekinetic+epotent));
    }
    if (myrank == 0) {
            printf("Completed after %d time steps.\n",ncur);
            printf("Global grid size: %d x %d.\n",N0, N1_all);
            printf("Number of ranks = %d\n",mysize);
    }    

    MPI_Finalize();

    free(psiold);
    free(psinew);
    free(psisouth);
    free(psinorth);
    free(vpot);
    return 0;
}

