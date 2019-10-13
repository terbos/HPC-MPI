#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define min(A,B) ((A)<(B) ? (A) : (B))
#define PI 3.1415926536
#define NSTEPS 10000
#define OUTFREQ 100
#define FILENAME_LEN 32

int main (int argc, char *argv[])
{
   double *psi, *new_psi, *old_psi, *buf;
   int rank, nprocs, mypos;
   int periods=0, reorder=1;
   int icheck, i, j, nsize, nlocal, nbeg, nend, istart, iend;
   int count_out=0, tag=111;
   int npts, left, right;
   int output_solution();
   double x, tau=0.05;
   MPI_Comm new_comm;
   MPI_Status status;

/* Initialise and find rank and number of processes */
   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);

/* Read in input parameters */
   if (rank==0) {
      	printf("\nGive the number of points => \n");
      	icheck = scanf ("%d", &npts);
   }
   
   MPI_Bcast (&npts, 1, MPI_INT, 0, MPI_COMM_WORLD);

/* Set up 1D topology */
   MPI_Cart_create (MPI_COMM_WORLD,1,&nprocs,&periods,reorder,&new_comm);
   MPI_Cart_coords (new_comm, rank, 1, &mypos);
   MPI_Cart_shift  (new_comm, 0, -1, &right, &left);
   
   nsize  = (npts-1)/nprocs + 1;
   psi = (double *)malloc(sizeof(double)*(nsize+2));
   new_psi = (double *)malloc(sizeof(double)*(nsize+2));
   old_psi = (double *)malloc(sizeof(double)*(nsize+2));
   buf = (double *)malloc(sizeof(double)*(nsize+2));

   nbeg   = mypos*nsize;
   nend   = min (nbeg+nsize-1, npts-1);
   nlocal = nend - nbeg + 1;
   printf("rank = %d, mypos = %d, left = %d, right = %d, nbeg = %d, nend = %d, nlocal = %d\n",rank,mypos,left,right,nbeg, nend, nlocal);

/* Initialise array */
   for(i=nbeg;i<=nend;i++){
          x = 2.0*PI*(double)(i)/(double)(npts-1);
          x = sin(x);
          psi[i+1-nbeg] = old_psi[i+1-nbeg] = x;
   }

/* Do the updates */
   istart = (mypos==0) ? 2 : 1;
   iend   = (mypos==nprocs-1) ? nlocal-1 : nlocal;
   for(j=0;j<NSTEPS;j++){
         MPI_Sendrecv (&psi[1],        1, MPI_DOUBLE, left,  tag, 
                       &psi[nlocal+1], 1, MPI_DOUBLE, right, tag,
                       new_comm, &status);
         MPI_Sendrecv (&psi[nlocal], 1, MPI_DOUBLE, right,  tag, 
                       &psi[0],      1, MPI_DOUBLE, left,   tag,
                       new_comm, &status);
         for(i=istart;i<=iend;i++){
            new_psi[i] = 2.0*psi[i] - old_psi[i] + tau*tau*(psi[i-1]-2.0*psi[i]+psi[i+1]);
         }
         for(i=1;i<=nlocal;i++){
            old_psi[i] = psi[i];
            psi[i]     = new_psi[i];
         }
         if(j%OUTFREQ == 0) {
                count_out++;
		output_solution(psi, buf, new_comm, rank, mypos, nprocs, nsize, nlocal, count_out);
         }
   }
/*
   if(rank==0){
         for(j=0;j<nprocs;j++){
             if (mypos!=j){
                 MPI_Cart_rank (new_comm, &j, &proc);
                 MPI_Recv (&old_psi[1], nsize, MPI_DOUBLE, proc, tag, new_comm, &status);
                 MPI_Get_count(&status, MPI_DOUBLE, &count);
                 for(i=1;i<=count;i++)
                     printf("%9.6f%c",old_psi[i],(i%8==0) ? '\n' : ' ');
             }
             else{
                 for(i=1;i<=nlocal;i++)
                     printf("%9.6f%c",psi[i],(i%8==0) ? '\n' : ' ');
             }
         }
    }
    else{
        MPI_Send(&psi[1], nlocal, MPI_DOUBLE, 0, tag, new_comm);
    }
*/


   MPI_Finalize ();
   return 0;
}

int output_solution(double *psi, double *buf, MPI_Comm new_comm, int rank, int mypos, int nprocs, int nsize, int nlocal, int count_out)
{
    int i, j, proc, count, tag=99;
    MPI_Status status;
    FILE *fp;
    const char base[] = "OUT/";
    char filename [ FILENAME_LEN ];

    if(rank==0){
         sprintf(filename, "%s%04d", base, count_out);
         fp = fopen(filename,"w");
         for(j=0;j<nprocs;j++){
             if (mypos!=j){
                 MPI_Cart_rank (new_comm, &j, &proc);
                 MPI_Recv (&buf[1], nsize, MPI_DOUBLE, proc, tag, new_comm, &status);
                 MPI_Get_count(&status, MPI_DOUBLE, &count);
                 for(i=1;i<=count;i++)
//                     printf("%9.6f%c",old_psi[i],(i%8==0) ? '\n' : ' ');
                     fprintf(fp,"[%9.6f %9.6f]\n",(double)(j*nsize+i),buf[i]);
             }
             else{
                 for(i=1;i<=nlocal;i++)
//                     printf("%9.6f%c",psi[i],(i%8==0) ? '\n' : ' ');
                     fprintf(fp,"[%9.6f %9.6f]\n",(double)(j*nsize+i),psi[i]);
             }
         }
         fclose(fp);
    }
    else{
        MPI_Send(&psi[1], nlocal, MPI_DOUBLE, 0, tag, new_comm);
    }
    return 0;
}
