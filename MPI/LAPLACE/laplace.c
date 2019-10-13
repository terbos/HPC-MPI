#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main (int argc, char *argv[])
{
   double *sbuf, *rbuf;
   double **phi, **oldphi;
   int    **mask;
   int    i, j, k;
   int    nptsx = 200, nptsy = 200;
   int    nsteps  = 5000;
   int    rank, nprocs, dims[2], periods[2], coords[2];
   int    myposx, myposy, nprocx, nprocy;
   int    up, down, left, right;
   int    bufsize, nsizex, nsizey, nlocalx, nlocaly;
   MPI_Comm new_comm;
   MPI_Status status;

/* Initialise and find rank and number of processes */
   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);

/* Work out number of processes in each direction of the process mesh */
   dims[0] = dims[1] = 0;
   MPI_Dims_create (nprocs, 2, dims);
   nprocy = dims[0];
   nprocx = dims[1];

/* Set up 2D topology */
   periods[0] = periods[1] = 0;
   MPI_Cart_create (MPI_COMM_WORLD, 2, dims, periods, 1, &new_comm);
   MPI_Cart_coords (new_comm, rank, 2, coords);
   myposy = coords[0];
   myposx = coords[1];

/* Determine neighbouring processes for communication shifts */
   MPI_Cart_shift  (new_comm, 0, -1, &up,  &down);
   MPI_Cart_shift  (new_comm, 1, -1, &right, &left);

/* Initialise arrays */
   nsizex = (nptsx-1)/nprocx + 1;
   nsizey = (nptsy-1)/nprocy + 1;
   bufsize = (nsizex>nsizey) ? nsizex : nsizey;
   sbuf = (double *)malloc((sizeof(double)*bufsize));
   rbuf = (double *)malloc((sizeof(double)*bufsize));
   phi = (double **)malloc((sizeof(double*)*(nsizey+2)));
   oldphi = (double **)malloc((sizeof(double*)*(nsizey+2)));
   mask = (int **)malloc((sizeof(int*)*(nsizey+2)));
   for (k=0;k<nsizey+2;k++){
        phi[k] = (double *)malloc(sizeof(double)*(nsizex+2));
        oldphi[k] = (double *)malloc(sizeof(double)*(nsizex+2));
        mask[k] = (int *)malloc(sizeof(int)*(nsizex+2));
   }
   nlocalx = (myposx==nprocx-1) ? nptsx-nsizex*(nprocx-1) : nsizex;
   nlocaly = (myposy==nprocy-1) ? nptsy-nsizey*(nprocy-1) : nsizey;

   printf("rank = %d: (myposx,myposy)=(%d,%d), (left,right,down,up)=(%d,%d,%d,%d), (nlocalx,nlocaly)=(%d,%d)\n",rank,myposx,myposy,left,right,down,up,nlocalx,nlocaly);

   setup_grid (phi, nptsx, nptsy, nprocx, nprocy, myposx, myposy, nlocalx, nlocaly, mask);
 
/* Iterate to find solution */
   for(k=1;k<=nsteps;k++){
         for(j=1;j<=nlocaly;j++)
            for(i=1;i<=nlocalx;i++)
               oldphi[j][i] = phi[j][i];
         MPI_Sendrecv (&oldphi[1][1],  nlocalx, MPI_DOUBLE,  down, 111,
            &oldphi[nlocaly+1][1], nlocalx, MPI_DOUBLE, up, 111,
              new_comm, &status);
         MPI_Sendrecv (&oldphi[nlocaly][1],nlocalx,MPI_DOUBLE, up,112,
            &oldphi[0][1], nlocalx, MPI_DOUBLE,down, 112,
              new_comm, &status);
         for(i=1;i<=nlocaly;i++) sbuf[i-1] = oldphi[i][1];
         MPI_Sendrecv (sbuf, nlocaly, MPI_DOUBLE, left, 113,
                       rbuf, nlocaly, MPI_DOUBLE, right, 113,
                       new_comm, &status);
         for(i=1;i<=nlocaly;i++) oldphi[i][nlocalx+1] = rbuf[i-1];   
         for(i=1;i<=nlocaly;i++) sbuf[i-1] = oldphi[i][nlocalx];
         MPI_Sendrecv (sbuf, nlocaly, MPI_DOUBLE, right, 114,
                       rbuf, nlocaly, MPI_DOUBLE, left, 114,
                       new_comm, &status);
         for(i=1;i<=nlocaly;i++) oldphi[i][0] = rbuf[i-1];
         for(j=1;j<=nlocaly;j++)
            for(i=1;i<=nlocalx;i++)
               if (mask[j][i]) phi[j][i] = 0.25*(oldphi[j][i-1] +
                  oldphi[j][i+1] + oldphi[j-1][i] + oldphi[j+1][i]);
   }
   output_array (phi, rank, nprocx, nprocy, nlocalx, nlocaly, nprocs , myposx, myposy, nptsx, nptsy, new_comm);
 
   MPI_Finalize();
   return 0;
}

setup_grid (phi, nptsx, nptsy, nprocx, nprocy, myposx, myposy, nlocalx, nlocaly, mask)
double  **phi;
int     nptsx, nptsy, nprocx, nprocy, myposx, myposy, nlocalx, nlocaly;
int     **mask;
{
    int i, j, global_x, global_y;

    for(j=0;j<=nlocaly+1;j++)
       for(i=0;i<=nlocalx+1;i++){
          phi[j][i]  = 0.0;
          mask[j][i] = 1;
       }

    if (myposy == 0)
       for(i=0;i<=nlocalx+1;i++) mask[1][i] = 0;

    if (myposy == nprocy-1)
       for(i=0;i<=nlocalx+1;i++) mask[nlocaly][i] = 0;

    if (myposx == 0)
       for(j=0;j<=nlocaly+1;j++) mask[j][1] = 0;

    if (myposx == nprocx-1)
       for(j=0;j<=nlocaly+1;j++) mask[j][nlocalx] = 0;

    for(j=1;j<=nlocaly;j++){
       global_y = nlocaly*myposy + j - 1;
       if (global_y == nptsy/2 || global_y == nptsy/2-1){
          for(i=1;i<=nlocalx;i++){
             global_x = nlocalx*myposx + i - 1;
             if (global_x == nptsx/2 || global_x == nptsx/2-1){
                mask[j][i] = 0;
                phi[j][i]  = 1.0;
             }
          }
       }
    }
}

output_array (phi, rank, nprocx, nprocy, nlocalx, nlocaly, nprocs, myposx, myposy, nptsx, nptsy, new_comm)
double **phi;
int    rank, nprocx, nprocy, nlocalx, nlocaly, nprocs, myposx, myposy, nptsx, nptsy;
MPI_Comm new_comm;
{
   int i, j, k=0, m, n;
   int jmax, nsizey, count, source;
   int coords[2];
   MPI_Status status;
   FILE *fp;

   
   if(rank==0){
	 fp = fopen("out.ps","w");
         fprintf(fp,"/picstr %d string def\n",nptsx);
         fprintf(fp,"50 50 translate\n");
         fprintf(fp,"%d %d scale\n",nptsx, nptsy);
         fprintf(fp,"%d %d 8 [%d 0 0 %d 0 %d] \n",nptsx, nptsy, nptsx, nptsy, -nptsx);
         fprintf(fp,"{currentfile 3 200 mul string readhexstring pop} bind false 3 colorimage\n");
   }

   nsizey = (nptsy-1)/nprocy + 1;
   for(m=0;m<nprocy;m++){
      jmax = (m==nprocy-1) ? nptsy - nsizey*(nprocy-1) : nsizey;
      for(j=1;j<=jmax;j++){
         for(n=0;n<nprocx;n++){
            if (myposx == n && myposy == m)
               if (rank!=0){
                      MPI_Send (&phi[j][1], nlocalx, MPI_DOUBLE, 0, 115,MPI_COMM_WORLD);
               }
               else
                  for(i=1;i<=nlocalx;i++,k++){
                     fprintf(fp,"%06x",RGBval(phi[j][i]));
                     if((k+1)%10==0) fprintf(fp,"\n");
//                     printf("%9.6f%c",phi[j][i],((k+1)%10==0) ? '\n' : ' ');
                  }
            else if (rank==0){
               coords[0] = m;
               coords[1] = n;
               MPI_Cart_rank (new_comm, coords, &source);
               MPI_Recv (&phi[0][1], nlocalx, MPI_DOUBLE, source, 115, MPI_COMM_WORLD, &status);
               MPI_Get_count(&status, MPI_DOUBLE, &count);
// need to check number received here
               for(i=1;i<=count;i++,k++){
                     fprintf(fp,"%06x",RGBval(phi[0][i]));
                     if((k+1)%10==0) fprintf(fp,"\n");
//                  printf("%9.6f%c",phi[1][i],((k+1)%10==0) ? '\n' : ' ');
               }
            }
         }
     }
   }
    if(rank==0) fclose(fp);
}

int RGBval(double x){
    int R, B, G, pow8 = 256;
    if(x<=0.5){
        B = (int)((1.0-2.0*x)*255.0);
        G = (int)(2.0*x*255.0);
	R = 0; 
    }
    else{
        B = 0;
        G = (int)((2.0-2.0*x)*255.0);
        R = (int)((2.0*x-1.0)*255.0);
    }
    return (B+(G+R*pow8)*pow8);
}
