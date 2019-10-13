#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
   double **phi, **oldphi;
   int    **mask;
   int    i, j, k;
   int    nptsx = 200, nptsy = 200;
   int    nsteps  = 500;

   phi = (double **)malloc(sizeof(double*)*nptsy);
   oldphi = (double **)malloc(sizeof(double*)*nptsy);
   mask = (int **)malloc(sizeof(int*)*nptsy);
   for (k=0;k<nptsy;k++){
        phi[k] = (double *)malloc(sizeof(double)*nptsx);
        oldphi[k] = (double *)malloc(sizeof(double)*nptsx);
        mask[k] = (int *)malloc(sizeof(int)*nptsx);
   }

   setup_grid (phi, nptsx, nptsy, mask);
 
/* Iterate to find solution */
   for(k=1;k<=nsteps;k++){
         for(j=0;j<nptsy;j++)
            for(i=0;i<nptsx;i++)
               oldphi[j][i] = phi[j][i];
         for(j=0;j<nptsy;j++)
            for(i=0;i<nptsx;i++)
               if (mask[j][i]) phi[j][i] = 0.25*(oldphi[j][i-1] +
                  oldphi[j][i+1] + oldphi[j-1][i] + oldphi[j+1][i]);
   }
   output_array (phi, nptsx, nptsy);
 
   return 0;
}

setup_grid (phi, nptsx, nptsy, mask)
double  **phi;
int     nptsx, nptsy;
int     **mask;
{
    int i, j, nx2, ny2;

    for(j=0;j<nptsy;j++)
       for(i=0;i<nptsx;i++){
          phi[j][i]  = 0.0;
          mask[j][i] = 1;
       }

    for(i=0;i<nptsx;i++) mask[0][i] = 0;

    for(i=0;i<nptsx;i++) mask[nptsy-1][i] = 0;

    for(j=0;j<nptsy;j++) mask[j][0] = 0;

    for(j=0;j<nptsy;j++) mask[j][nptsx-1] = 0;

    nx2 = nptsx/2;
    ny2 = nptsy/2;
    mask[ny2][nx2] = 0;
    mask[ny2][nx2-1] = 0;
    mask[ny2-1][nx2] = 0;
    mask[ny2-1][nx2-1] = 0;
    phi[ny2][nx2]  = 1.0;
    phi[ny2][nx2-1]  = 1.0;
    phi[ny2-1][nx2]  = 1.0;
    phi[ny2-1][nx2-1]  = 1.0;
}

output_array (phi, nptsx, nptsy)
double **phi;
int    nptsx, nptsy;
{
   int i, j, k=0;
   FILE *fp;

   
   fp = fopen("outSEQ.ps","w");
   fprintf(fp,"/picstr %d string def\n",nptsx);
   fprintf(fp,"50 50 translate\n");
   fprintf(fp,"%d %d scale\n",nptsx, nptsy);
   fprintf(fp,"%d %d 8 [%d 0 0 %d 0 %d] \n",nptsx, nptsy, nptsx, nptsy, -nptsx);
   fprintf(fp,"{currentfile 3 200 mul string readhexstring pop} bind false 3 colorimage\n");

   for(j=0;j<nptsy;j++){
        for(i=0;i<nptsx;i++,k++){
             fprintf(fp,"%06x",RGBval(phi[j][i]));
             if((k+1)%10==0) fprintf(fp,"\n");
        }
   }
   fclose(fp);
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
