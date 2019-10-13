#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 128

main()
{
    double psi[MAX_POINTS];
    double new_psi[MAX_POINTS];
    double old_psi[MAX_POINTS];
    double pi = 3.141592654;
    double tau = 0.05;
    double x;
    int    i, j;
    int    npoints = 128;
    int    nsteps  = 100;

    for(i=0;i<npoints;i++){
        x = 2.0*pi*(double)i/(double)(npoints-1);
        x = sin(x);
        psi[i] = old_psi[i] = x;
    }

    for(j=0;j<nsteps;j++){
        for(i=1;i<npoints-1;i++){
            new_psi[i] = 2.0*psi[i] - old_psi[i] +
                tau*tau*(psi[i-1]-2.0*psi[i]+psi[i+1]);
        }
        for(i=0;i<npoints;i++){
            old_psi[i] = psi[i];
            psi[i]     = new_psi[i];
        }
    }

    for(i=0;i<npoints;i++)
        printf("%9.6f%c",psi[i],((i+1)%8==0) ? '\n' : ' ');
 
    exit(0);
}
