#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"

void disp_diffuse(double next_pt[NB][3], double x[NB][3], double* dx, double* dy, double* dz){
    int i;
    
    for(i=0; i<NB; i++){
        //if(i ==NB-2){printf("\n\ndx[%d] = %lf; dy[%d] = %lf; dz[%d] = %lf\n\n",i, dx[i], i, dy[i], i, dz[i]);}
        next_pt[i][0] = x[i][0] + dx[i]; next_pt[i][1] = x[i][1] + dy[i]; next_pt[i][2] = x[i][2] + dz[i];
    }

}