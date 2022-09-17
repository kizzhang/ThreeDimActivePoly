  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include "ran2.h"  
  #include "myparam.h"

void populate_randx_mat(double randx_mat[TS][NB][dd]);
void populate_randx_mat(double randx_mat[TS][NB][dd]){
    long dum[NB][dd];
    int i,j,k;
    
    for(i=0; i<NB; i++){
            for(j=0; j<dd; j++){
                dum[i][j] = (100+100*i+j);
            }
    }

    for(i=0; i<TS; i++){
        for(j=0; j<NB; j++){
            for(k=0;k<3;k++){
                randx_mat[i][j][0] = ran2(&dum[j][k], 'g');
                randx_mat[i][j][1] = ran2(&dum[j][k], 'g');
                randx_mat[i][j][2] = ran2(&dum[j][k], 'g');
                //printf("\n rand_mat[%d][%d][%d] = %lf\n\n",i,j,0, randx_mat[i][j][0]);
            }
        }
    }    
    
}