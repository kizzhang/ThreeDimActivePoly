#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"  
#include "myparam.h"

void populate_randv_mat(double randv_mat[TS][NB][2]);
void populate_randv_mat(double randv_mat[TS][NB][2]){
    long dum[NB][2];
    int i,j,k;
    
    for(i=0; i<NB; i++){
            for(j=0; j<2; j++){
                dum[i][j] = -1*(100+100*i+j);
            }
    }

    for(i=0; i<TS; i++){
        for(j=0; j<NB; j++){
            for(k=0;k<2;k++){
                randv_mat[i][j][0] = ran2(&dum[j][k], 'g');
                randv_mat[i][j][1] = ran2(&dum[j][k], 'g');
                //printf("\n rand_mat[%d][%d][%d] = %lf\n\n",i,j,0, rand_mat[i][j][0])
            }
        }
    }    
}