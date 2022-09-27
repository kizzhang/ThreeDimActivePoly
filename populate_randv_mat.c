#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"  
#include "myparam.h"

void populate_randv_mat(long randv_mat[NB][2]);
void populate_randv_mat(long randv_mat[NB][2]){
    int i,j,k;
    
    for(i=0; i<NB; i++){
            for(j=0; j<2; j++){
                randv_mat[i][j] = -1*(NB*10+100*i+j);
            }
    }
}