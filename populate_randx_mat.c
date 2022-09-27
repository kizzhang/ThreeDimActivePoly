  #include <math.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include "ran2.h"  
  #include "myparam.h"

void populate_randx_mat(long randx_mat[NB][dd]);
void populate_randx_mat(long randx_mat[NB][dd]){
    int i,j;
    
    for(i=0; i<NB; i++){
            for(j=0; j<dd; j++){
                randx_mat[i][j] = -1* (NB*20+100*i+j);
            }
    }
    
}