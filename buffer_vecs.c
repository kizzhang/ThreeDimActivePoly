#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"


void buffer_vecs(double** dir_vec, double** v_vec, double** angle, double** dir_buffer, double** v_buffer, double** a_buffer){
    int i;

    for(i=0; i<NB; i++){
        dir_buffer[i][0] = dir_vec[i][0]; dir_buffer[i][1] = dir_vec[i][1]; dir_buffer[i][2] = dir_vec[i][2];
        v_buffer[i][0] = v_vec[i][0]; v_buffer[i][1] = v_vec[i][1]; v_buffer[i][2] = v_vec[i][2];
        a_buffer[i][0] = angle[i][0]; v_buffer[i][1] = angle[i][1];
    }
}