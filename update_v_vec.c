#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"
#include "diffuse_v_vec.h"

void update_v_vec(double dir_vec[NB][3], double dir_buffer[NB][3], double angle[NB][2], double abuffer[NB][2], double randv_mat[NB][2]);
void diffuse_v_vec(double* dir_vec, double* dir_buffer, double* angle, double* abuffer, double dtheta, double dphi);

void update_v_vec(double dir_vec[NB][3], double dir_buffer[NB][3], double angle[NB][2], double abuffer[NB][2], double randv_mat[NB][2]){
    int i; 
    double dtheta[NB];
    double dphi[NB];
    
    printf("\nYES!!\n");
    for(i=0; i<NB; i++){
        dtheta[i] = randv_mat[i][0]* sqrt(2*rot_diff*time_step);
        dphi[i] = randv_mat[i][1] * sqrt(2*rot_diff*time_step);
        /*printf("v_vec[%d][0] is %lf\n", i, v_vec[i][0]);
        printf("v_vec[%d][0] is %lf\n", i, v_vec[i][1]);
        printf("v_vec[%d][2] is %lf\n", i, v_vec[i][2]);
        printf("angle[%d][0] is %lf\n", i, angle[i][0]);
        printf("angle[%d][1] is %lf\n", i, angle[i][1]);*/

        diffuse_v_vec(dir_vec[i], dir_buffer[i], angle[i], abuffer[i], dtheta[i], dphi[i]);
        
        
        //printf("dtheta[%d] is %lf\n", i, dtheta[i]);
        //printf("dphi[%d] is %lf\n", i, dphi[i]);
    }
    
}