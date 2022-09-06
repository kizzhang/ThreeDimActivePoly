#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"
#include "diffuse_v_vec.h"

void update_v_vec(double dir_vec[][3], double dir_buffer[][3], double v_vec[][3], double vbuffer[][3], double angle[][2], double abuffer[][2]);
void diffuse_v_vec(double* dir_vec, double* dir_buffer, double* v_vec, double* vbuffer, double* angle, double* abuffer, double dtheta, double dphi);

void update_v_vec(double dir_vec[][3], double dir_buffer[][3], double v_vec[][3], double vbuffer[][3], double angle[][2], double abuffer[][2]){
    int i; 
    long dum = -100;
    long *dvar;
    double dtheta[NB];
    double dphi[NB];
    
    dvar=&dum;

    for(i=0; i<NB; i++){
        dtheta[i] = ran2(dvar, 'g') * sqrt(2*rot_diff*time_step);
        dphi[i] = ran2(dvar, 'g') * sqrt(2*rot_diff*time_step);
        diffuse_v_vec(dir_vec[i], dir_buffer[i], v_vec[i], vbuffer[i], angle[i], abuffer[i], dtheta[i], dphi[i]);
    }
}