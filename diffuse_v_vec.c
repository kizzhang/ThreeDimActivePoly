#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"

void diffuse_v_vec(double* dir_vec, double* dir_buffer, double* v_vec, double* vbuffer, double* angle, double* abuffer, double dtheta, double dphi);


void diffuse_v_vec(double* dir_vec, double* dir_buffer, double* v_vec, double* vbuffer, double* angle, double* abuffer, double dtheta, double dphi){
    double norm, v_norm;

    double dir_x = *dir_buffer;
    double dir_y = *(dir_buffer+1);
    double dir_z = *(dir_buffer+2);

    double unnorm_x, unnorm_y, unnorm_z;

    double theta = *abuffer;
    double phi = *(abuffer+1);

    double vx = *vbuffer;
    double vy = *(vbuffer+1);
    double vz = *(vbuffer+2);

    unnorm_x = dir_x - sin(angle[0]) * dphi + cos(angle[0]) * cos(angle[1]) * dtheta; //- 2 * rot_diff * time_step * cos(angle[0]) * sin(angle[1]);
    unnorm_y = dir_y + cos(angle[0]) * dphi + sin(angle[0]) * cos(angle[1]) * dtheta; //- 2 * rot_diff * time_step * sin(angle[0]) * cos(angle[1]); 
    unnorm_z = dir_z - sin(angle[1]) * dtheta; //- 2 * rot_diff * time_step * cos(angle[1]);
    norm = sqrt(pow(unnorm_x,2) + pow(unnorm_y,2) + pow(unnorm_z, 2));

    dir_vec[0] = unnorm_x / norm; dir_vec[1] = unnorm_y / norm; dir_vec[2] = unnorm_z / norm;  

    v_norm  = sqrt(pow(dir_vec[0],2) + pow(dir_vec[1],2) + pow(dir_vec[2], 2));
    
    //printf("\n Norm of dir_vec is %lf",  v_norm);

    v_vec[0] = walk_speed * dir_vec[0];
    v_vec[1] = walk_speed * dir_vec[1];
    v_vec[2] = walk_speed * dir_vec[2];

    angle[0] = atan(v_vec[1] / v_vec[0]);
    angle[1] = acos(v_vec[2]);

    /*
    *(v_vec+time_i) = *(v_vec+time_i) + sin(dtheta) / dtheta * (-1. * dthetaz * *(v_vec+time_i+ NB) + dthetay * *(v_vec+time_i+ NB*2))
                + (1-cos(dtheta)) / (dtheta * dtheta) * ( (-1. * dthetaz*dthetaz - dthetay*dthetay)*(*v_vec + time_i) +(dthetax*dthetay) * *(v_vec+time_i+ NB) +
                (dthetax*dthetaz)* *(v_vec+time_i+ NB*2) );


    *(v_vec+time_i+ NB) = *(v_vec+time_i+ NB) + sin(dtheta) / dtheta * (dthetaz * *(v_vec+time_i) - dthetax * *(v_vec+time_i+ NB*2))
                + (1-cos(dtheta)) / (dtheta * dtheta) * ( (dthetax*dthetay)*(*v_vec + time_i) +(-1.*dthetaz*dthetaz - dthetax*dthetax) * *(v_vec+time_i+ NB) +
                (dthetay*dthetaz)* *(v_vec+time_i+ NB*2) );

    *(v_vec+time_i+ NB*2) = *(v_vec+time_i+ NB*2) + sin(dtheta) / dtheta * (-1.*dthetay * *(v_vec+time_i) + dthetax * *(v_vec+time_i+ NB))
                + (1-cos(dtheta)) / (dtheta * dtheta) * ( (dthetax*dthetaz)*(*v_vec + time_i) +(dthetay*dthetaz) * *(v_vec+time_i+ NB) +
                (-1.* dthetax*dthetax -dthetay*dthetay)* *(v_vec+time_i+ NB*2) );*/
}

