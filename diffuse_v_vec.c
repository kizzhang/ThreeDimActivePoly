#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"

void diffuse_v_vec(double* dir_vec, double* dir_buffer, double* angle, double* abuffer, double dtheta, double dphi);

void diffuse_v_vec(double* dir_vec, double* dir_buffer, double* angle, double* abuffer, double dtheta, double dphi){
    double norm;

    double dir_x = *dir_buffer;
    double dir_y = *(dir_buffer+1);
    double dir_z = *(dir_buffer+2);

    double unnorm_x, unnorm_y, unnorm_z;

    double theta = *abuffer;
    double phi = *(abuffer+1);

    unnorm_x = dir_x - sin(angle[0]) * dphi + cos(angle[0]) * cos(angle[1]) * dtheta; //- 2 * rot_diff * time_step * cos(angle[0]) * sin(angle[1]);
    unnorm_y = dir_y + cos(angle[0]) * dphi + sin(angle[0]) * cos(angle[1]) * dtheta; //- 2 * rot_diff * time_step * sin(angle[0]) * cos(angle[1]); 
    unnorm_z = dir_z - sin(angle[1]) * dtheta; //- 2 * rot_diff * time_step * cos(angle[1]);
    norm = sqrt(pow(unnorm_x,2) + pow(unnorm_y,2) + pow(unnorm_z, 2));

    dir_vec[0] = unnorm_x / norm; dir_vec[1] = unnorm_y / norm; dir_vec[2] = unnorm_z / norm;  
    
    angle[0] = atan(dir_vec[1] / dir_vec[0]);
    angle[1] = acos(dir_vec[2]);
    
}

