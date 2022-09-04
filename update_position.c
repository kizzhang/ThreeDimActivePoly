#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"
#include "distance_matrix.h"

void distance_matrix(double *X, double* distMat, const int row);



void update_position(double* dx, double* dy, double* dz, double** x, double** v_vec, long* dvar){
    int num; 
    int i,j;
    double force_field[NB][dd];
    double dist_matrix[NB][NB];
    double sigma = 1;
    double epsilon = 1; //k_b * T 
    double K_fene = 27*epsilon/pow(sigma,2);
    double R_0 = 1.5*sigma;

    distance_matrix(x, dist_matrix, NB);

    for(i=0; i<NB; i++){
        for(j=i+1; j<NB; j++){
            force_field[i][0] = force_field[i][0] + (x[i][0]-x[j][0]) / dist_matrix[i][j] * (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)) +\
                                K_fene * dist_matrix[i][j] / ( 1 - pow(dist_matrix[i][j]/R_0, 2)));
            force_field[i][1] = force_field[i][1] + (x[i][1]-x[j][1]) / dist_matrix[i][j] * (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)) +\
                                K_fene * dist_matrix[i][j] / ( 1 - pow(dist_matrix[i][j]/R_0, 2)));
            force_field[i][2] = force_field[i][2] + (x[i][2]-x[j][2]) / dist_matrix[i][j] * (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)) +\
                                K_fene * dist_matrix[i][j] / ( 1 - pow(dist_matrix[i][j]/R_0, 2)));
        }
    }

    dx[0] = k_spring/gamma * (x[1][0] - 2 * x[0][0]) * time_step + force_field[0][0] * time_step\
            + a_force * v_vec[0][0] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dy[0] = k_spring/gamma * (x[1][1] - 2 * x[0][1]) * time_step + force_field[0][1] * time_step\
            + a_force * v_vec[0][1] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dz[0] = k_spring/gamma * (x[1][2] - 2 * x[0][2]) * time_step + force_field[0][2] * time_step\
            + a_force * v_vec[0][2] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');

    for(num=1; num<NB-1; num++){
        dx[num] = k_spring/gamma * (x[num+1][0] + x[num-1][0] - 2 * x[num][0]) * time_step  \
                + a_force * v_vec[num][0] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g')\
                + force_field[num][0] * time_step;
        dy[num] = k_spring/gamma * (x[num+1][1] + x[num-1][1] - 2 * x[num][1]) * time_step \
               + a_force * v_vec[num][1] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g')\
               + force_field[num][1] * time_step;
        dz[num] = k_spring/gamma * (x[num+1][2] + x[num-1][2] - 2 * x[num][2]) * time_step \
               + a_force * v_vec[num][2] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g')\
               + force_field[num][2] * time_step;
    }

    dx[NB-1] = k_spring/gamma * (x[NB-2][0] - 2 * x[NB-1][0]) * time_step + force_field[NB-1][0] * time_step\
                + a_force * v_vec[NB-1][0] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dy[NB-1] = k_spring/gamma * (x[NB-2][1] - 2 * x[NB-1][1]) * time_step + force_field[NB-1][1] * time_step\
                + a_force * v_vec[NB-1][1] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dz[NB-1] = k_spring/gamma * (x[NB-2][2] - 2 * x[NB-1][2]) * time_step + force_field[NB-1][2] * time_step\
                + a_force * v_vec[NB-1][2] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');

}