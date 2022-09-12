#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"
#include "distance_matrix.h"

void distance_matrix(double *X, double* distMat, const int row);

void update_position(double* dx, double* dy, double* dz, double x[][3], double dir_vec[][3]){
    int num; 
    int i,j;
    double force_field[NB][dd];
    double dist_matrix[NB][NB];
    double sigma = 0.2;
    double k_spring = 4*rot_diff/(3*pos_diff); // CHANGE LATER
    double gamma = 0.0038; // CHANGE LATER
    double a_force = 0.0038; // CHANGE LATER
    double epsilon = 1; //k_b * T 
    double K_fene = 27*epsilon/pow(sigma,2);
    double R_0 = 1.5*sigma;
    long dum = -100;
    long *dvar;

    dvar=&dum;
    
    /*
    printf("x[0][0] = %lf; x[0][1] = %lf; x[0][2] = %lf", x[0][0], x[0][1], x[0][2]);
    printf("\n");
    printf("x[6][0] = %lf; x[6][1] = %lf; x[6][2] = %lf", x[6][0], x[6][1], x[6][2]);
    printf("\n");*/

    distance_matrix(*x, *dist_matrix, NB);
 /*
    for(i=0; i<NB; i++){
        for(j=0; j<NB; j++){
            printf("%lf ", dist_matrix[i][j]);
            if (j == NB-1) {
                printf("\n");
            }
        }
    } 
    printf("\n\n");*/
    //printf("K_fene = %lf", K_fene);
    for(i=0; i<NB; i++){
        for(j=0; j<3; j++){
           force_field[i][j] = 0;
        }
    }

    for(i=0; i<NB; i++){
        for(j=i; j<NB; j++){
            if(i==j){continue;}

            if(dist_matrix[i][j] <= pow(2,1/6) * sigma){
                force_field[i][0] = force_field[i][0] + (x[i][0]-x[j][0]) / dist_matrix[i][j] *\
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)));
                force_field[i][1] = force_field[i][1] + (x[i][1]-x[j][1]) / dist_matrix[i][j] *\
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)));
                force_field[i][2] = force_field[i][2] + (x[i][2]-x[j][2]) / dist_matrix[i][j] *\
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)));
                             
            }   
            force_field[i][0] = force_field[i][0] + (x[i][0]-x[j][0]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) ));
            force_field[i][1] = force_field[i][1] + (x[i][1]-x[j][1]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) ));
            force_field[i][2] = force_field[i][2] + (x[i][2]-x[j][2]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) ));
           /* if(i==5){
                printf("\n\n(x[%d][2]-x[%d][2]) = %lf \n\n",i,j,(x[i][2]-x[j][2]));
                //printf("\n\ndist_mat[5][%d] = %lf \n\n",j,dist_matrix[5][j]);
                printf("\n\nforce_field[%d][0] = %lf\n",i, force_field[i][0]);
                printf("force_field[%d][1] = %lf\n",i, force_field[i][1]);    
                printf("force_field[%d][2] = %lf\n\n",i, force_field[i][2]);} */
            /*printf("\n\nforce_field[%d][0] = %lf\n",i, force_field[i][0]);
            printf("force_field[%d][1] = %lf\n",i, force_field[i][1]);    
            printf("force_field[%d][2] = %lf\n\n",i, force_field[i][2]);   */
        }
    }

    //printf("dx[6] = %lf; dy[6] = %lf; dz[6] = %lf\n\n", dx[6], dy[6], dz[6]);
    //printf("dir_vec[6] = %lf; dir_vec[6] = %lf; dir_vec[6] = %lf\n\n", dir_vec[6], dir_vec[6], dir_vec[6]);
    dx[0] = -1 * k_spring/gamma * (x[1][0] - 2 * x[0][0]) * time_step + force_field[0][0] * time_step\
            + a_force * dir_vec[0][0] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dy[0] = -1 * k_spring/gamma * (x[1][1] - 2 * x[0][1]) * time_step + force_field[0][1] * time_step\
            + a_force * dir_vec[0][1] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dz[0] = -1 * k_spring/gamma * (x[1][2] - 2 * x[0][2]) * time_step + force_field[0][2] * time_step\
            + a_force * dir_vec[0][2] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    
    for(num=1; num<NB-1; num++){
        dx[num] = -1 * k_spring/gamma * (x[num+1][0] + x[num-1][0] - 2 * x[num][0]) * time_step  \
                + a_force * dir_vec[num][0] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g')\
                + force_field[num][0] * time_step;
        dy[num] = -1 * k_spring/gamma * (x[num+1][1] + x[num-1][1] - 2 * x[num][1]) * time_step \
                + a_force * dir_vec[num][1] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g')\
                + force_field[num][1] * time_step;
        dz[num] = -1 * k_spring/gamma * (x[num+1][2] + x[num-1][2] - 2 * x[num][2]) * time_step \
                + a_force * dir_vec[num][2] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g')\
                + force_field[num][2] * time_step;
    }

    //printf("dx[6] = %lf; dy[6] = %lf; dz[6] = %lf\n\n", dx[6], dy[6], dz[6]);

    dx[NB-1] = -1 * k_spring/gamma * (x[NB-2][0] - 2 * x[NB-1][0]) * time_step + force_field[NB-1][0] * time_step\
                + a_force * dir_vec[NB-1][0] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dy[NB-1] = -1 * k_spring/gamma * (x[NB-2][1] - 2 * x[NB-1][1]) * time_step + force_field[NB-1][1] * time_step\
                + a_force * dir_vec[NB-1][1] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');
    dz[NB-1] = -1 * k_spring/gamma * (x[NB-2][2] - 2 * x[NB-1][2]) * time_step + force_field[NB-1][2] * time_step\
                + a_force * dir_vec[NB-1][2] * time_step + sqrt(2*pos_diff*time_step) * ran2(dvar, 'g');

}