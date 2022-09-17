#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"
#include "distance_matrix.h"

void distance_matrix(double *X, double* distMat, const int row);

void update_position(double* dx, double* dy, double* dz, double x[NB][3], double dir_vec[NB][3], double randx_mat[NB][dd]){
    printf("\nYES!!\n");
    int num; 
    int i,j;
    double force_field[NB][dd];
    double dist_matrix[NB][NB];
    double sigma = 0.8;
    double k_spring = 4*rot_diff/(3*pos_diff); 
    double gamma = 0.038; // CHANGE LATER
    double a_force = 0.038; // CHANGE LATER
    double epsilon = 1e-7; //k_b * T 
    double K_fene =  27*epsilon/pow(sigma,2);
    double R_0 = 1.5*sigma;

    distance_matrix(*x, *dist_matrix, NB);
    printf("\nYES!!\n");
    /*
    printf("\n\n");
    for(i=0; i<NB; i++){
        for(j=0; j<NB; j++){
            printf("%lf ", dist_matrix[i][j]);
            if (j == NB-1) {
                printf("\n");
            }
        }
    } 
    printf("\n\n");
    */

    //printf("K_fene = %lf", K_fene);
    for(i=0; i<NB; i++){
        for(j=0; j<3; j++){
           force_field[i][j] = 0;
        }
    }

    for(i=0; i<NB; i++){
        for(j=i; j<NB; j++){
            if(i==j){continue;}
            /*
            if(dist_matrix[i][j] <= pow(2,1/6) * sigma){
                force_field[i][0] = force_field[i][0] + epsilon * (x[i][0]-x[j][0]) / dist_matrix[i][j] *\
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)));
                force_field[i][1] = force_field[i][1] + epsilon * (x[i][1]-x[j][1]) / dist_matrix[i][j] *\
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)));
                force_field[i][2] = force_field[i][2] + epsilon * (x[i][2]-x[j][2]) / dist_matrix[i][j] *\
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7)));
                printf("\n\nforce_field[%d][0] = %lf\n\n",i, force_field[i][0]);
                printf("\n\nCOLLIDED!!\n\n");
            }   
            force_field[i][0] = force_field[i][0] + (x[i][0]-x[j][0]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) ));
            force_field[i][1] = force_field[i][1] + (x[i][1]-x[j][1]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) ));
            force_field[i][2] = force_field[i][2] + (x[i][2]-x[j][2]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) ));
            */
            //printf("\n\nforce_field[%d][0] = %lf\n\n",i,j, force_field[i][0]);
            /*if(dist_matrix[i][j] <= pow(2,1/6)*sigma && abs(force_field[i][0])>1){
                printf("\n\n(x[%d][0]-x[%d][0]) = %lf \n\n",i,j,(x[i][0]-x[j][0]));
                printf("\n\nRepulsion x = %lf \n",epsilon * (x[i][0]-x[j][0]) / dist_matrix[i][j] * \
                                        (48 * (pow(sigma,12) / pow(dist_matrix[i][j],13) - 0.5 * pow(sigma,6) / pow(dist_matrix[i][j],7))));
                printf("\n FENE force x = %lf \n\n",(x[i][0]-x[j][0]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) )));
                printf("\n\ndist_mat[%d][%d] = %lf \n\n",i,j,dist_matrix[i][j]);
                printf("\n\nforce_field[%d][0] = %lf\n",i, force_field[i][0]);
                printf("force_field[%d][1] = %lf\n",i, force_field[i][1]);    
                printf("force_field[%d][2] = %lf\n\n",i, force_field[i][2]);} */
            /*printf("\n\nforce_field[%d][0] = %lf\n",i, force_field[i][0]);
            printf("force_field[%d][1] = %lf\n",i, force_fi1eld[i][1]);    
            printf("force_field[%d][2] = %lf\n\n",i, force_field[i][2]);   */
        }
    }

    dx[0] = 1/gamma * ( -1 * k_spring * (x[1][0] - 2 * x[0][0]) * time_step + force_field[0][0] * time_step\
            + a_force * dir_vec[0][0] * time_step) + sqrt(2*pos_diff*time_step) * randx_mat[0][0];
    dy[0] = 1/gamma * ( -1 * k_spring * (x[1][1] - 2 * x[0][1]) * time_step + force_field[0][1] * time_step\
            + a_force * dir_vec[0][1] * time_step) + sqrt(2*pos_diff*time_step) * randx_mat[0][1];
    dz[0] = 1/gamma * ( -1 * k_spring * (x[1][2] - 2 * x[0][2]) * time_step + force_field[0][2] * time_step\
            + a_force * dir_vec[0][2] * time_step) + sqrt(2*pos_diff*time_step) * randx_mat[0][2];

    for(num=1; num<NB-1; num++){
        dx[num] = 1/gamma * ( -1 * k_spring * (x[num+1][0] + x[num-1][0] - 2 * x[num][0]) * time_step  \
                + a_force * dir_vec[num][0] * time_step + force_field[num][0] * time_step)\
                + sqrt(2*pos_diff*time_step) * randx_mat[num][0];
        dy[num] = 1/gamma * ( -1 * k_spring * (x[num+1][1] + x[num-1][1] - 2 * x[num][1]) * time_step \
                + a_force * dir_vec[num][1] * time_step + force_field[num][1] * time_step)\
                + sqrt(2*pos_diff*time_step) * randx_mat[num][1];
        dz[num] = 1/gamma * ( -1 * k_spring* (x[num+1][2] + x[num-1][2] - 2 * x[num][2]) * time_step \
                + a_force * dir_vec[num][2] * time_step + force_field[num][2] * time_step)\
                + sqrt(2*pos_diff*time_step) * randx_mat[num][2];
    }

    dx[NB-1] = 1/gamma * ( -1 * k_spring * (x[NB-2][0] - 2 * x[NB-1][0]) * time_step + force_field[NB-1][0] * time_step\
                + a_force * dir_vec[NB-1][0] * time_step) + sqrt(2*pos_diff*time_step) * randx_mat[NB-1][0];
    dy[NB-1] = 1/gamma * ( -1 * k_spring * (x[NB-2][1] - 2 * x[NB-1][1]) * time_step + force_field[NB-1][1] * time_step\
                + a_force * dir_vec[NB-1][1] * time_step) + sqrt(2*pos_diff*time_step) * randx_mat[NB-1][1];
    dz[NB-1] = 1/gamma * ( -1 * k_spring * (x[NB-2][2] - 2 * x[NB-1][2]) * time_step + force_field[NB-1][2] * time_step\
                + a_force * dir_vec[NB-1][2] * time_step) + sqrt(2*pos_diff*time_step) * randx_mat[NB-1][2];
   // printf("\n\ndx[%d] = %lf; dy[%d] = %lf; dz[%d] = %lf\n\n", NB-1, dx[NB-1], NB-1, dy[NB-1], NB-1, dz[NB-1]);
    //printf("\n\nx[%d][0] - 2 * x[%d][0]")
}