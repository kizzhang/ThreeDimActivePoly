#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"
#include "distance_matrix.h"

void distance_matrix(double *X, double* distMat, const int row);

double force_field[NB][dd];
double dist_matrix[NB][NB];
double x_mean[NB][dd];
void update_position(double* dx, double* dy, double* dz, double x[NB][3], double dir_vec[NB][3], long randx_mat[NB][dd]){
    int num; 
    int i,j;
    double k_spring = 3*KbT/pow(b_khun,2); 
    double a_force = 0; // 0.038; // CHANGE LATER
    double epsilon = KbT; //k_b * T 
    double K_fene =  27*epsilon/pow(sigma,2);
    double R_0 = 1.5*sigma;
    double rlogarg;
    double LJForce;

    distance_matrix(*x, *dist_matrix, NB);
    
    for(i=0; i<NB; i++){
        for(j=0; j<3; j++){
           force_field[i][j] = 0;
        }
    }

    for(i=0; i<NB; i++){
        
        //printf("\n\nforce_field[%d][0] = %lf\n\n",i, force_field[i][0]);
        for(j=i; j<NB; j++){
            if(i==j){continue;}
            
            if(dist_matrix[i][j] <= pow(2.,1./6.) * sigma){
                LJForce = epsilon / dist_matrix[i][j] * (48. * (pow(sigma,12.) / pow(dist_matrix[i][j],13.) - 0.5 * pow(sigma,6.) / pow(dist_matrix[i][j],7.)));

                force_field[i][0] = force_field[i][0] + (x[i][0]-x[j][0]) * LJForce;
                force_field[i][1] = force_field[i][1] + (x[i][1]-x[j][1]) * LJForce;
                force_field[i][2] = force_field[i][2] + (x[i][2]-x[j][2]) * LJForce;
                //printf("\n\nj = %d; LJ force[%d][0] = %lf\n\n",j,i,   force_field[i][0]);
                //printf("\n\nCOLLIDED!!\n\n");
                force_field[j][0] = force_field[j][0] - (x[i][0]-x[j][0]) * LJForce;
                force_field[j][1] = force_field[j][1] - (x[i][1]-x[j][1]) * LJForce;
                force_field[j][2] = force_field[j][2] - (x[i][2]-x[j][2]) * LJForce;
                if(LJForce * time_step > 5){
                    printf("ERROR: LJ TOO LARGE! Failed for value %lf with distance = %lf\n\n",LJForce,dist_matrix[i][j]);
                    exit(0);
                }
            }   
        }
           // for(int z=0; z<NB; z++){printf("\n\nAND!!: for i = %d: LJ force[%d][0] = %lf\n\n",i,z, force_field[z][0]);};
        if(i > 0){
            rlogarg = 1 - pow(dist_matrix[i][i-1]/R_0, 2);

           if(rlogarg < 0.1){
                printf("\n\nFENE bond too long: atom i = %d and atom j = %d for r = %lf; rlogarg = %lf \n\n", i,i-1,dist_matrix[i][i-1],rlogarg);
                if (rlogarg <= -3.0){
                    printf("\n\nERROR: BAD BOND!\n\n");
                    exit(0);
                }
                rlogarg = 0.1;
            }

            force_field[i][0] = force_field[i][0] - (x[i][0]-x[i-1][0]) * (K_fene / rlogarg);
            force_field[i][1] = force_field[i][1] - (x[i][1]-x[i-1][1]) * (K_fene / rlogarg);
            force_field[i][2] = force_field[i][2] - (x[i][2]-x[i-1][2]) * (K_fene / rlogarg);
                
            //printf("\n\nFENE[%d][%d] = %lf\n\n",i,j, (x[i][0]-x[j][0]) * (K_fene / ( 1 - pow(dist_matrix[i][j]/R_0, 2) )));

            force_field[i-1][0] = force_field[i-1][0] + (x[i][0]-x[i-1][0]) * (K_fene / rlogarg);
            force_field[i-1][1] = force_field[i-1][1] + (x[i][1]-x[i-1][1]) * (K_fene / rlogarg);
            force_field[i-1][2] = force_field[i-1][2] + (x[i][2]-x[i-1][2]) * (K_fene / rlogarg);
        }

    }
    
    // ADD OTHER TWO FORCES 
    force_field[0][0] = force_field[0][0] + k_spring * (x[1][0] - x[0][0]) + a_force * dir_vec[0][0];
    force_field[0][1] = force_field[0][1] + k_spring * (x[1][1] - x[0][1]) + a_force * dir_vec[0][1];
    force_field[0][2] = force_field[0][2] + k_spring * (x[1][2] - x[0][2]) + a_force * dir_vec[0][2] ;

    dx[0] = 1/gamma * ( force_field[0][0] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[0][0],'g');
    dy[0] = 1/gamma * ( force_field[0][1] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[0][1],'g');
    dz[0] = 1/gamma * ( force_field[0][2] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[0][2],'g');

    for(num=1; num<NB-1; num++){
        force_field[num][0] = force_field[num][0] + k_spring * (x[num+1][0] + x[num-1][0] - 2 * x[num][0]) + a_force * dir_vec[num][0];
        force_field[num][1] = force_field[num][1] + k_spring * (x[num+1][1] + x[num-1][1] - 2 * x[num][1]) + a_force * dir_vec[num][1];
        force_field[num][2] = force_field[num][2] + k_spring * (x[num+1][2] + x[num-1][2] - 2 * x[num][2]) + a_force * dir_vec[num][2];
       
        dx[num] = 1/gamma * (force_field[num][0] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[num][0],'g');
        dy[num] = 1/gamma * (force_field[num][1] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[num][1],'g');
        dz[num] = 1/gamma * (force_field[num][2] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[num][2],'g');
    }

    force_field[NB-1][0] = force_field[NB-1][0] + k_spring * (x[NB-2][0] - x[NB-1][0]) + a_force * dir_vec[NB-1][0];
    force_field[NB-1][1] = force_field[NB-1][1] + k_spring * (x[NB-2][1] - x[NB-1][1]) + a_force * dir_vec[NB-1][1];
    force_field[NB-1][2] = force_field[NB-1][2] + k_spring * (x[NB-2][2] - x[NB-1][2]) + a_force * dir_vec[NB-1][2];

    dx[NB-1] = 1/gamma * (force_field[NB-1][0] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[NB-1][0], 'g');
    dy[NB-1] = 1/gamma * (force_field[NB-1][1] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[NB-1][1], 'g');
    dz[NB-1] = 1/gamma * (force_field[NB-1][2] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[NB-1][2], 'g');


    for(i=0; i<NB; i++){
        x_mean[i][0] = x[i][0] + dx[i];
        x_mean[i][1] = x[i][1] + dy[i];
        x_mean[i][2] = x[i][2] + dz[i];
    }
    
    // REAL FORCE 
    distance_matrix(*x_mean, *dist_matrix, NB);

    for(i=0; i<NB; i++){
        //printf("\n\nforce_field[%d][0] = %lf\n\n",i, force_field[i][0]);
        for(j=i; j<NB; j++){
            if(i==j){continue;}
            
            if(dist_matrix[i][j] <= pow(2.,1./6.) * sigma){
                LJForce = epsilon / dist_matrix[i][j] * (48. * (pow(sigma,12.) / pow(dist_matrix[i][j],13.) - 0.5 * pow(sigma,6.) / pow(dist_matrix[i][j],7.)));

                force_field[i][0] = force_field[i][0] + (x_mean[i][0]-x_mean[j][0]) * LJForce;
                force_field[i][1] = force_field[i][1] + (x_mean[i][1]-x_mean[j][1]) * LJForce;
                force_field[i][2] = force_field[i][2] + (x_mean[i][2]-x_mean[j][2]) * LJForce;
                //printf("\n\nj = %d; LJ force[%d][0] = %lf\n\n",j,i,   force_field[i][0]);
                //printf("\n\nCOLLIDED!!\n\n");
                force_field[j][0] = force_field[j][0] - (x_mean[i][0]-x_mean[j][0]) * LJForce;
                force_field[j][1] = force_field[j][1] - (x_mean[i][1]-x_mean[j][1]) * LJForce;
                force_field[j][2] = force_field[j][2] - (x_mean[i][2]-x_mean[j][2]) * LJForce;
                if(LJForce * time_step > 5){
                    printf("ERROR: LJ TOO LARGE FOR X BAR: Failed for value %lf with distance = %lf\n\n",LJForce,dist_matrix[i][j]);
                    exit(0);
                }
            }   
        }
           // for(int z=0; z<NB; z++){printf("\n\nAND!!: for i = %d: LJ force[%d][0] = %lf\n\n",i,z, force_field[z][0]);};
        if(i > 0){
            rlogarg = 1 - pow(dist_matrix[i][i-1]/R_0, 2);

           if(rlogarg < 0.1){
                printf("\n\nFENE bond too long: atom i = %d and atom j = %d for r = %lf; rlogarg = %lf \n\n", i,i-1,dist_matrix[i][i-1],rlogarg);
                if (rlogarg <= -3.0){
                    printf("\n\nERROR: BAD BOND!\n\n");
                    exit(0);
                }
                rlogarg = 0.1;
            }

            force_field[i][0] = force_field[i][0] - (x_mean[i][0]-x_mean[i-1][0]) * (K_fene / rlogarg);
            force_field[i][1] = force_field[i][1] - (x_mean[i][1]-x_mean[i-1][1]) * (K_fene / rlogarg);
            force_field[i][2] = force_field[i][2] - (x_mean[i][2]-x_mean[i-1][2]) * (K_fene / rlogarg);
                
            force_field[i-1][0] = force_field[i-1][0] + (x_mean[i][0]-x_mean[i-1][0]) * (K_fene / rlogarg);
            force_field[i-1][1] = force_field[i-1][1] + (x_mean[i][1]-x_mean[i-1][1]) * (K_fene / rlogarg);
            force_field[i-1][2] = force_field[i-1][2] + (x_mean[i][2]-x_mean[i-1][2]) * (K_fene / rlogarg);
        }
    }
    
    force_field[0][0] = force_field[0][0] + k_spring * (x_mean[1][0] - x_mean[0][0]) + a_force * dir_vec[0][0];
    force_field[0][1] = force_field[0][1] + k_spring * (x_mean[1][1] - x_mean[0][1]) + a_force * dir_vec[0][1];
    force_field[0][2] = force_field[0][2] + k_spring * (x_mean[1][2] - x_mean[0][2]) + a_force * dir_vec[0][2] ;

    dx[0] = 1/gamma * (0.5 * force_field[0][0] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[0][0],'g');
    dy[0] = 1/gamma * (0.5 * force_field[0][1] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[0][1],'g');
    dz[0] = 1/gamma * (0.5 * force_field[0][2] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[0][2],'g');

    for(num=1; num<NB-1; num++){
        force_field[num][0] = force_field[num][0] + k_spring * (x_mean[num+1][0] + x_mean[num-1][0] - 2 * x_mean[num][0]) + a_force * dir_vec[num][0];
        force_field[num][1] = force_field[num][1] + k_spring * (x_mean[num+1][1] + x_mean[num-1][1] - 2 * x_mean[num][1]) + a_force * dir_vec[num][1];
        force_field[num][2] = force_field[num][2] + k_spring * (x_mean[num+1][2] + x_mean[num-1][2] - 2 * x_mean[num][2]) + a_force * dir_vec[num][2];
       
        dx[num] = 1/gamma * (0.5 * force_field[num][0] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[num][0],'g');
        dy[num] = 1/gamma * (0.5 * force_field[num][1] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[num][1],'g');
        dz[num] = 1/gamma * (0.5 * force_field[num][2] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[num][2],'g');
    }

    force_field[NB-1][0] = force_field[NB-1][0] + k_spring * (x_mean[NB-2][0] - x_mean[NB-1][0]) + a_force * dir_vec[NB-1][0];
    force_field[NB-1][1] = force_field[NB-1][1] + k_spring * (x_mean[NB-2][1] - x_mean[NB-1][1]) + a_force * dir_vec[NB-1][1];
    force_field[NB-1][2] = force_field[NB-1][2] + k_spring * (x_mean[NB-2][2] - x_mean[NB-1][2]) + a_force * dir_vec[NB-1][2];

    dx[NB-1] = 1/gamma * (0.5 * force_field[NB-1][0] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[NB-1][0], 'g');
    dy[NB-1] = 1/gamma * (0.5 * force_field[NB-1][1] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[NB-1][1], 'g');
    dz[NB-1] = 1/gamma * (0.5 * force_field[NB-1][2] * time_step) + sqrt(2*pos_diff*time_step) * ran2(&randx_mat[NB-1][2], 'g');

}