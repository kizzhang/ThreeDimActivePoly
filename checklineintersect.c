#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int check_line_intersect(double* point1, double* point2, double* point3, double* point4);
void Cross_product(double vect_A[], double vect_B[],double cross_P[]);
double Dot_product(double* vec1, double* vec2);

void Cross_product(double vect_A[], double vect_B[],double cross_P[])
{
    //double * myArray = malloc(sizeof(int) * 3);

    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}
int check_line_intersect(double* point1, double* point2, double* point3, double* point4){
    double coplane_det;
    double s;
    double a_vec[3];
    double b_vec[3];
    double c_vec[3];
    double cross_vec1[3];
    double cross_vec2[3];
    for(int i=0;  i <3;i++){
        a_vec[i] = point2[i] - point1[i];
        b_vec[i] = point4[i] - point3[i];
        c_vec[i] = point3[i] - point1[i];
    }
    Cross_product(a_vec,b_vec,cross_vec1);
    
    coplane_det = Dot_product(c_vec, cross_vec1);
    if(coplane_det != 0){return 0;}
    Cross_product(c_vec,b_vec,cross_vec2); 
    s = Dot_product(cross_vec2,cross_vec1) / Dot_product(cross_vec1,cross_vec1);
    if(s > 0 && s < 1){
        printf("s is %lf",s);
        return 1;
    }
    else{
        return 0;
    }
}
