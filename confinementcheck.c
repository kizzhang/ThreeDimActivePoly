#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myparam.h"

void RaySphereIntersection(double* pt_in, double* pt_out, double* intersection);
double Dot_product(double* vec1, double* vec2);
void Reflection_from_sphere(double* pt_in, double* next_pt, double* reflection);


void RaySphereIntersection(double* pt_in, double* pt_out, double* intersection)
{
    double A;
    double B;
    double C;
    double det;
    double t;
    
    A = pow((*pt_in), 2) + pow(*(pt_in+1), 2) + pow(*(pt_in+2), 2) - pow(radius,2);
    C = pow((*pt_in)-(*pt_out), 2) + pow(*(pt_in+1)-*(pt_out+1), 2) + pow(*(pt_in+2) -*(pt_out+2), 2);
    B = pow((*pt_out), 2) + pow(*(pt_out+1), 2) + pow(*(pt_out+2), 2) - A - C -  pow(radius,2);
   
    det = pow(B,2.) - 4.*A*C;
    
    if(det<0){
        printf("SOMETHING IS WRONG: det is %lf\n",det);
        for(int ii =0; ii<3; ii++){
            printf("pt_in[%d] is %lf\n", ii, pt_in[ii]);
            printf("pt_out[%d] is %lf\n", ii, pt_out[ii]);
        }
        
        printf("ERROR: For reflection from boundary, there is no intersections!\n");
        exit(-1);
    }

    else if(det>0){
        t = (-1*B + sqrt(det)) / 2 / C;
        if(t < 0){t = (-1*B - sqrt(det)) / 2 / C;}
        
        *intersection = (*pt_in) * (1-t) + t * (*pt_out);
        *(intersection+1) = *(pt_in+1) * (1-t) + t * *(pt_out+1);
        *(intersection+2) = *(pt_in+2) * (1-t) + t * *(pt_out+2);
    }
    else{
        printf("Something is wrong: Touching the sphere!\n");
        exit(-1);
    }
}

double Dot_product(double* vec1, double* vec2){
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

void Reflection_from_sphere(double* pt_in, double* next_pt, double* reflection){
    double a[3];
    double pt_new[3];
    double pt_out[3];
    double intersect[3]={0,0,0};
    int passed =0;

    pt_new[0] = pt_in[0]; pt_new[1] = pt_in[1]; pt_new[2] = pt_in[2];
    while(sqrt(pow(next_pt[0],2) + pow(next_pt[1],2) + pow(next_pt[2],2)) > radius){
        RaySphereIntersection(pt_new, next_pt, intersect);

        a[0] = next_pt[0]-intersect[0]; a[1] = next_pt[1]-intersect[1]; a[2] = next_pt[2]-intersect[2];
        
        next_pt[0] = next_pt[0] - 2 * Dot_product(a, intersect) / Dot_product(intersect, intersect) * intersect[0];
        next_pt[1] = next_pt[1] - 2 * Dot_product(a, intersect) / Dot_product(intersect, intersect) * intersect[1];
        next_pt[2] = next_pt[2] - 2 * Dot_product(a, intersect) / Dot_product(intersect, intersect) * intersect[2];

        pt_new[0] = intersect[0]; pt_new[1] = intersect[1]; pt_new[2] = intersect[2];

        reflection[0] = a[0] - 2 * Dot_product(a, intersect) / Dot_product(intersect, intersect) * intersect[0];
        reflection[1] = a[1] - 2 * Dot_product(a, intersect) / Dot_product(intersect, intersect) * intersect[1];
        reflection[2] = a[2] - 2 * Dot_product(a, intersect) / Dot_product(intersect, intersect) * intersect[2];
        passed = 1;
    }
    if(passed == 0){
        printf("ERROR: Reflection from surface is not processed. Check the radius of the sphere.\n");
    	exit(-1);
    }
    double norm = sqrt(Dot_product(reflection, reflection));
    
    reflection[0] = reflection[0] / norm;
    reflection[1] = reflection[1] / norm;
    reflection[2] = reflection[2] / norm;
    /*printf("reflection vector is is vx = %lf \t", reflection[0]);
    printf("reflection vector is is vx = %lf \t", reflection[1]);
    printf("reflection vector is is vx = %lf \n", reflection[2]);*/
 
}

