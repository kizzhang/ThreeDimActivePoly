#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"

void initialize_polymer(double x[][3], double v_vec[][3], double dir_vec[][3], double angle[][2]){
    double theta, phi;
    int i; 
    long dum = -100;
    long *dvar;
    dvar=&dum;
    
    for(i=0; i < NB; i++){
        x[i][0] = (double) i; // A chain on x axis spaced by 1
        x[i][1] = 0.;
        x[i][2] = 0.;

        phi =  ran2(dvar,'u')* 2. * PI;
        theta = acos(1.-2.*ran2(dvar,'u'));
 
        angle[i][0] = phi;
        angle[i][1] = theta;

        dir_vec[i][0] = sin(theta) * cos(phi);
        dir_vec[i][1] = sin(theta) * sin(phi);
        dir_vec[i][2] = cos(theta);

        v_vec[i][0] = walk_speed * sin(theta) * cos(phi);
        v_vec[i][1] = walk_speed * sin(theta) * sin(phi);
        v_vec[i][2] = walk_speed * cos(theta);
    }
}