#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ran2.h"
#include "myparam.h"


void adopt_new_pos(double** x, double** next_pt);


void adopt_new_pos(double** x, double** next_pt){
    int i;


    for(i=0; i<NB; i++){
        x[i][0] = next_pt[i][0];
        x[i][1] = next_pt[i][1];
        x[i][2] = next_pt[i][2];
    }
}