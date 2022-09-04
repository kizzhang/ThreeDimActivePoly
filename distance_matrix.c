#include <stdio.h>
#include <math.h>
#include "myparam.h"

double seuclidean_distance(double *u, double *v, const int n);
int cdist_seuclidean(double *A, double *B, double *dm, const int num_rowsA, const int num_rowsB, const int num_cols);

int cdist_seuclidean(double *A, double *B,
                 double *dm, const int num_rowsA, const int num_rowsB,
                 const int num_cols)
{
    int i, j;

    for (i = 0; i < num_rowsA; ++i){
        double *u = A + (num_cols * i);
        for (j = 0; j < num_rowsB; ++j, ++dm) {
            double *v = B + (num_cols * j);
            *dm = seuclidean_distance(u, v, num_cols);
        }
    }
    return 0;
}

double seuclidean_distance(double *u, double *v, const int n)
{
    double s = 0.0;
    int i;

    for (i = 0; i < n; ++i) {
        const double d = u[i] - v[i];
        s += (d * d) ;
    }
    return sqrt(s);
}


void distance_matrix(double *X, double* distMat, const int row) {
    const int dim = 3;
    int i, j;

    cdist_seuclidean(X, X, distMat, row, row, dim);
     
}


int main () {
    int i,j;
    const int col = 3;
    double dist_Mat[test][test];

    double myArray[test][3];// = { {1, 1, 1}, {2 , 2, 2}, {3, 3, 3}, {4, 4, 4} };
    printf("11myArray is \n");
    for(i=0; i<col; i++){
      myArray[0][i] = 1;
      myArray[1][i] = 2;
      myArray[2][i] = 3;
      myArray[3][i] = 4;
    }
    
    distance_matrix(*myArray, *dist_Mat, test);

    for (i = 0; i < test; i++) {
      for (j = 0; j < test; j++) {

        printf("%f ", dist_Mat[i][j]);
        if (j == test-1) {
         printf("\n");
        }
      }
    }
}