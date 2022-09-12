#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myparam.h"

extern double knottable[Nknottypetotal][3];

double aglobal[MAXcross][MAXcross];

int flag_close_type_global;


int    cal_cross_point(double A[2], double B[2], double C[2], double D[2], double k[2]);
void   sort_pos(double p[MAXcross][2], int crosstype[MAXcross], int Ncross);
double cal_array_determinant( int N);

void   add_ends(double x[NBmax][dd], int Lchain, double xnew[NBmax][dd]);

double knot_Apoly(double x[NBmax][3], int N, double tvalue);
int    get_knottype(double x[NBmax][3], int N );
double max(double a[], int N);

int get_knottype_of_entire_chain(double x[NBmax][dd], int Lchain);
int simplify_conf(double xin[NBmax][dd], int Lin, double xsimple[NBmax][dd]);
int hull_ends(double x[NBmax][dd],  int L, double point1[dd], double point2[dd] );


int  add_ends_by_qhull(double x[NBmax][dd], int Lchain, double xnew[NBmax][dd] ) {
    int flag_close_type, i, d;
    double point1[dd], point2[dd];

    // printf("OK2\n");

    flag_close_type =  hull_ends(x, Lchain,  point1, point2) ;

    //printf("type %5d\n", flag_close_type);

    if(flag_close_type == 1) {
       for(d=0;d<dd;d++) {
          xnew[0][d]         =  point1[d];
          xnew[Lchain+1][d]  =  point2[d];
       }

       for(i=0;i<Lchain;i++) for(d=0;d<dd;d++) xnew[i+1][d] = x[i][d];
    } else { // it is shorter to close two ends by a line
       for(i=0;i<Lchain;i++) for(d=0;d<dd;d++) xnew[i][d] = x[i][d];
       for(d=0;d<dd;d++) xnew[Lchain][d] = x[0][d];
    }

    return flag_close_type;


}


int get_knottype_of_entire_chain(double x[NBmax][dd], int Lchain) {
   int j,k,d,  Lchain_simple, knottype;
   double xnew[NBmax][dd], xsimple[NBmax][dd];

   

   if(flag_connect_ends_by_a_line == 1) {
      x[NB][0] = x[0][0];  x[NB][1] = x[0][1]; x[NB][2] = x[0][2]; // used for a ring polymer
      printf("Finished!\n\n");
      knottype = get_knottype(x, NB+1);

   } else {
       // printf("OK1\n");
       if(Lchain<=5) return 0; 

       flag_close_type_global  = add_ends_by_qhull(x, Lchain, xnew );
       //printf("flag_close_type %5d\n", flag_close_type_global);
      
       if(flag_close_type_global==0) {
           Lchain_simple  =  simplify_conf(xnew, Lchain+1,xsimple);
           knottype       =  get_knottype(xsimple, Lchain_simple);
       } else {
           Lchain_simple  =  simplify_conf(xnew, Lchain+2,xsimple);
           knottype       =  get_knottype(xsimple, Lchain_simple);
       }
   }
   
   return knottype;

}




int get_knottype(double x[NBmax][3], int N ){
   int i, knottype;
   double topA, topB, Apoly, small_value;

   small_value = 0.001; // to allow numerical error in calculating Apoly

   topA  = knot_Apoly(x, N, -0.5);
   topB  = knot_Apoly(x, N, -2.0);
   Apoly = topA * topB;

   // printf("Apoly %10.3f\n", Apoly);

   knottype = 999999; // unknown knotype
   for(i=0;i<Nknottypetotal;i++) {
      if(fabs( Apoly - knottable[i][2]) < small_value) {
         knottype = (int)(knottable[i][1] + 0.1);
	 break;
      }
   }
   return knottype;
}

double knot_Apoly(double x[NBmax][3], int N, double tvalue) {
   int ierr, d;
   double topl;
   int flag_cross,i,j,m, icross,Ncross, cr[MAXcross], crosstype[MAXcross];
   double A[2],B[2],C[2],D[2],k[2], p[MAXcross][2], pAB, pCD, zAB,zCD, vAB[2], vCD[2], zABCD;
   double t1,t2,t3,t4;


   icross = 0;
   for(i=0;i<(N-1);i++) {
      A[0] = x[i][0];    A[1] = x[i][1];
      B[0] = x[i+1][0];  B[1] = x[i+1][1];
      for(d=0;d<2;d++)  vAB[d]  = B[d] -A[d];
      
      for(j=i+2;j<(N-1);j++) {
         if(flag_connect_ends_by_a_line==1 || flag_close_type_global ==0) { if(i==0 && (j+1)==(N-1)) continue; } // used for a ring polymer 

         C[0] = x[j][0];    C[1] = x[j][1];
	 D[0] = x[j+1][0];  D[1] = x[j+1][1];

	 flag_cross = cal_cross_point(A,B,C,D,k);
      
	 if(flag_cross == 0) continue;
	 pAB = i + k[0]; // contour position
	 pCD = j + k[1];

	 zAB  = x[i][2] + k[0] * (x[i+1][2] - x[i][2]);
	 zCD  = x[j][2] + k[1] * (x[j+1][2] - x[j][2]);

         for(d=0;d<2;d++) vCD[d] = D[d] - C[d];
	 zABCD = vAB[0]*vCD[1] - vAB[1]*vCD[0];

	 if(zAB > zCD) {
	    p[icross][0] = pAB; // records the contour position of the under-going point of a crossing
	    p[icross][1] = pCD; // records the contour position of the over-pass point of a crossing
            if(zABCD < 0) crosstype[icross] = 1; 
	    else          crosstype[icross] = 2;
	 }else{
	    p[icross][0] = pCD;
	    p[icross][1] = pAB;
	    if(zABCD < 0) crosstype[icross] = 2;
	    else          crosstype[icross] = 1;
	 }
	 icross++;
	 if(icross > MAXcross) {fprintf(stderr,"increase MAXcross %10d\n",MAXcross); exit(-1);}
      }
   }
   Ncross = icross;

   if(Ncross<3) return 1; // not knotted

   sort_pos(p,crosstype,Ncross); // sort crossings by undergoing positions

   //get cr(i) =k; record the pair (i,k); i is the index of crossing point; k is the index of generator overpassing i-th crossing point
   //note that generator #j is between #(j-1) and #j crossing points
   //the first generator is between the last #(Ncross-1) and first #(0) crossing points,  
   for(i=0;i<Ncross;i++) { 
      cr[i] = 0;
      for(j=0;j<Ncross;j++)  {
         if(p[i][1] > p[j][0]) cr[i] = j+1; // find the last undergoing segment p[j][0] which is less than overpassing point p[i][1]
      }
      if(cr[i]==Ncross) cr[i] = 0;
   }


   for(i=0;i<Ncross;i++) 
      for(j=0;j<Ncross;j++)
          aglobal[i][j] = 0;

   for(i=0;i<Ncross;i++) {
      j = cr[i];
      if(j==i || j==(i+1))  { 
         aglobal[i][i] = -1; 
	 aglobal[i][i+1] = 1;
      } else if(crosstype[i]==1){
         aglobal[i][i]   = 1;
	 aglobal[i][i+1] = -tvalue;
	 aglobal[i][j]   =  tvalue-1;
      } else if(crosstype[i]==2) {
         aglobal[i][i]   = -tvalue;
	 aglobal[i][i+1] = 1;
	 aglobal[i][j]   =  tvalue-1;
      } else{
         fprintf(stderr,"crosstype is %5d\n",crosstype[i]);
	 exit(-1);
      }
   }


   topl = cal_array_determinant(Ncross-1);
   return topl;

}

void sort_pos(double p[MAXcross][2], int crosstype[MAXcross], int Ncross) {
  int i, imin,m;
  double amin, q[MAXcross][2], crosstype_sort[MAXcross], maxnum=1.e10;

  m = 0;
  while(m < Ncross) {
    amin = p[0][0]; imin = 0;
    for(i=1;i<Ncross;i++) {
        if(p[i][0] < amin) {
	   amin = p[i][0];
	   imin = i;
	}
    }
    q[m][0]     = p[imin][0]; q[m][1] = p[imin][1];  crosstype_sort[m] = crosstype[imin];
    p[imin][0]  = maxnum;
    m++;
  }

  for(i=0;i<Ncross;i++) {
    p[i][0] = q[i][0]; p[i][1] = q[i][1];  crosstype[i] = crosstype_sort[i];
  }
     
}

int cal_cross_point(double A[2], double B[2], double C[2], double D[2], double k[2]) {
  // determine the cross-point of line AB and line CD
    double a,b,c,d,e,f, kAB, kCD;
    a    = B[0] - A[0];
    b    = C[0] - D[0];
    c    = C[0] - A[0];

    d    = B[1] - A[1];
    e    = C[1] - D[1];
    f    = C[1] - A[1];

    kAB  = (c*e - b*f)/(a*e-b*d);
    kCD  = (a*f - c*d)/(a*e-b*d);

    k[0] = kAB;
    k[1] = kCD;
    if(kAB>0. && kAB<1. && kCD>0. && kCD <1.)
       return 1; // two lines cross each other
    else
       return 0; // no cross
}

// I don't remember where I got this code to calculate determinant. This code is simple, efficient and appears to be correct
double cal_array_determinant(int N){
   int i,j,m,n,s,t,k=1;
   double f=1,c,x,sn;

   for (i=0,j=0;i<N&&j<N;i++,j++) {
       if(aglobal[i][j]==0) {
          for (m=i;aglobal[m][j]==0;m++);
	  if (m==N) {
	     sn=0;
	     return sn;
	  } else
	  for (n=j;n<N;n++) {
	     c=aglobal[i][n];
	     aglobal[i][n]=aglobal[m][n];
	     aglobal[m][n]=c;
	  }
	  k*=(-1);
      }
      for(s=N-1;s>i;s--) {
         x=aglobal[s][j];   
	 for (t=j;t<N;t++)   aglobal[s][t]-=aglobal[i][t]*(x/aglobal[i][j]);
      }
   }
   
   for (i=0;i<N;i++)  f*=aglobal[i][i];
   sn=k*f;
   return sn;
 
} 

// extend from the center to each end, and add a point along the extension beyond the current conformation such 
// that the connection between two new ends has no crossing with other segments.

void add_ends(double x[NBmax][dd], int Lchain, double xnew[NBmax][dd]) {
   int j, d, Lchainnew;
   double xcen[dd],dis1, dis2, xdeviation[NBmax][dd], xdev_max[dd], xcen_to_end1[dd], xcen_to_end2[dd];
   double frescale1[dd],frescale2[dd], frescale_max1, frescale_max2;


   for(d=0;d<dd;d++) xcen[d] = 0;

   for(j=0;j<Lchain;j++) for(d=0;d<dd;d++)  xcen[d]       += x[j][d];

   for(d=0;d<dd;d++) xcen[d] /= Lchain;

   for(d=0;d<dd;d++) xcen_to_end1[d] = x[0][d] - xcen[d];
   for(d=0;d<dd;d++) xcen_to_end2[d] = x[Lchain-1][d] - xcen[d];

   for(j=0;j<Lchain;j++) for(d=0;d<dd;d++) xdeviation[j][d] = x[j][d] - xcen[d];

   for(d=0;d<dd;d++) xdev_max[d] = fabs(xdeviation[0][d]);

   for(d=0;d<dd;d++)
      for(j=0;j<Lchain;j++) {
         if( fabs(xdeviation[j][d]) > xdev_max[d] ) xdev_max[d] = fabs(xdeviation[j][d]);
      }

   for(d=0;d<dd;d++) frescale1[d] = 2*xdev_max[d]/(fabs(xcen_to_end1[d]));
   for(d=0;d<dd;d++) frescale2[d] = 2*xdev_max[d]/(fabs(xcen_to_end2[d]));

   frescale_max1 = max(frescale1,3);
   frescale_max2 = max(frescale2,3);


   Lchainnew = Lchain + 2;
   for(j=0;j<Lchain;j++)    for(d=0;d<dd;d++) xnew[j+1][d]   = x[j][d];
   for(d=0;d<dd;d++)  xnew[0][d]        = xcen[d] + frescale_max1 * xcen_to_end1[d];
   for(d=0;d<dd;d++)  xnew[Lchain+1][d] = xcen[d] + frescale_max2 * xcen_to_end2[d];

/*
   printf("old conformation\n");
   for(j=0;j<Lchain;j++) printf("%10.3f %10.3f %10.3f\n",x[j][0],  x[j][1], x[j][2]);
   printf("xcen         %10.3f %10.3f %10.3f\n", xcen[0], xcen[1],xcen[2]);
   printf("xcen_to_end1 %10.3f %10.3f %10.3f\n", xcen_to_end1[0], xcen_to_end1[1],xcen_to_end1[2]);
   printf("xcen_to_end2 %10.3f %10.3f %10.3f\n", xcen_to_end2[0], xcen_to_end2[1],xcen_to_end2[2]);
   printf("xdev_max     %10.3f %10.3f %10.3f\n", xdev_max[0], xdev_max[1],xdev_max[2]);
   printf("frescale_max 1st end %10.3f; 2nd end %10.3f\n", frescale_max1, frescale_max2);
   printf("\n");


   for(j=0;j<Lchainnew;j++) printf("%10.3f %10.3f %10.3f\n",xnew[j][0],  xnew[j][1], xnew[j][2]);
   exit(0);
*/
}

double max(double a[], int N) {
  int i;
  double amax = a[0];
  for(i=1;i<N;i++) if(a[i]>amax) amax = a[i];
  return amax;
}



