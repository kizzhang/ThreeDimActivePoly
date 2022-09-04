#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myparam.h"

double pos[NBmax][dd];
int    simplify_conf(double xin[NBmax][dd], int Lin, double xsimple[NBmax][dd]);
int    rm_node(int nodetoindex[NBmax], int Nin, int nodetoindex_new[NBmax]);
int    judge_line_plane_cross(double planeijk[4],double x1[dd], double x2[dd], double A[dd], double B[dd], double C[dd]);
int    judge_point_in_triangle(double P[dd], double A[dd], double B[dd], double C[dd]);
void   cal_plane(double I[dd], double J[dd], double K[dd], double planeijk[4]);
void   cross_product(double a[dd],double b[dd], double res[dd]);
double dot_product(double a[dd], double b[dd]);


int simplify_conf(double xin[NBmax][dd], int Lin, double xsimple[NBmax][dd]) {
   int i,j, d, Nin,  flag_success_remove_node, nodetoindex[NBmax], nodetoindex_new[NBmax];
   
   Nin                       =  Lin;
   flag_success_remove_node  =  1;
   for(i=0;i<Nin;i++)  nodetoindex[i] = i;
   for(i=0;i<Nin;i++)  for(d=0;d<dd;d++) pos[i][d] = xin[i][d];

   while( Nin >= 6 && flag_success_remove_node==1) {
      flag_success_remove_node = rm_node(nodetoindex, Nin, nodetoindex_new);
      if(flag_success_remove_node==1) {
         Nin -= 1;
         for(i=0;i<Nin;i++) nodetoindex[i] = nodetoindex_new[i];
      }
   }

   for(i=0;i<Nin;i++) {
      j      = nodetoindex[i];
      for(d=0;d<dd;d++) xsimple[i][d] = pos[j][d];
   }

   return Nin;
}

int rm_node(int nodetoindex[NBmax], int Nin, int nodetoindex_new[NBmax]) {
   int inode, jnode,line1, inode_remove,line2, tri1, tri2, tri3;
   int flag_success_remove_node,  flag_cross;
   double planeijk[4];

   for(inode=2;inode<Nin;inode++) { 
       tri1 =  nodetoindex[inode-2];
       tri2 =  nodetoindex[inode-1];
       tri3 =  nodetoindex[inode];

       cal_plane( pos[tri1], pos[tri2], pos[tri3], planeijk );

       for(jnode=1;jnode<Nin;jnode++) {
          line1 = nodetoindex[jnode-1];
          line2 = nodetoindex[jnode];
          if(line1==tri1 || line1==tri2 || line1==tri3 || line2==tri1 || line2==tri2 || line2==tri3) continue;
          flag_cross = judge_line_plane_cross(planeijk,pos[line1],pos[line2],pos[tri1],pos[tri2],pos[tri3]);
	  // printf("======= line %3d %3d; triangle %3d %3d %3d flag_cross %3d\n", line1,line2,tri1,tri2, tri3,flag_cross);
          if(flag_cross==1) break; 
        }    
	
	if(flag_cross == 0) break; // this triangle #inode has no crossing with any segment
   }
   inode_remove = inode - 1;

   if(flag_cross == 0) {
       // printf("remove inode_remove %5d index_in_chain %5d\n", inode_remove, nodetoindex[inode_remove]);
       jnode    = 0;
       for(inode=0;inode<Nin;inode++) {
	   if(inode==inode_remove) continue;
	   nodetoindex_new[jnode] = nodetoindex[inode];
	   jnode++;
       }
       if(jnode != (Nin-1)) {fprintf(stderr,"Nout should be %d, but it is %d\n",Nin-1,jnode); exit(-1);}
       flag_success_remove_node = 1;
   } else {
       flag_success_remove_node = 0;
   }

   return flag_success_remove_node;
}

int judge_line_plane_cross(double planeijk[4],double x1[dd], double x2[dd], double A[dd], double B[dd], double C[dd]) {
   int  d;
   double  d1, d2, lambda, dis, cpot[dd];


   d1 = planeijk[3]; for(d=0;d<dd;d++) d1 += planeijk[d] * x1[d];
   d2 = planeijk[3]; for(d=0;d<dd;d++) d2 += planeijk[d] * x2[d];

   
   if(d1>0. && d2>0.) return 0; // two points are located in the same side of plane
   if(d1<0. && d2<0.) return 0;

   dis = fabs(d1-d2);

   lambda = fabs(d1)/dis;
   for(d=0;d<dd;d++) cpot[d] = x1[d] + (x2[d] - x1[d]) * lambda;
   if( judge_point_in_triangle(cpot, A, B, C) == 1 ) return 1; // line-triangle crossing point is inside triangle
   return 0; 
}

int judge_point_in_triangle(double P[dd], double A[dd], double B[dd], double C[dd]) {
   int d;
   double  v0[dd], v1[dd], v2[dd], dot00, dot01, dot02, dot11, dot12, inverDeno, u,v;
   // use a method copied from internet; this method has been validated

   for(d=0;d<dd;d++) v0[d] = C[d] - A[d];
   for(d=0;d<dd;d++) v1[d] = B[d] - A[d];
   for(d=0;d<dd;d++) v2[d] = P[d] - A[d];

   dot00 = dot_product(v0, v0);
   dot01 = dot_product(v0, v1);
   dot02 = dot_product(v0, v2);
   dot11 = dot_product(v1, v1);
   dot12 = dot_product(v1, v2);
   inverDeno = 1. / (dot00 * dot11 - dot01 * dot01) ;
   u = (dot11 * dot02 - dot01 * dot12) *inverDeno ;
   if (u < 0 || u > 1) return 0;
   v = (dot00 * dot12 - dot01 * dot02) *inverDeno ;
   if (v < 0 || v > 1) return 0;
   if((u+v) <= 1.) return 1;
   return 0;
}

void cal_plane(double I[dd], double J[dd], double K[dd], double planeijk[4]) {
   int d;
   double  vij[dd], vik[dd], vtmp[dd];

   for(d=0;d<dd;d++) vij[d] = J[d] - I[d];
   for(d=0;d<dd;d++) vik[d] = K[d] - I[d];
   cross_product(vij, vik, vtmp);
   for(d=0;d<dd;d++) planeijk[d] = vtmp[d];
   planeijk[3]  = 0.;
   for(d=0;d<dd;d++)  planeijk[3] += I[d]*planeijk[d];
   planeijk[3]  = -planeijk[3];

}


void cross_product(double a[dd],double b[dd], double res[dd]) {
   res[0] = a[1]*b[2] - a[2]*b[1];
   res[1] = a[2]*b[0] - a[0]*b[2];
   res[2] = a[0]*b[1] - a[1]*b[0];
}

double dot_product(double a[dd], double b[dd]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}



