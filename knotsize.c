
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myparam.h"

extern int flag_close_type_global;


int    get_knottype_of_subchain(double x[NBmax][dd], int ileft, int iright);
void   get_knotcore(double x[NBmax][dd], int Lchain, int knottype, int result[3]) ;
void   get_knotcore_circular(double x[NBmax][dd], int Lchain, int knottype_init, int knotsize[3] );
void   get_knotcore_faster(double x[NBmax][dd], int Lchain, int knottype, int result[3]); 
int    get_knottype(double x[NBmax][3], int N );
void   add_ends(double x[NBmax][dd], int Lchain, double xnew[NBmax][dd]);
int    simplify_conf(double xin[NBmax][dd], int Lin, double xsimple[NBmax][dd]);
int    add_ends_by_qhull(double x[NBmax][dd], int Lchain, double xnew[NBmax][dd] ) ;


void get_knotcore_circular(double x[NBmax][dd], int Lchain,  int knottype_init, int knotsize[3] ) {
   int d, ind, ind_final, ishift, i, j, k,ib, Lknot[4][3], Lknot_final[3], knotsize_tmp[3];
   double xshift[NBmax][dd];

   for(ind=0; ind<4; ind++) { // because I don't know whether to cut the circular chain before knot size determination, I try 4 sites
      ishift = (int)((double)(ind * Lchain)/4.);
      for(ib=0;ib<Lchain;ib++) {
         k = ib + ishift;
         if(k>=Lchain) k -= Lchain;
         for(d=0;d<dd;d++) xshift[ib][d] = x[k][d];
       }

       if(flag_get_knotcore_faster == 1) {
          get_knotcore_faster(xshift, Lchain, knottype_init, knotsize_tmp);
       } else {
          get_knotcore(xshift, Lchain, knottype_init, knotsize_tmp);
       }

       for(j=0;j<3;j++) Lknot[ind][j] = knotsize_tmp[j];

    }


    ind_final = 0;
    for(ind=0;ind<4;ind++) if(Lknot[ind][2] < Lknot[ind_final][2]) ind_final = ind;

    for(j=0;j<3;j++)   Lknot_final[j] = Lknot[ind_final][j];


    Lknot_final[0] =   Lknot_final[0] +  (int)((double)(ind_final * Lchain)/4.);
    Lknot_final[1] =   Lknot_final[1] +  (int)((double)(ind_final * Lchain)/4.);

    if(Lknot_final[0] > Lchain) { Lknot_final[0] -= Lchain;   Lknot_final[1] -= Lchain;  }

    for(j=0;j<3;j++) knotsize[j] = Lknot_final[j];


}


int get_knottype_of_subchain(double x[NBmax][dd], int ileft, int iright) {
   int j,k,d, Lchain, Lchain_simple, knottype;
   double  xcore[NBmax][dd], xnew[NBmax][dd], xsimple[NBmax][dd];


   k = 0;
   for(j=(ileft-1);j<iright;j++) {
      for(d=0;d<dd;d++) {
         xcore[k][d] = x[j][d];
      }
      k++;
   }
   Lchain = k;

   if(Lchain<=5) return 0;

   flag_close_type_global  = add_ends_by_qhull(xcore, Lchain, xnew );


   if(flag_close_type_global==0) {
        Lchain_simple  =  simplify_conf(xnew, Lchain+1,xsimple);
        knottype       =  get_knottype(xsimple, Lchain_simple);
   } else {
        Lchain_simple  =  simplify_conf(xnew, Lchain+2,xsimple);
        knottype       =  get_knottype(xsimple, Lchain_simple);
   }

   return knottype;

}


void get_knotcore(double x[NBmax][dd], int Lchain, int knottype, int result[3]) {
   int d, Atmp, ileft, iright, icen, Lleft,Lright,Lleft2,Lright2, knotsizeA, knotsizeB;


   // Part A1: determine left boundary first
   ileft = 1;
   Atmp = get_knottype_of_subchain(x, ileft,Lchain);
   while( Atmp == knottype) { ileft++; Atmp  = get_knottype_of_subchain(x, ileft,Lchain);  }
   Lleft = ileft-1;
   if(Lleft<1) Lleft = 1;


   // Part A2: determine right boundary
   ileft = Lleft; iright = Lchain;
   Atmp = get_knottype_of_subchain(x,Lleft,Lchain);
   while( Atmp == knottype) { iright--; Atmp  = get_knottype_of_subchain(x,Lleft,iright); }
   Lright = iright + 1;
   if(Lright>Lchain) Lright = Lchain;


   // Part B1: determine right boundary first
   iright = Lchain;
   Atmp = get_knottype_of_subchain(x,1,iright);
   while( Atmp == knottype) { iright--; Atmp  = get_knottype_of_subchain(x,1,iright); }
   Lright2 = iright + 1;
   if(Lright2>Lchain) Lright2 = Lchain;


   // Part B2: determine left boundary
   ileft = 1;
   Atmp = get_knottype_of_subchain(x,ileft,Lright2);
   while( Atmp == knottype) { ileft++; Atmp  = get_knottype_of_subchain(x,ileft,Lright2); }
   Lleft2 = ileft - 1;
   if(Lleft2<1) Lleft2 = 1;

   knotsizeA = Lright  - Lleft  + 1;
   knotsizeB = Lright2 - Lleft2 + 1;

   if(knotsizeA < knotsizeB) {
      result[0] = Lleft;
      result[1] = Lright;
      result[2] = knotsizeA;
   } else {
      result[0] = Lleft2;
      result[1] = Lright2;
      result[2] = knotsizeB;
   }
}


void get_knotcore_faster(double x[NBmax][dd], int Lchain, int knottype, int result[3]) {
   int d, Atmp, ileft, iright, icen, Lleft,Lright,Lleft2,Lright2, knotsizeA, knotsizeB;


   // Part A1: determine left boundary first
   ileft = 1; iright = Lchain; icen = (int)((ileft+iright)/2.);  Atmp  = get_knottype_of_subchain(x,icen,Lchain);
   while((iright-ileft)>2) {
      if(Atmp != knottype) { iright = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,icen,Lchain); }
      else               { ileft  = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,icen,Lchain); }
   }
   Atmp = get_knottype_of_subchain(x,ileft,Lchain);
   while( Atmp == knottype) { ileft  = ileft + 1; Atmp  = get_knottype_of_subchain(x,ileft,Lchain); }
   Lleft = ileft-1;

   // Part A2: determine right boundary
   ileft = Lleft; iright = Lchain; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,Lleft,icen);
   while((iright-ileft)>2) {
      if(Atmp != knottype) { ileft  = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,Lleft,icen);}
      else               { iright = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,Lleft,icen);}
   }
   Atmp = get_knottype_of_subchain(x,Lleft,iright);
   while( Atmp == knottype) { iright = iright - 1; Atmp  = get_knottype_of_subchain(x,Lleft,iright); }
   Lright = iright + 1;
   // Part B1: determine right boundary first
   ileft = 1; iright = Lchain; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,1,icen);
   while((iright-ileft)>2) {
      if(Atmp != knottype) { ileft  = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,1,icen);}
      else               { iright = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,1,icen);}
   }
   Atmp = get_knottype_of_subchain(x,1,iright);
   while( Atmp == knottype) { iright = iright - 1; Atmp  = get_knottype_of_subchain(x,1,iright); }
   Lright2 = iright + 1;

   // Part B2: determine left boundary
   ileft = 1; iright = Lright2; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,icen,Lright2);
   while((iright-ileft)>2) {
      if(Atmp != knottype) { iright = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,icen,Lright2);}
      else               { ileft  = icen; icen = (int)((ileft+iright)/2.); Atmp = get_knottype_of_subchain(x,icen,Lright2);}
   }
   Atmp = get_knottype_of_subchain(x,ileft,Lright2);
   while( Atmp == knottype) { ileft  = ileft + 1; Atmp  = get_knottype_of_subchain(x,ileft,Lright2); }
   Lleft2 = ileft - 1;


   knotsizeA = Lright  - Lleft  + 1;
   knotsizeB = Lright2 - Lleft2 + 1;


   if(knotsizeA < knotsizeB) {
      result[0] = Lleft;
      result[1] = Lright;
      result[2] = knotsizeA;
   } else {
      result[0] = Lleft2;
      result[1] = Lright2;
      result[2] = knotsizeB;
   }

}



