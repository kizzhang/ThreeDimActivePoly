/*******************************************************************************  
   subroutine returns generated random numbers between 0 and 1  
   [status: checked + confirmed]
********************************************************************************
   Pseudorandom number generator source code: rand2 from Numerical Recipe in C 
*******************************************************************************/
#include "ran2.h"
#include <math.h>
#include <stdio.h>
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI 3.141592653589793

float ran2_gen(long *idum){
/*******************************************************************************
Long period (> 2 x 10^18) random number generator of L�Ecuyer with Bays-Durham 
shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 
1.0 (exclusive of the endpoint values). Call with idum a negative integer to 
initialize; thereafter, do not alter idum between successive deviates in a 
sequence. RNMX should approximate the largest floating value that is less than 
1. 
*******************************************************************************/
   int j;
   long k;
   static long idum2=123456789;
   static long iy=0;
   static long iv[NTAB];
   float temp;
   

   if (*idum <= 0) {                 		// Initialize.
      if (-(*idum) < 1) *idum=1;   			// Be sure to prevent idum = 0.
      else *idum = -(*idum);
      idum2=(*idum);
      for (j=NTAB+7;j>=0;j--) {   			// Load the shuffle table (after 8 warm-ups).
         k=(*idum)/IQ1;
         *idum=IA1*(*idum-k*IQ1)-k*IR1;
         if (*idum < 0) *idum += IM1;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ1; 				           // Start here when not initializing.
   *idum=IA1*(*idum-k*IQ1)-k*IR1;		   // Compute without idum=(IA1*idum) % IM1 
   if (*idum < 0) *idum += IM1;            // overflows by Schrage�s method.
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2; 		   // Compute idum2=(IA2*idum) % IM2 likewise.
   if (idum2 < 0) idum2 += IM2;
   j=iy/NDIV; 				               // Will be in the range 0..NTAB-1.
   iy=iv[j]-idum2; 				           // Here idum is shuffled, idum and idum2 are
   iv[j] = *idum; 				           // combined to generate output.
   if (iy < 1) iy += IMM1;
   if ((temp=AM*iy) > RNMX) return RNMX;   // Because users don�t expect endpoint values.
   else return temp;
}

double ran2(long *idum, char option){
   float ran_num;
   ran_num = ran2_gen(idum);

   if(option ==  'u'){
      return ran_num;
   }
   else if(option == 'g'){
      float ran_num2 = ran2_gen(idum);

      return cos(2*PI*ran_num)*sqrt(-2.*log(ran_num2));
   }  
}

