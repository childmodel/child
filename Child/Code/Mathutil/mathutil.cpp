/**************************************************************************/
/**
**  @file mathutil.cpp
**  @brief Special math routines not in math libraries. Most or all
**         from Numerical Recipes in C by Press et al.
**
**  $Id: mathutil.cpp,v 1.4 2003-07-15 17:24:57 childcvs Exp $
*/
/**************************************************************************/

#include "mathutil.h"


/*********************************************************\
**  ran3
**
**  Random number generator from Numerical Recipes.
**  Returns a uniform random number between 0.0 and 1.0.
**  Set idum to any negative value to initialize or
**  reinitialize the sequence.
**
**  Parameters: idum - random seed
**
\*********************************************************/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
   static int inext,inextp;
   static long ma[56];
   static int iff=0;
   long mj,mk;
   int i,ii,k;
   
   if (*idum < 0 || iff == 0) {
      iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
          for (i=1;i<=55;i++) {
             ma[i] -= ma[1+(i+30) % 55];
             if (ma[i] < MZ) ma[i] += MBIG;
          }
      inext=0;
      inextp=31;
      *idum=1;
   }
   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;
   return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#include "mt19937ar-cok.cpp"
