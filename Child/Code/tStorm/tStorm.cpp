/*
**  tStorm.cpp
** 
**  Functions for tStorm objects.
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm. It is essentially an implementation of the model of
**  P. Eagleson, 1978b, Water Resources Research. Its services include
**  reading the necessary parameters from a tInputFile, generating a new      
**  storm, and reporting its various values.
**    If you want to provide an option for NOT having storms vary
**  randomly, you can do so by setting optVariable to zero on initialization.
**    The GammaDev() function is provided for future reference; it is not
**  actually used in version 1.0.
**    tStorm objects could be easily modified (or inherited from) to use
**  different distributions. They can also be modified to create objects
**  for other random processes such as river flows, etc.
**    The random number generation routine ran3() from Numerical Recipes
**  is included in this file.
**
**  Version 1.0, Greg Tucker, November 1997.
**  $Id: tStorm.cpp,v 1.1 1998-01-14 20:46:05 gtucker Exp $
*/
#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <math.h>

#include "../Definitions.h"
#include "../Classes.h"
#include "../GlobalFns.h"
#include "../tInputFile/tInputFile.h"
#include "../tRunTimer/tRunTimer.h"

#include "tStorm.h"

//double ran3( long * );  // Random number generator from Numerical Recipes in C


/*
**  tStorm::tStorm:  Constructor for storms. The default constructor
**                   assigns a value of unity to storm depth, duration,
**                   and interstorm duration.
*/
tStorm::tStorm( int optvar )
{
   optVariable = optvar;
   pMean = 1.0;
   stdurMean = 1.0;
   istdurMean = 1.0;
   p = 1.0;
   stdur = 1.0;
   istdur = 1.0;
   srand( 0 );
}
/*
**  tStorm::tStorm:  The constructor that's really used assigns values for
**                   mean storm depth, duration, and interstorm duration,
**                   initializes current depth, etc, to the mean values,
**                   and initializes the random number generator.
*/
tStorm::tStorm( float mp, float ms, float mis, unsigned sd, int optvar )
{
   optVariable = optvar;
   pMean = mp;
   stdurMean = ms;
   istdurMean = mis;
   p = pMean;
   stdur = stdurMean;
   istdur = istdurMean;
   seed = sd;
   srand( seed );
}


/*
**  tStorm::tStorm
**
**  Alternative constructor that reads parameters directly from a tInputFile
**  object. Reads option for variable storms (normally this is "yes"---that's
**  the point of these objects---but a user may wish to switch off variation
**  as a test), mean values for rainfall intensity, duration, and interstorm
**  period, and a random seed to initialize the random number generator.
*/
tStorm::tStorm( tInputFile &infile )
{
   optVariable = infile.ReadItem( optVariable, "OPTVAR" );
   pMean = infile.ReadItem( pMean, "PMEAN" );
   stdurMean = infile.ReadItem( stdurMean, "STDUR" );
   istdurMean = infile.ReadItem( istdurMean, "ISTDUR" );
   p = pMean;
   stdur = stdurMean;
   istdur = istdurMean;
   seed = infile.ReadItem( seed, "SEED" );
   srand( seed );
}


/*
**  GenerateStorm
**
**  Generates a new storm by drawing new values of p, stdur, and istdur from
**  an exponential distribution and updating the random seed.
*/
void tStorm::GenerateStorm()
{
   if( optVariable )
   {
      p = pMean*ExpDev( &seed );
      stdur = stdurMean*ExpDev( &seed );
      istdur = istdurMean*ExpDev( &seed );
      srand( seed );
   }
}


/*
**  tStorm::ExpDev:  Finds a random number with an exponential distribution
**                   (adapted from Numerical Recipes).
*/
float tStorm::ExpDev( long *idum )
{
    float dum;

    do
        dum = ran3( idum );
    while( dum==0.0 );
    return -log(dum);
}


/*
**  GetStormDuration
**
**  Returns the storm duration.
*/
float tStorm::GetStormDuration()
{
   return stdur;
}

/*
**  InterstormDur
**
**  Returns the interstorm duration.
*/
float tStorm::InterstormDur()
{
   return istdur;
}

/*
**  GetRainrate
**
**  Returns the rainfall rate.
*/
float tStorm::GetRainrate()
{
   return p;
}


/*
**  GammaDev
**
**  Returns a random variable drawn from a Gamma distribution with parameter m.
**
**  (Note: not actually called; provided for future use).
*/
double tStorm::GammaDev(double m, long * idum)
{
  double x, y,z, c,t,b,u,w,v;
  if (m<1)
    {
      c = 1/m;
      t = 0.07 + 0.75*sqrt(1-m);
      b = 1 + exp(-t)*m/t;
      int accept = 0;
      while (accept == 0)
        {
          u = ran3(idum);
          w = ran3(idum);
          v = b *u;
          if (v<=1)
            {
              x  = t * pow(v, c);
              accept = ((w<=((2-x)/(2+x))) || (w<=exp(-x)));
            }
          else
            {
              x = -log(c*t*(b-v));
              y = x/t;
              accept = (((w*(m + y - m*y)) <= 1) || (w<= ( pow(y, (m-1)))));
            }
        }
    }
  else
    {
      b = m-1;
      c = 3*m - 0.75;
      int accept = 0;
      while (accept == 0)
        {
          u = ran3(idum); v = ran3(idum);
          w = u* ( 1-u);
          y = sqrt(c/w) * (u - 0.5);
          x = b + y;
          if ( x>= 0)
            {
              z = 64*( pow(w,3))*v*v;
              accept = (z <= ( 1 - 2*y*y/x)) || ( log(z) <= (2*(b*log(x/b) - y)));
            }
        }
    }
  return x;
}


/*
**  ran3:  Random number generator from Numerical Recipes.
*/
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



