/**************************************************************************/
/**
**  @file tStorm.cpp
**  @brief Functions for class tStorm.
**
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm. It is essentially an implementation of the model of
**  P. Eagleson, 1978b, Water Resources Research. Its services include
**  reading the necessary parameters from a tInputFile, generating a new      
**  storm, and reporting its various values.
**
**  $Id: tStorm.cpp,v 1.31 2003-10-15 09:20:09 childcvs Exp $
*/
/**************************************************************************/


#include <math.h>
#include <string.h>
#include "../Mathutil/mathutil.h"
#include "tStorm.h"
#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif



/**************************************************************************\
**
**  tStorm::tStorm:  Constructor for storms. The default constructor
**                   assigns a value of unity to storm depth, duration,
**                   and interstorm duration. (Note:
**                   this constructor does not allow option for sinusoidal
**                   variation in means).
**
\**************************************************************************/
tStorm::tStorm( bool optvar )
  :
  rand(0),
  stdurMean(1.0),
  istdurMean(1.0),
  pMean(1.0),
  p(1.0),
  stdur(1.0),
  istdur(1.0),
  p0(-1.),
  stdur0(-1.),
  istdur0(-1.),
  pdev(-1.),
  stdurdev(-1.),
  istdurdev(-1.),
  twoPiLam(-1.),
  endtm(1.0e9),
  optVariable(optvar),
  optSinVar(false)
{
   //srand( 0 );
}


/**************************************************************************\
**
**  tStorm::tStorm:  The constructor that's really used assigns values for
**                   mean storm depth, duration, and interstorm duration,
**                   initializes current depth, etc, to the mean values,
**                   and initializes the random number generator. (Note:
**                   this constructor does not allow option for sinusoidal
**                   variation in means).
**
\**************************************************************************/
tStorm::tStorm( double mp, double ms, double mis, tRand *rand_, bool optvar, double et )
  :
  rand(rand_),
  stdurMean(ms),
  istdurMean(mis),
  pMean(mp),
  p(pMean),
  stdur(stdurMean),
  istdur(istdurMean),
  p0(-1.),
  stdur0(-1.),
  istdur0(-1.),
  pdev(-1.),
  stdurdev(-1.),
  istdurdev(-1.),
  twoPiLam(-1.),
  endtm(et),
  optVariable(optvar),
  optSinVar(false)
{
}


/**************************************************************************\
**
**  tStorm::tStorm
**
**  Alternative constructor that reads parameters directly from a tInputFile
**  object. Reads option for variable storms (normally this is "yes"---that's
**  the point of these objects---but a user may wish to switch off variation
**  as a test), mean values for rainfall intensity, duration, and interstorm
**  period, and a random seed to initialize the random number generator.
**    Also reads an option for long-term sinusoidal variations in the mean
**  values, and if the option is selected, reads the relevant parameters.
**  Variables p0, stdur0, and istdur0 are the mean values of the means;
**  pdev, stdurdev, and istdurdev are the range of variation (e.g., if pMean
**  were to fluctuate between 5 and 10, p0 would be 7.5 and pdev 2.5).
**
**  Modifications:
**   - 3/00 initialization now includes creation of ".storm" file for storm
**     history (GT)
**
\**************************************************************************/
tStorm::tStorm( const tInputFile &infile, tRand *rand_ ) :
  rand(rand_)
{
   // Read + set parameters for storm intensity, duration, and spacing
   {
     int tmp_;
     tmp_ = infile.ReadItem( tmp_, "OPTVAR" );
     optVariable = BOOL(tmp_ != 0);
   }
   pMean = infile.ReadItem( pMean, "PMEAN" );
   stdurMean = infile.ReadItem( stdurMean, "STDUR" );
   istdurMean = infile.ReadItem( istdurMean, "ISTDUR" );
   p = pMean;
   stdur = stdurMean;
   istdur = istdurMean;

 
   endtm = infile.ReadItem( endtm, "RUNTIME" );
   double help;

   help = infile.ReadItem( help, "OPTREADINPUT" );
   if(help>0){
      help = infile.ReadItem( help, "INPUTTIME" );
      endtm += help;
   }


   // Handle option for sinuidoil variation in means
   {
     int tmp_;
     tmp_ = infile.ReadItem( tmp_, "OPTSINVAR" );
     optSinVar = BOOL(tmp_ != 0);
   }
   if( optSinVar )
   {
      p0 = pMean;
      stdur0 = stdurMean;
      istdur0 = istdurMean;
      twoPiLam = (2.0*PI)/(infile.ReadItem( twoPiLam, "PERIOD" ));
      pdev = infile.ReadItem( pdev, "MAXPMEAN" ) - pMean;
      if( pdev<0 ) cerr << "Warning: MAXPMEAN < PMEAN !";
      else if( pdev > pMean ) cerr << "Warning: MINPMEAN < 0 !";
      stdurdev = infile.ReadItem( stdurdev, "MAXSTDURMN" ) - stdurMean;
      if( stdurdev<0 ) cerr << "Warning: MAXSTDURMN < STDURMN !";
      else if( stdurdev > stdurMean ) cerr << "Warning: MINSTDURMN < 0 !";
      istdurdev = infile.ReadItem( istdurdev, "MAXISTDURMN" ) - istdurMean;
      if( istdurdev<0 ) cerr << "Warning: MAXISTDURMN < ISTDURMN !";
      else if( istdurdev > istdurMean ) cerr << "Warning: MINISTDURMN < 0 !";
   }

   // If variable storms used, create a file for writing them
   if( optVariable )
   {
      char fname[87];
#define THEEXT ".storm"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      stormfile.open( fname );
      if( !stormfile.good() )
          cerr << "Warning: unable to create storm data file '" 
               << fname << "'\n";
   }
}


/**************************************************************************\
**
**  tStorm::GenerateStorm
**
**  Generates a new storm by drawing new values of p, stdur, and istdur from
**  an exponential distribution and updating the random seed.
**    If the minp parameter is greater than zero, the function will keep
**  picking storms until it finds one with p>minp. The total elapsed time,
**  including the rejected storms and their associated interstorm periods,
**  is stored istdur.
**
**  Inputs:      minp -- minimum value of rainfall rate p (default 0)
**               mind -- minimum storm depth to produce runoff (default 0)
**               tm -- current time in simulation
**  Members updated:  p, stdur, istdur take on new random values (if optVar)
**                    pMean, stdurMean, istdurMean adjusted (if optSinVar)
**  Assumptions:  pMean > 0
**
**  Modifications:
**   - changed AND to OR in while loop, GT 5/99
**   - added to while loop an additional check to see if time has run out,
**     NG 2/00
**   - added output of time, storm intensity, & duration to ".storm" file
**     GT 3/00
**
\**************************************************************************/
void tStorm::GenerateStorm( double tm, double minp, double mind )
{
   // If option for sinusoidal variation is on, adjust the means
   // Also set values for current storm to the means, in case option for
   // random storms is off.
   if( optSinVar )
   {
      double sinfn = sin( tm*twoPiLam );
      pMean = p0 + pdev*sinfn;
      //cout << "pMean = " << pMean << endl;
      stdurMean = stdur0 + stdurdev*sinfn;
      istdurMean = istdur0 + istdurdev*sinfn;
      p = pMean;
      stdur = stdurMean;
      istdur = istdurMean;
   }
   
   // If option for random storms is on, pick a storm at random.
   // Keep picking and accumulating time until the storm depth or intensity
   // is greater than the minimum value needed to produce runoff.

   if( optVariable )
   {
      stdur = 0;
      istdur = 0;
      do
      {
         p = pMean*ExpDev();
         istdur += istdurMean*ExpDev() + stdur;
         stdur = stdurMean*ExpDev();
	 if(0) { // Debug
	   cout << "P " << p << "  ST " << stdur << "  IST " << istdur
		<< "  DP " << p*stdur << " minp " << minp << " mind " <<mind <<
	     endl;
	 }
      } while( (p<=minp || (p*stdur)<=mind) && (tm+istdur+stdur<endtm) );
      stormfile << istdur << " " << p << " " << stdur << endl;
   }
}


/**************************************************************************\
**
**  tStorm::ExpDev:  Finds a random number with an exponential distribution
**                   (adapted from Numerical Recipes).
**
\**************************************************************************/
double tStorm::ExpDev() const
{
    double dum;

    do
        dum = rand->ran3();
    while( dum==0.0 );
    return -log(dum);
}


/**************************************************************************\
**
**  tStorm "get" routines: return various variables
**
\**************************************************************************/

double tStorm::getStormDuration() const { return stdur; }
double tStorm::interstormDur() const { return istdur; }
double tStorm::getRainrate() const { return p; }
double tStorm::getMeanStormDur() const {return stdurMean;}
double tStorm::getMeanInterstormDur() const {return istdurMean;}
double tStorm::getMeanPrecip() const {return pMean;}

/**************************************************************************\
**
**  GammaDev
**
**  Returns a random variable drawn from a Gamma distribution with parameter m.
**
**  (Note: not actually called; provided for future use).
**
\**************************************************************************/
double tStorm::GammaDev(double m) const
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
          u = rand->ran3();
          w = rand->ran3();
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
          u = rand->ran3(); v = rand->ran3();
          w = u* ( 1-u);
          y = sqrt(c/w) * (u - 0.5);
          x = b + y;
          if ( x>= 0)
            {
              z = 64*( pow(w,3.))*v*v;
              accept = (z <= ( 1 - 2*y*y/x)) || ( log(z) <= (2*(b*log(x/b) - y)));
            }
        }
    }
  return x;
}

