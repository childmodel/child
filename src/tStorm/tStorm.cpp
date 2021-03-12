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
**  $Id: tStorm.cpp,v 1.36 2004-06-16 13:37:42 childcvs Exp $
*/
/**************************************************************************/


#include <math.h>
#include <string.h>
#include "../Mathutil/mathutil.h"
#include <iostream>
#include <fstream>

#include "tStorm.h"

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
  p(1.0),
  stdur(1.0),
  istdur(1.0),
  endtm(1.0e9),
  optVariable(optvar)
{
   //srand( 0 );
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
tStorm::tStorm( const tInputFile &infile, tRand *rand_, 
		bool no_write_mode /* = false */ ) :
  rand(rand_)
{
   // Read + set parameters for storm intensity, duration, and spacing
   optVariable = infile.ReadBool( "OPTVAR" );
   if( !no_write_mode )
     {
       infile.WarnObsoleteKeyword("PMEAN", "ST_PMEAN");
       infile.WarnObsoleteKeyword("STDUR", "ST_STDUR");
       infile.WarnObsoleteKeyword("ISTDUR", "ST_ISTDUR");
       infile.WarnObsoleteKeyword("OPTSINVAR", "ST_PMEAN");
       infile.WarnObsoleteKeyword("PERIOD", "ST_PMEAN");
       infile.WarnObsoleteKeyword("START_CYCLE_TIME", "ST_PMEAN");
       infile.WarnObsoleteKeyword("MAXPMEAN", "ST_PMEAN");
       infile.WarnObsoleteKeyword("MAXSTDURMN", "ST_STDUR");
       infile.WarnObsoleteKeyword("MAXISTDURMN", "ST_ISTDUR");
     }
   infile.ReadItem( p_ts, "ST_PMEAN");
   infile.ReadItem( stdur_ts, "ST_STDUR");
   infile.ReadItem( istdur_ts, "ST_ISTDUR");

   p = p_ts.calc(0.);
   stdur = stdur_ts.calc(0.);
   istdur = istdur_ts.calc(0.);

   endtm = infile.ReadItem( endtm, "RUNTIME" );
   double help;

   help = infile.ReadItem( help, "OPTREADINPUT" );
   if(help>0){
      help = infile.ReadItem( help, "INPUTTIME" );
      endtm += help;
   }

   // If variable storms used, create a file for writing them
   if( optVariable && !no_write_mode )
   {
      char fname[87];
#define THEEXT ".storm"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      stormfile.open( fname );
      if( !stormfile.good() )
          std::cerr << "Warning: unable to create storm data file '"
		    << fname << "'\n";
   }
}

tStorm::tStorm( const tStorm& orig )
   :  stormfile(),
      rand(orig.rand),
      p_ts(orig.p_ts),
      stdur_ts(orig.stdur_ts),
      istdur_ts(orig.istdur_ts),
      p(orig.p),
      stdur(orig.stdur),
      istdur(orig.istdur),
      endtm(orig.endtm),
      optVariable(orig.optVariable)
{}
/**************************************************************************\
**
**  tStorm::TurnOnOutput, TurnOffOutput
**
**  Open output file so output will be directed to it, 
**  or close output file so there won't be output.
**
**  SL, 10/2010
**
\**************************************************************************/
void tStorm::TurnOnOutput( const tInputFile& infile )
{
     // If variable storms used, create a file for writing them
  if( !stormfile.good() && optVariable )
   {
      char fname[87];
#define THEEXT ".storm"
      infile.ReadItem( fname, sizeof(fname)-sizeof(THEEXT), "OUTFILENAME" );
      strcat( fname, THEEXT );
#undef THEEXT
      stormfile.open( fname );
      if( !stormfile.good() )
          std::cerr << "Warning: unable to create storm data file '"
		    << fname << "'\n";
   }
}

void tStorm::TurnOffOutput()
{
  if( stormfile.good() )
    stormfile.close();
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

   p = p_ts.calc(tm);
   stdur = stdur_ts.calc(tm);
   istdur = istdur_ts.calc(tm);

   // If option for random storms is on, pick a storm at random.
   // Keep picking and accumulating time until the storm depth or intensity
   // is greater than the minimum value needed to produce runoff.

   if( optVariable )
   {
      const double pMean = p;
      const double stdurMean = stdur;
      const double istdurMean = istdur;

      stdur = 0;
      istdur = 0;
      do
      {
         p = pMean*ExpDev();
         istdur += istdurMean*ExpDev() + stdur;
         stdur = stdurMean*ExpDev();
	 if(0) { // Debug
	   std::cout << "P " << p << "  ST " << stdur << "  IST " << istdur
		     << "  DP " << p*stdur << " minp " << minp
		     << " mind " <<mind << std::endl;
	 }
      } while( (p<=minp || (p*stdur)<=mind) && (tm+istdur+stdur<endtm) );
      if( stormfile.good() )
	stormfile << istdur << " " << p << " " << stdur << std::endl;
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

/**************************************************************************\
**
**  tStorm "set" routines: set various variables
**
\**************************************************************************/
void tStorm::setRainrate( double pMeanNew ) 
{
  std::stringstream ss;
  ss << pMeanNew;
  p_ts.reconfigure( ss.str().c_str() );
}


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
      bool accept = false;
      while (!accept)
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
      bool accept = false;
      while (!accept)
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

