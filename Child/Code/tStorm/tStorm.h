/*
**  tStorm.h
** 
**  Header for tStorm objects.
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
**
**  Version 1.0, Greg Tucker, November 1997.
**  $Id: tStorm.h,v 1.4 1998-01-27 23:34:12 stlancas Exp $
*/

#ifndef TSTORM_H
#define TSTORM_H

/*#include "../tRunTimer/tRunTimer.h"*/
#include "../tInputFile/tInputFile.h"


class tStorm
{
  public:
   tStorm( int optVariable=1 );
   tStorm( float, float, float, unsigned, int optvar=1 );
   tStorm( tInputFile & );
   void  GenerateStorm();
   float GetStormDuration();
   float InterstormDur();
   float GetRainrate();
   float getMeanStormDur() const;
   float getMeanInterstormDur() const;
   float getMeanPrecip() const;
   
  private:
   float ExpDev( long * );
   double GammaDev(double, long*);
   
   int optVariable;   // Flag indicating whether storms are random or not
   float stdurMean;   // Mean duration
   float istdurMean;  // Mean time between storms
   float pMean;       // Mean rainfall intensity
   float p;           // Actual rainfall intensity for the current storm
   float stdur;       // Actual storm duration
   float istdur;      // Actual time between storms
   long  seed;        // Random seed
};


#endif
