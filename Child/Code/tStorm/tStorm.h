/**************************************************************************\
*
**  tStorm.h: Header for class tStorm
**
**  A tStorm object generates random storms assuming an exponential
**  distribution of rainfall intensity, storm duration, and time to the
**  next storm. It is essentially an implementation of the model of
**  P. Eagleson, 1978b, Water Resources Research. Its services include
**  reading the necessary parameters from a tInputFile, generating a new      
**  storm, and reporting its various values.
**    If you want to provide an option for NOT having storms vary
**  randomly, you can do so by setting optVariable to zero on initialization.
**    The GammaDev() function is provided for future reference; it is not
**  actually used in the current version.
**    At the user's option, the storm parameters can also be varied 
**  sinusoidally with time to simulate long-term climatic fluctuations.
**  (This is done by GenerateStorm).
**    tStorm objects could be easily modified (or inherited from) to use
**  different distributions. They can also be modified to create objects
**  for other random processes such as river flows, etc.
**
**  Created by Greg Tucker, November 1997.
**
**  $Id: tStorm.h,v 1.13 1999-04-01 17:31:11 gtucker Exp $
\**************************************************************************/

#ifndef TSTORM_H
#define TSTORM_H

#include <stdlib.h>
#include "../tInputFile/tInputFile.h"

#define kSecperyear 31536000

class tStorm
{
public:
    tStorm( int optVariable=1 );
    tStorm( double, double, double, unsigned, int optvar=1 );
    tStorm( tInputFile & );
    void  GenerateStorm( double tm, double minp=0.0, double mind=0.0 );
    double getStormDuration();
    double interstormDur();
    double getRainrate();
    double getMeanStormDur() const;
    double getMeanInterstormDur() const;
    double getMeanPrecip() const;
    double getOptVar() const;
   
private:
    double ExpDev( long * );
    double GammaDev(double, long*);
   
    int optVariable;   // Flag indicating whether storms are random or not
    int optSinVar;     // Option for sinusoidal variation in storm params
    double stdurMean;   // Mean duration
    double istdurMean;  // Mean time between storms
    double pMean;       // Mean rainfall intensity
    double p;           // Actual rainfall intensity for the current storm
    double stdur;       // Actual storm duration
    double istdur;      // Actual time between storms
    double p0;         // Climatological mean: the "weather" means can 
    double stdur0;     //  themselves vary over geologic time; p0, etc, are
    double istdur0;    //  the means of the means so to speak.
    double pdev;       // Absolute magnitude of deviation from p0, etc, under
    double stdurdev;   //  sinusoidal variation.
    double istdurdev;
    double twoPiLam;   // Parameter for sinusoidal variation: 2pi / period
    long  seed;        // Random seed
};


inline double tStorm::getOptVar() const {return optVariable;}

#endif
