//-*-c++-*- 

/************************************************************************/
/**
**  @file tUplift.h
**  @brief Header file for tUplift objects.
**
**  tUplift objects manage data and functions related to all types of
**  deformation and baselevel (including sea level) change. Thus, the
**  name "tUplift" is a bit misleading in that the relative motions that
**  the tUplift class models need not be solely vertical, but could for
**  example include strike-slip deformation. Moreover, the class
**  models _relative_ motion, which of course need not be of tectonic
**  origin but could instead reflect sealevel or other baselevel changes.
**  Since we have no good term that encompasses all these possibilities,
**  "tUplift" seems the most convenient catch-all.
**
**  The constructor reads in a code for the type of boundary
**  condition desired, and then reads in parameters appropriate to that
**  type. The "uplift" is then applied via a general DoUplift()
**  function, which in turn calls the appropriate function (e.g.,
**  UpliftUniform()) to implement the desired behavior.
**
**  Major modifications:
**    - added StrikeSlip and FoldPropErf functions (gt, May 2000)
**    - added FaultBendFold function (srm, August 2002)
**    - added UpliftRateMap functions (GT, June 2006)
**    - added PropagatingFront function (GT & BY, Julu 2007)
**    - added time series rate variable rate_ts, and implemented it
**      for two uplift functions (which now take current time as a
**      parameter)
**
**  $Id: tUplift.h,v 1.26 2008/07/09 16:35:34 childcvs Exp $
*/
/************************************************************************/

#ifndef TUPLIFT_H
#define TUPLIFT_H

#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"
#include "../tTimeSeries/tTimeSeries.h"

class tUplift
{
public:
    tUplift( const tInputFile &infile );
    void DoUplift( tMesh<tLNode> *mp, double delt, double current_time );
    double getDuration() const;
    double getRate() const;
private:
   void UpliftUniform( tMesh<tLNode> *mp, double delt, double currentTime );
   void BlockUplift( tMesh<tLNode> *mp, double delt, double currentTime );
   void StrikeSlip( tMesh<tLNode> *mp, double delt ) const;
   void FoldPropErf( tMesh<tLNode> *mp, double delt );
   void CosineWarp2D( tMesh<tLNode> *mp, double delt );
   void PropagatingFold( tMesh<tLNode> *mp, double delt ) const;
   void TwoSideDifferential( tMesh<tLNode> *mp, double delt ) const;
   void FaultBendFold( tMesh<tLNode> *mp, double delt ) const;
   void FaultBendFold2( tMesh<tLNode> *mp, double delt ) const;
   void NormalFaultTiltAccel( tMesh<tLNode> *mp, double delt, double currentTime ) const;
   void LinearUplift( tMesh<tLNode> *mp, double delt );
   void PowerLawUplift( tMesh<tLNode> *mp, double delt );
   void UpliftRateMap( tMesh<tLNode> *mp, double delt, double currentTime );
   void PropagatingFront( tMesh<tLNode> *mp, double delt, double currentTime );
   
private:
   typedef enum {
       kNoUplift = 0,
       k1,
       k2,
       k3,
       k4,
       k5,
       k6,
       k7,
       k8,
       k9,
       k10,
       k11,
       k12,
	   k13
   } tUplift_t;
   
   static tUplift_t DecodeType(int);
   
   tUplift_t typeCode;    // Code for the type of uplift desired
   double duration;       // Duration of uplift
   double rate;           // Rate of uplift
   double rate2;          // Second rate (e.g., second structure)
   tTimeSeries rate_ts;   // Rate of uplift as a time series
   double faultPosition;  // Position of fault (y-location)
   double positionParam1; // Another position parameter
   double slipRate;       // Slip rate for strike-slip motion and fault prop
   double foldParam;      // Parameter used in folding calculation
   double foldParam2;     // Another one
   double deformStartTime1; // Parameter for onset of uplift/deformation
   double flatDepth;      // Depth below surface (at faultPos) of thrust flat
   double rampDip;	   // Dip of thrust ramp (in degrees)
   double kinkDip;	   // Dip of axial surface in hangingwall that
   //	initiates at lower end of ramp.
   double upperKinkDip;   // Dip of axial surface that initiates at upper
   //	end of ramp.
   double meanElevation;  // Mean elevation of surface at t0.
   double decayParam;  //decay rate for powerlaw
   double timeParam1;     // Timing parameter
   double width;  //y-size of mesh, needed for setting up linear uplift change
   double dupdy;  //change in uplift with y for linear uplift change
   int optincrease; //for linear uplift change - 0=>dec to mtn frnt, 1=>inc to mtn frnt
   char mUpliftMapFilename[120]; // Name for files containing uplift maps
   int miNumUpliftMaps;     // Number of uplift maps
   tArray<double> mUpliftMapTimes; // Time corresponding to each uplift map
   int miCurUpliftMapNum;  // Current uplift map number
   double mdNextUpliftMapTime; // Time at which to read next map
   double mdUpliftFrontGradient; // Horizontal gradient (dy/dx) of propagating uplift front
   
private:
   tUplift();
};

#endif








