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
**
**  $Id: tUplift.h,v 1.19 2004-04-16 18:32:54 childcvs Exp $
*/
/************************************************************************/

#ifndef TUPLIFT_H
#define TUPLIFT_H

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"


class tUplift
{
public:
    tUplift( const tInputFile &infile );
    void DoUplift( tMesh<tLNode> *mp, double delt );
    double getDuration() const;
    double getRate() const;
private:
    void UpliftUniform( tMesh<tLNode> *mp, double delt );
    void BlockUplift( tMesh<tLNode> *mp, double delt );
    void StrikeSlip( tMesh<tLNode> *mp, double delt ) const;
    void FoldPropErf( tMesh<tLNode> *mp, double delt );
    void CosineWarp2D( tMesh<tLNode> *mp, double delt );
    void PropagatingFold( tMesh<tLNode> *mp, double delt ) const;
    void TwoSideDifferential( tMesh<tLNode> *mp, double delt ) const;
    void FaultBendFold( tMesh<tLNode> *mp, double delt ) const;
    void FaultBendFold2( tMesh<tLNode> *mp, double delt ) const;

private:
    int typeCode;          // Code for the type of uplift desired
    double duration;       // Duration of uplift
    double rate;           // Rate of uplift
    double rate2;          // Second rate (e.g., second structure)
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

private:
    tUplift();
};

#endif








