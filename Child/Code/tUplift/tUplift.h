//-*-c++-*- 

/************************************************************************\
**
**  tUplift.h:  Header file for tUplift objects.
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
**
**  $Id: tUplift.h,v 1.13 2002-09-23 12:10:44 arnaud Exp $
\************************************************************************/

#ifndef TUPLIFT_H
#define TUPLIFT_H

#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"


class tUplift
{
public:
    tUplift( tInputFile &infile );
    void DoUplift( tMesh<tLNode> *mp, double delt );
    double getDuration();
    double getRate() const;
private:
    void UpliftUniform( tMesh<tLNode> *mp, double delt );
    void BlockUplift( tMesh<tLNode> *mp, double delt );
    void StrikeSlip( tMesh<tLNode> *mp, double delt );
    void FoldPropErf( tMesh<tLNode> *mp, double delt );
    void CosineWarp2D( tMesh<tLNode> *mp, double delt );
    void PropagatingFold( tMesh<tLNode> *mp, double delt );
    void TwoSideDifferential( tMesh<tLNode> *mp, double delt );

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

};

#endif








