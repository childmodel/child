//-*-c++-*- 

/**************************************************************************/
/**
**  @file tEolian.h
**  @brief header file for class tEolian
**
**  Loess! Class tEolian is a simple little "process class" that simulates
**  deposition of wind-borne sediment across the landscape.
**    In this initial version, the deposition rate is uniform in space
**  and time. It is implemented by simply depositing a given depth of
**  material at each active node during each time step (the depth is
**  equal to the rate times time-step duration). If multiple grain sizes
**  are used, deposited material consists of 100% of the finest
**  size-fraction, which is assumed to the first size-fraction.
**    Capabilities: a constructor that sets up parameters, and a
**  DepositLoess() function that does the deposition.
**    Note that a more advanced version could use a polynomial surface
**  to represent deposition rates: rate(x,y) = x*(ax*x+bx)+y*(ay*y+by)+c.
**
**  (Created 2/99 by GT)
**
**  $Id: tEolian.h,v 1.10 2004/04/16 18:35:47 childcvs Exp $
*/
/**************************************************************************/

#ifndef TEOLIAN_H
#define TEOLIAN_H

#include <assert.h>
#include "../tMesh/tMesh.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"


/**************************************************************************/
/**
**  @class tEolian
**
**  Class tEolian contains two parameters: the loess deposition rate,
**  and an array of deposition depths for each grain size fraction.
**  All but the first (assumed to be the finest) size-fractions are
**  initialized to zero and never changed.
**
*/
/**************************************************************************/
class tEolian
{
    tEolian();
public:
    tEolian( const tInputFile &infile );
    void DepositLoess( tMesh<tLNode> *mp, double delt, double ctime );

private:
    double loessDepRate;
    tArray<double> depositDepth;
};

#endif
