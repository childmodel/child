/**************************************************************************\
**
**  tFloodplain.h: header file for class tFloodplain and its helper class
**                 tFloodNode.
**
**  Class tFloodplain simulates overbank deposition in the CHILD model,
**  using a modified version of Howard's (1996) diffusive overbank
**  deposition model.
**
**  (Created 1/99 by GT)
**
**  $Id: tFloodplain.h,v 1.1 1999-01-25 16:36:22 gtucker Exp $
\**************************************************************************/

#ifndef TFLOODPLAIN_H
#define TFLOODPLAIN_H

#include <assert.h>
#include "../tGrid/tGrid.h"
#include "../tGridList/tGridList.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"

#define kVeryFar 1.0e12


/**************************************************************************\
**** class tFloodplain ****************************************************
**
**
\**************************************************************************/
class tFloodplain
{

public:
    tFloodplain( tInputFile &infile, tGrid<tLNode> *gp );
    void DepositOverbank( double precip, double delt, double ctime );
    
private:
    double fpmu;
    double fplamda;
    double kdb;
    double event_min;
    double drarea_min;
    double mqs;
    double mqbmqs;
    tGrid<tLNode> *gridPtr;
};


/**************************************************************************\
**** class tFloodNode *****************************************************
**
**
\**************************************************************************/
class tFloodNode
{
    friend class tFloodplain;
    
private:
    tLNode *nodePtr;
    double wsh;
};

#endif
