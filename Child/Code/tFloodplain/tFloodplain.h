/**************************************************************************\
**
**  tFloodplain.h: header file for class tFloodplain and its helper class
**                 tFloodNode.
**
**  Class tFloodplain simulates overbank deposition in the CHILD model,
**  using a modified version of Howard's (1996) diffusive overbank
**  deposition model. Implementation of the model is described in
**  tFloodplain.cpp. The class consists of parameters for the model,
**  a pointer to the mesh, a constructor to read in parameters, and
**  a DepositOverbank function to implement deposition during a given
**  storm event. Helper class tFloodNode stores info for a "flood node"
**  (one containing a channel big enough to generate a sedimentologically
**  significant flood), and is used to construct a list of flood nodes.
**
**  (Created 1/99 by GT)
**
**  $Id: tFloodplain.h,v 1.4 1999-04-05 14:54:38 gtucker Exp $
\**************************************************************************/

#ifndef TFLOODPLAIN_H
#define TFLOODPLAIN_H

#include <assert.h>
#include "../tMesh/tMesh.h"
#include "../tMeshList/tMeshList.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"

#define kVeryFar 1.0e12


/**************************************************************************\
**** class tFloodplain ****************************************************
**
**  Class tFloodplain contains the parameters for the modified Howard
**  floodplain deposition model, as well as a constructor to read in
**  parameters from an input file and a DepositOverbank function to
**  implement the model.
**
**  Modifications:
**    - 3/13/99 GT added deparr variable which is used by DepositOverbank
**              (could be a local var but this saves having to reinit
**              it with every call to tFloodplain)
**
\**************************************************************************/
class tFloodplain
{

public:
    tFloodplain( tInputFile &infile, tMesh<tLNode> *mp );
    void DepositOverbank( double precip, double delt, double ctime );
    
private:
    double fpmu;            // "mu" parameter of Howard model
    double fplamda;         // "lamda" (distance-decay) parameter
    double kdb;             // depth-disch coeff (lumped; see tFloodplain.cpp)
    double event_min;       // bankfull event precip rate
    double drarea_min;      // min drainage area for a "flood node"
    double mqs;             // depth-disch at-a-station exponent
    double mqbmqs;          // bankfull minus at-a-station exponents
    tArray<double> deparr;  // depth deposited (# grn size; all but 1st=0)
    tMesh<tLNode> *meshPtr; // ptr to mesh
};


/**************************************************************************\
**** class tFloodNode *****************************************************
**
**  Basically just a data record used to construct a list of flood nodes.
**  Contains a pointer to the node and the water surface height.
**
\**************************************************************************/
class tFloodNode
{
    friend class tFloodplain;
    
private:
    tLNode *nodePtr; // ptr to flood node
    double wsh;      // water surface height at flood node
};

#endif
