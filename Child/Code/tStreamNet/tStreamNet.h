/**************************************************************************\
**
**  tStreamNet.h: Header file for class tStreamNet and related class
**                tInlet.
**
**  The tStreamNet class contains data and functions related to flow 
**  routing across the landscape surface. Data include hydrologic
**  parameters such as infiltration capacity, soil transmissivity, etc,
**  as well as flags that indicate which of several different options
**  the user desires (e.g., type of runoff generation, whether uniform
**  infiltration-excess, TOPMODEL-type saturation excess, etc). Actions
**  include setting drainage directions, computing total drainage areas,
**  computing overland flow discharge, and filling "lakes".
**
**  Class tInlet is used to model the entry of a river at an edge of the
**  model mesh.
**
**  $Id: tStreamNet.h,v 1.22 1998-06-04 21:28:10 gtucker Exp $
\**************************************************************************/

#ifndef TSTREAMNET_H
#define TSTREAMNET_H

#include "../Definitions.h"
#include "../Classes.h"
#include "../tUplift/tUplift.h"
#include "../Erosion/erosion.h"
#include "../GridElements/gridElements.h"
#include "../tGrid/tGrid.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"
#include "../tStorm/tStorm.h"

#define kHortonian 0       // Option for uniform infilt-excess runoff
#define kSaturatedFlow1 1  // Option for sat-excess runoff w/ return flow
#define kSaturatedFlow2 2  // Option for sat-excess runoff w/o return flow
#define kConstSoilStore 3  // Option for "bucket"-type flow generation
#define kSecperyear 31536000  // No. of seconds in one year


/** class tInlet *************************************************************/
class tInlet
{
    friend class tStreamNet;
public:
    tInlet();
    tInlet( tGrid< tLNode > *, tInputFile & );
    ~tInlet();
   void FindNewInlet();
   double getInSedLoad() const;
   double getInSedLoad( int );
   tArray< double > getInSedLoadm() const;
   void setInSedLoad( double );
   void setInSedLoad( int, double );
   double getInDrArea() const;
   void setInDrArea( double );
   tLNode *getInNodePtr();
   void setInNodePtr( tLNode * );
private:
    tLNode *innode;
    double inDrArea;
   double inSedLoad;
   tArray< double > inSedLoadm; // incoming sediment load if multi-sizes
    tGrid< tLNode > *gridPtr;
};



/** class tStreamNet *********************************************************/
class tStreamNet
{
   friend class tStreamTransport;
   friend class tStreamMeander;
public:
   tStreamNet();
   tStreamNet( tGrid< tLNode > &, tStorm &, tInputFile & );
   ~tStreamNet();
   void ResetGrid( tGrid< tLNode > & );
   const tGrid< tLNode > *getGridPtr() const;
   tGrid< tLNode > *getGridPtrNC();
   const tStorm *getStormPtr() const;
   tStorm *getStormPtrNC();
   int getFlowGenOpt() const;
   int getFillLakesOpt() const;
   double getRainRate() const;
   double getTransmissivity() const;
   double getInfilt() const;
   double getSoilStore() const;
   double getInDrArea() const;
   double getInSedLoad() const;
   tArray< double > getInSedLoadm() const;
   tLNode *getInletNodePtr() const;
   tLNode *getInletNodePtrNC();
   double getMndrDirChngProb() const;
   void setFlowGenOpt( int );
   void setFillLakesOpt( int );
   void setRainRate( double );
   void setTransmissivity( double );
   void setInfilt( double );
   void setInDrArea( double );
   void setInSedLoad( double );
   void setInSedLoadm( int, double );
   void setInletNodePtr( tLNode * );
   void setMndrDirChngProb( double );
   void UpdateNet();
   void UpdateNet( tStorm & );
   void CheckNetConsistency();
   void CalcSlopes();
   void InitFlowDirs();
   void FlowDirs();
   void DrainAreaVoronoi();
   void RouteFlowArea( tLNode *, double );
   void RouteRunoff( tLNode *, double, double );
   void SetVoronoiVertices();
   void MakeFlow();
   void FlowUniform();
   void FlowSaturated1();
   void FlowSaturated2();
   void FlowBucket();
   void FillLakes();
   int FindLakeNodeOutlet( tLNode * );
   void SortNodesByNetOrder();
   int DamBypass( tLNode * );
   
protected:
   tGrid< tLNode > * gridPtr;
   tStorm *stormPtr;
   int flowgen;
   int filllakes;
   double rainrate;
   double trans;
   double infilt;
   double soilStore;        /* soil water storage, depth equiv (m) */
   //double inDrArea;
   tInlet inlet;
   double mndrDirChngProb;
};

#endif
