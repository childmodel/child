/**************************************************************************\
**
**  tStreamNet.h: Header file for class tStreamNet.
**
**  tStreamNet objects contain data and functions related to flow routing
**  and sediment transport across the landscape surface.
**
**  $Id: tStreamNet.h,v 1.20 1998-05-05 19:46:17 gtucker Exp $
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

#define kSaturatedFlow 1   // Option for saturation-excess flow generation
#define kConstSoilStore 2  // Option for "bucket"-type flow generation

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
   void SetVoronoiVertices();
   void MakeFlow();
   void FlowUniform();
   void FlowSaturated();
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
