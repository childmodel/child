/**************************************************************************\
**
**  tStreamNet.h: Header file for class tStreamNet.
**
**  tStreamNet objects contain data and functions related to flow routing
**  and sediment transport across the landscape surface.
**
**  $Id: tStreamNet.h,v 1.16 1998-04-22 20:34:44 nmgaspar Exp $
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
   tArray< double > getInSedLoadm() const;
   void setInSedLoad( double );
   void setInSedLoadm( int, double );
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
   double getInDrArea() const;
   double getInSedLoad() const;
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
   void setInletNodePtr( tLNode * );
   void setMndrDirChngProb( double );
   void UpdateNet();
   void UpdateNet( tStorm & );
   void CheckNetConsistency();
   void CalcSlopes();
   //Xvoid CalcVAreas();
   void InitFlowDirs();
   void FlowDirs();
   void DrainAreaVoronoi();
   //Xdouble VoronoiArea( tLNode * );
   void RouteFlowArea( tLNode *, double );
   void SetVoronoiVertices();
   void MakeFlow();
   void FlowUniform();
   void FlowSaturated();
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
   //double inDrArea;
   tInlet inlet;
   double mndrDirChngProb;
};

#endif
