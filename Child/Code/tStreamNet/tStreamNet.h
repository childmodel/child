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
**  $Id: tStreamNet.h,v 1.25 1998-07-20 22:09:42 gtucker Exp $
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

double DistanceToLine( double x2, double y2, double a, double b, double c );
double DistanceToLine( double x2, double y2, tNode *p0, tNode *p1 );

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
   void UpdateNet( double time );
   void UpdateNet( double time, tStorm & );
   void CheckNetConsistency();
   void CalcSlopes();
   void InitFlowDirs();
   void FlowDirs();
   void DrainAreaVoronoi();
   void RouteFlowArea( tLNode *, double );
   void RouteRunoff( tLNode *, double, double );
   void SetVoronoiVertices();
   void MakeFlow( double tm );
   void FlowUniform();
   void FlowSaturated1();
   void FlowSaturated2();
   void FlowBucket();
   void FillLakes();
   int FindLakeNodeOutlet( tLNode * );
   void SortNodesByNetOrder();
   int DamBypass( tLNode * );
   //find hydraulic and channel geometries, respectively;
   //FindHydrGeom is contingent upon current storm conditions
   //and storm variability;
   //FindChanGeom is based on the 1.5-yr storm event,
   //or the mean rainrate if no variability:   
   void FindChanGeom();
   void FindHydrGeom();
   
protected:
   tGrid< tLNode > * gridPtr;
   tStorm *stormPtr;
   int flowgen;
   int filllakes;
   int optrainvar;  //flag w/ 1=>varying storms=>hydraulic geom != chan. geom.
   double kwds, ewds, ewstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. width
   double kdds, edds, edstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. depth
   double knds, ends, enstn;//coeff's & exp's for dwnstrm & at-a-stn hydr. roughness
   double klambda, elambda; //coeff & exp for downstrm bank roughness length
   double rainrate;
   double trans;
   double infilt;
   double soilStore;        /* soil water storage, depth equiv (m) */
   //double inDrArea;
   tInlet inlet;
   double mndrDirChngProb;
    int optSinVarInfilt;  // opt for sinusoidal variation in infilt cap
    double infilt_dev;    // max +/- variation from mean infilt cap
    double infilt0;    // mean infilt cap
    double twoPiLam;   // Parameter for sinusoidal variation: 2pi / period

};

#endif






