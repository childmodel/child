/**************************************************************************\
**
**  tStreamNet.h: Header file for class tStreamNet.
**
**  tStreamNet objects contain data and functions related to flow routing
**  and sediment transport across the landscape surface.
**
**  $Id: tStreamNet.h,v 1.6 1998-01-29 20:06:53 stlancas Exp $
\**************************************************************************/

#ifndef TSTREAMNET_H
#define TSTREAMNET_H

#include "../Definitions.h"
#include "../Classes.h"
#include "../Erosion/erosion.h"
#include "../tGrid/tGrid.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"
#include "../tStorm/tStorm.h"
 
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
   void setFlowGenOpt( int );
   void setFillLakesOpt( int );
   void setRainRate( double );
   void setTransmissivity( double );
   void setInfilt( double );
   void UpdateNet();
   void UpdateNet( tStorm & );
   void CalcSlopes();
   void CalcVAreas();
   void InitFlowDirs();
   void FlowDirs();
   void DrainAreaVoronoi();
   double VoronoiArea( tLNode * );
   void RouteFlowArea( tLNode *, double );
   void SetVoronoiVertices();
   void MakeFlow();
   void FlowUniform();
   void FlowSaturated();
   void FillLakes();
   int FindLakeNodeOutlet( tLNode * );
   void SortNodesByNetOrder();
   void ErodeDetachLim( double dtg );
   
protected:
   
   tGrid< tLNode > * gridPtr;
   tStorm *stormPtr;
   int flowgen;
   int filllakes;
   double rainrate;
   double trans;
   double infilt;    
   tBedErodePwrLaw bedErode;
};

#endif
