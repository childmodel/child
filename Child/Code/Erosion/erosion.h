/***************************************************************************\
**
**  erosion.h
**
**  Header file for objects related to sediment transport and detachment.
**
**  This file includes a group of sediment transport objects (tSedTransX)
**  and bed material (bedrock, cohesive sediment, vegetation-covered
**  regolith, etc) detachment objects (tDetachX), where X describes
**  the particular type of transport/detachment function the object
**  implements.
**
**  tSedTrans objects contain data and functions needed to compute
**  runoff-driven sediment transport capacity at a point, while
**  tDetach objects do the same for bed material detachment. Sediment
**  transport/erosion data and routines can be quite simple, so why group
**  them into their own objects? The reason is that it allows for maximum
**  flexibility. The idea is that there can many different types of
**  transport and erosion objects defined, ranging from a simple power-law
**  to more complex formulas involving thresholds, multiple grain-sizes,
**  etc. Each object knows all the data it needs, and knows how to read
**  those data in from an input file and how to compute the transport
**  capacity for a node. (For performance reasons, inheritance is not
**  used; rather, each transport function is its own object type.)
**  The idea is that if you want to modify the code to add your own
**  function type, you can do so by simply creating a new transport or
**  erosion object and then defining "tSedTrans" or "tDetachment"
**  as your particular object type in the file tStreamNet.h, for example:
**
**    #define tSedTrans tSedTransNewImproved
**
**  The only requirements are that each transport object must include a
**  constructor that reads in all the relevant parameters and a
**  TransportCapacity() function that returns the total transport capacity;
**  for detachment objects, a constructor and a DetachCapacity() function
**  must be provided.
**
**  Note that these functions assume that each tLNode includes both
**  an elevation and a "potential change in elevation", dz, that
**  is used to store intermediate solutions in a numerical solution
**  algorithm. Thus, slope is always computed as (zi+dzi - zj+dzj)/dx
**  (possibly with dz=0 of course).
**
**    Created 1/98 gt
**
**  $Id: erosion.h,v 1.4 1998-01-29 18:53:51 stlancas Exp $
\***************************************************************************/

#ifndef EROSION_H
#define EROSION_H

#include "../Definitions.h"
#include "../Classes.h"
#include "../tInputFile/tInputFile.h"
#include "../tLNode/tLNode.h"

/***************************************************************************\
**  class tSedTransPwrLaw
**
**  Manages data and routines to compute sediment transport capacity as a
**  simple power function of slope and total discharge (channel width and
**  flow depth are implicit in the power-law derivation):
**    Qs = kf Q^m S^n
\***************************************************************************/
class tSedTransPwrLaw
{
  public:
   tSedTransPwrLaw( tInputFile &infile );
   float TransCapacity( tLNode * n );

  private:
   float kf;  // Transport capacity coefficient
   float mf;  // Exponent on total discharge
   float nf;  // Exponent on slope
};


/***************************************************************************\
**  class tBedErodePwrLaw
**
**  Assumes bedrock detachment proportional to a power function of slope
**  and total discharge. Regolith is considered infinitely detachable, so
**  that the total detachable depth of material over a duration dt is
**  equal to the thickness of any regolith (alluvium) present plus
**    Dc = kb Q^mb S^nb dt
\***************************************************************************/
class tBedErodePwrLaw
{
  public:
   tBedErodePwrLaw( tInputFile &infile );
     //Computes depth of potential erosion at node n over time interval dt
   float DetachCapacity( tLNode * n, float dt );
     //Returns an estimate of maximum stable & accurate time step size
   float SetTimeStep( tLNode * n );

  private:
   float kb;  // Erosion coefficient
   float mb;  // Exponent on total discharge
   float nb;  // Exponent on slope
};



#endif
