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
**  $Id: erosion.h,v 1.22 1999-04-04 21:34:26 gtucker Exp $
\***************************************************************************/

#ifndef EROSION_H
#define EROSION_H

#include "../Definitions.h"
#include "../Classes.h"
#include "../tArray/tArray.h"
#include "../tInputFile/tInputFile.h"
#include "../tLNode/tLNode.h"
#include "../tUplift/tUplift.h"
#include "../tStreamNet/tStreamNet.h"
#include "../tRunTimer/tRunTimer.h"

#define tSedTrans tSedTransPwrLaw


/***************************************************************************\
**  class tEquilibCheck
**
**  Enables dynamic equilibrium checking, both short- and a specified long-
**  term. The idea is to find the rate of total volume change over the mesh.
**  With meandering, this will never be zero over the short term, but we
**  should be able to find an average over the long term.
**
**  Needs to look at the mesh; can either make it a template or just use
**  the mesh of tLNodes. Do the latter...
\***************************************************************************/
class tEquilibCheck
{
public:
    tEquilibCheck();
    tEquilibCheck( tMesh< tLNode > &, tRunTimer & );
    tEquilibCheck( tMesh< tLNode > &, tRunTimer &, tInputFile & );
    ~tEquilibCheck();
    double getLongTime() const; //get the interval for long-term change
    void setLongTime( double ); //set the interval for long-term change
    const tMesh< tLNode > *getMeshPtr() const;
    tMesh< tLNode > *getMeshPtrNC();
    void setMeshPtr( tMesh< tLNode > & );
    void setMeshPtr( tMesh< tLNode > * );
    const tRunTimer *getTimePtr() const;
    tRunTimer *getTimePtrNC();
    void setTimePtr( tRunTimer & );
    void setTimePtr( tRunTimer * );
    double getLongRate() const;
    double getShortRate() const;
    double FindIterChngRate(); //find the change rate since last call
    double FindLongTermChngRate(); //find change rate over pre-set interval
    double FindLongTermChngRate( double ); //find rate over given interval
private:
    tMesh< tLNode > *meshPtr; //ptr to tMesh
    tRunTimer *timePtr; //ptr to tRunTimer
    double longTime; //'long' time interval
    tList< tArray< double > > massList; //linked list of arrays: (time, mesh mass)
                                        //'mass' is misnomer--actually mean elev.
    double longRate;
    double shortRate;
};


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
   double TransCapacity( tLNode * n );
   double TransCapacity( tLNode *n, int i, double weight);

  private:
   double kf;  // Transport capacity coefficient
   double mf;  // Exponent on total discharge
   double nf;  // Exponent on slope
};


/**************************************************************************\
**  class tSedTransWilcock
**
**  Manages data and routines to compute sediment transport capacity
**  of a sand a gravel class (two grain sizes) using the sediment transport
**  formula and critical shear stress rules developed by P. Wilcock (1997)
\**************************************************************************/
class tSedTransWilcock
{
public:
   tSedTransWilcock( tInputFile &infile );
   double TransCapacity( tLNode * n ); // returns total volumetric load
   double TransCapacity( tLNode *n, int i, double weight);
   //returns total volumetric load
   
private:
   double taudim;
   double refs;
   double refg;
   double lowtaucs;
   double lowtaucg;
   double sandb;
   double hightaucs;
   double hightaucg;
   double sands;
   double gravb;
   double gravs;
   tArray< double > grade;
   
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
   double DetachCapacity( tLNode * n, double dt );
   //Computes rate of potential erosion of layer i at node n 
   double DetachCapacity( tLNode * n, int i );
     //Computes rate of erosion at node n
   double DetachCapacity( tLNode * n );
     //Returns an estimate of maximum stable & accurate time step size
   double SetTimeStep( tLNode * n );

  private:
   double kb;  // Erosion coefficient
   double mb;  // Exponent on total discharge
   double nb;  // Exponent on slope
};


/***************************************************************************\
**  class tErosion
**
**  Manages data and routines related to various aspects of erosion.
**
**  (class added by gt 3/98; routines were previously collected under
**  tStreamNet).
\***************************************************************************/
class tErosion
{
public:
    tErosion( tMesh< tLNode > *, tInputFile & );
    void ErodeDetachLim( double dtg );
    void ErodeDetachLim( double dtg, tUplift * );
    void StreamErode( double dtg, tStreamNet * );
    void StreamErodeMulti( double dtg, tStreamNet *, double time);
    void DetachErode( double dtg, tStreamNet *, double time);
    double TransportCapacity(tLNode * n );
    void Diffuse( double dtg, int detach );
    void UpdateExposureTime( double dtg);

private:
    tMesh<tLNode> *meshPtr;    // ptr to mesh
    tBedErodePwrLaw bedErode;  // bed erosion object
    tSedTransWilcock sedTrans; // sediment transport object 
    double kd;                 // Hillslope transport (diffusion) coef

};


#endif
