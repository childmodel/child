//-*-c++-*-

/***************************************************************************/
/**
 **  @file erosion.h
 **  @brief Header file for objects related to sediment transport and
 **         detachment.
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
 **  constructor that reads in all the relevant parameters and a constructor
 **  and a DetachCapacity() function must be provided.
 **
 **  Note that these functions assume that each tLNode includes both
 **  an elevation and a "potential change in elevation", dz, that
 **  is used to store intermediate solutions in a numerical solution
 **  algorithm. Thus, slope is always computed as (zi+dzi - zj+dzj)/dx
 **  (possibly with dz=0 of course).
 **
 **    Created 1/98 gt
 **
 **    Modifications:
 **     - added an erosion-rate-based adaptive mesh capability, with
 **       a new DensifyMesh function and supporting data member
 **       mdMeshAdaptMaxFlux (gt 2/2000)
 **     - 2/02 new class tSedTransPwrLawMulti to handle multi-size
 **       transport in power-law (excess-shear stress) formulation
 **       (GT)
 **     - Added variants of power law detachment and transport formulae
 **       of form tau^p - tauc^p rather than (tau - tauc)^p (GT 4/02)
 **     - Added Bridge-Dominic form of Bagnold bedload transport formula
 **       (single-size) (GT 5/02)
 **     - Added codes to go along with erosion & transport options, to
 **       enable checking against user-specified options (GT 7/02)
 **
 **  $Id: erosion.h,v 1.56 2005/07/21 19:30:24 childcvs Exp $
 */
/***************************************************************************/

#ifndef EROSION_H
#define EROSION_H

#include "../Definitions.h"
#include "../Classes.h"
#include "../tArray/tArray.h"
#include "../tArray/tArray2.h"
#include "../tInputFile/tInputFile.h"
#include "../tLNode/tLNode.h"
#include "../tUplift/tUplift.h"
#include "../tStreamNet/tStreamNet.h"
#include "../tRunTimer/tRunTimer.h"
#include "../tVegetation/tVegetation.h"

/***************************************************************************/
/*
**  @class tEquilibCheck
**
**  Enables dynamic equilibrium checking, both short- and a specified long-
**  term. The idea is to find the rate of total volume change over the mesh.
**  With meandering, this will never be zero over the short term, but we
**  should be able to find an average over the long term.
**
**  Needs to look at the mesh; can either make it a template or just use
**  the mesh of tLNodes. Do the latter...
*/
/***************************************************************************/
class tEquilibCheck
{
  tEquilibCheck(const tEquilibCheck&);
  tEquilibCheck& operator=(const tEquilibCheck&);
public:
  tEquilibCheck();
  tEquilibCheck( tMesh< tLNode > &, tRunTimer & );
  tEquilibCheck( tMesh< tLNode > &, tRunTimer &, const tInputFile & );
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
  tList< tArray2< double > > massList; //linked list of arrays: (time, mesh mass)
  //'mass' is misnomer--actually mean elev.
  double longRate;
  double shortRate;
};

/***************************************************************************/
/**
 **  @class tSedTrans
 **
 **  Abstract Base Class for transport law
 */
/***************************************************************************/
class tSedTrans
{
public:
  virtual ~tSedTrans() {}
  virtual double TransCapacity( tLNode *n ) = 0;
  virtual double TransCapacity( tLNode *n, int i, double weight) = 0;
};

/***************************************************************************/
/**
 **  @class tSedTransPwrLaw
 **
 **  Manages data and routines to compute sediment transport capacity as a
 **  simple power function of slope and total discharge (channel width and
 **  flow depth are implicit in the power-law derivation):
 **    Qs = kf ( tau - tauc ) ^ pf,   tau = kt (Q/W)^mb S^nf
 */
/***************************************************************************/
class tSedTransPwrLaw : public tSedTrans
{
  tSedTransPwrLaw();
public:
  tSedTransPwrLaw( const tInputFile &infile );
  double TransCapacity( tLNode * n );
  double TransCapacity( tLNode *n, int i, double weight);

private:
  double kf;  // Transport capacity coefficient
  double kt;  // Shear stress coefficient
  double mf;  // Exponent on total discharge
  double nf;  // Exponent on slope
  double pf;  // Excess shear exponent
  double tauc; // Entrainment threshold
};

/***************************************************************************/
/**
 **  @class tSedTransPwrLaw2
 **
 **  Manages data and routines to compute sediment transport capacity as a
 **  simple power function of slope and total discharge (channel width and
 **  flow depth are implicit in the power-law derivation):
 **    Qs = kf ( tau^pf - tauc^pf ),  tau = kt (Q/W)^mf S^nf
 */
/***************************************************************************/
class tSedTransPwrLaw2 : public tSedTrans
{
  tSedTransPwrLaw2();
public:
  tSedTransPwrLaw2( const tInputFile &infile );
  double TransCapacity( tLNode * n );
  double TransCapacity( tLNode *n, int i, double weight);

private:
  double kf;  // Transport capacity coefficient
  double kt;  // Shear stress coefficient
  double mf;  // Exponent on total discharge
  double nf;  // Exponent on slope
  double pf;  // Excess shear exponent
  double tauc; // Entrainment threshold
};

/***************************************************************************/
/**
 **  @class tSedTransPwrLawSimp
 **
 **  Manages data and routines to compute sediment transport capacity as a
 **  simple power function of slope and total discharge 
 **  This is the ultra simple version put in for use with the sediment-flux
 **  detachment models.
 **    Qc = kf Q^mf S^nf
 */
/***************************************************************************/
class tSedTransPwrLawSimp : public tSedTrans
{
public:
  tSedTransPwrLawSimp( const tInputFile &infile );
  double TransCapacity( tLNode * n );
  double TransCapacity( tLNode *n, int i, double weight);

private:
  double kf;  // Transport capacity coefficient
  double mf;  // Exponent on total discharge
  double nf;  // Exponent on slope
};


/***************************************************************************/
/**
 **  @class tSedTransBridgeDom
 **
 **  Manages data and routines to compute sediment transport capacity
 **  using the Bridge and Dominic (1984) version of the Bagnold bedload
 **  transport formula.
 **
 */
/***************************************************************************/
class tSedTransBridgeDom : public tSedTrans
{
  tSedTransBridgeDom();
public:
  tSedTransBridgeDom( const tInputFile &infile );
  double TransCapacity( tLNode * n );
  double TransCapacity( tLNode *n, int i, double weight);

private:
  double kf;  // Transport capacity coefficient
  double kt;  // Shear stress coefficient
  double mf;  // Exponent on total discharge
  double nf;  // Exponent on slope
  double tauc; // Entrainment threshold
  double sqrtTauc;  // Threshold value of U_* rho^0.5
};


/***************************************************************************/
/**
 **  @class tSedTransPwrLawMulti
 **
 **  Manages data and routines to compute sediment transport capacity for
 **  multiple grain size fractions, using an excess shear stress formulation
 **  a la Meyer-Peter & Mueller. Uses Komar-style hiding & protrusion
 **  function.
 **
 */
/***************************************************************************/
class tSedTransPwrLawMulti : public tSedTrans
{
  tSedTransPwrLawMulti();
public:
  tSedTransPwrLawMulti( const tInputFile &infile );
  double TransCapacity( tLNode * n );
  double TransCapacity( tLNode *n, int lyr, double weight );

private:
  double kf;  // Transport capacity coefficient
  double kt;  // Shear stress coefficient
  double mf;  // Exponent on total discharge
  double nf;  // Exponent on slope
  double pf;  // Excess shear exponent
  tArray<double> mdGrndiam;  // Grain diameters
  tArray<double> mdTauc; // Entrainment threshold
  int miNumgrnsizes;  // No. grain size classes
  double mdHidingexp; // Hiding/protrusion exponent
};


/**************************************************************************/
/**
 **  @class tSedTransWilcock
 **
 **  Manages data and routines to compute sediment transport capacity
 **  of a sand a gravel class (two grain sizes) using the sediment transport
 **  formula and critical shear stress rules developed by P. Wilcock (1997)
 */
/**************************************************************************/
class tSedTransWilcock : public tSedTrans
{
  tSedTransWilcock();
public:
  tSedTransWilcock( const tInputFile &infile );
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

/************************************************************************/
/**
 ** @class tSedTransMineTailings
 **
 ** Manages data and routines to compute sediment transport capacity
 ** using the parameters and equation from Willgoose and Riley (1998).
 ** This study was performed on mine tailings in an Australian Uranium
 ** mine.  Don't have a critical shear stress method, so just use that of
 ** of Wilcock for sand and gravel.
 **
 ** added 04/2000 ng
 */
/************************************************************************/
class tSedTransMineTailings : public tSedTrans
{
  tSedTransMineTailings();
public:
  tSedTransMineTailings( const tInputFile &infile );
  double TransCapacity( tLNode * n ); // returns total volumetric load
  double TransCapacity( tLNode *n, int i, double weight);

private:
  //all of these are just the same as for tSedTransWilcock since using
  //the same critical shear stress method
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

/***************************************************************************/
/**
 **  @class tBedErode
 **
 **  Abstract Base Class for detachment law
 */
/***************************************************************************/
class tBedErode
{
public:
  virtual ~tBedErode() {}
  //Computes depth of potential erosion at node n over time interval dt
  virtual double DetachCapacity( tLNode * n, double dt ) = 0 ;
  //Computes rate of potential erosion of layer i at node n
  virtual double DetachCapacity( tLNode * n, int i ) = 0 ;
  //Computes rate of erosion at node n
  virtual double DetachCapacity( tLNode * n ) = 0 ;
  //Returns an estimate of maximum stable & accurate time step size
  virtual double SetTimeStep( tLNode * n ) = 0 ;
};

/***************************************************************************/
/**
 **  @class tBedErodePwrLaw
 **
 **  Assumes bedrock detachment proportional to a power function of slope
 **  and total discharge. Regolith is considered infinitely detachable, so
 **  that the total detachable depth of material over a duration dt is
 **  equal to the thickness of any regolith (alluvium) present plus
 **    Dc = kb Q^mb S^nb dt
 */
/***************************************************************************/
class tBedErodePwrLaw : public tBedErode
{
  tBedErodePwrLaw();
public:
  tBedErodePwrLaw( const tInputFile &infile );
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
  double kt;  // Shear stress (or stream power) coefficient
  double mb;  // Exponent on total discharge
  //double ma;  // Exponent on drainage area (can differ from mb!)
  double nb;  // Exponent on slope
  double pb;  // Exponent on excess erosion capacity (e.g., shear stress)
  double taucd;  // Erosion threshold
};


/***************************************************************************/
/**
 **  @class tBedErodePwrLaw2
 **
 **  This is a variation of tBedErodePwrLaw that differs in the following
 **  respect:
 **    tBedErodePwrLaw:  erorate ~ ( tau - taucrit ) ^ pb
 **    tBedErodePwrLaw2: erorate ~ tau^pb - taucrit^pb
 **
 **  In other words, this function computes a more analytically tractable
 **  form of the power-law erosion equation.
 **
 **  Created: April 2002 (GT)
 **
 */
/***************************************************************************/
class tBedErodePwrLaw2 : public tBedErode
{
  tBedErodePwrLaw2();
public:
  tBedErodePwrLaw2( const tInputFile &infile );
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
  double kt;  // Shear stress (or stream power) coefficient
  double mb;  // Exponent on total discharge
  //double ma;  // Exponent on drainage area (can differ from mb!)
  double nb;  // Exponent on slope
  double pb;  // Exponent on excess erosion capacity (e.g., shear stress)
  double taucd;  // Erosion threshold
};

/***************************************************************************/
/**
 **  @class tBedErodeAParabolic1
 **
 **  This function considers a variable erodibility as a function of the
 **  ratio of incoming sediment load to transport capacity.
 **  fqs = variable erodibility = f(Qs/Qc) following the function described
 **  by Whipple and Tucker 2002, with a linear beginning to account for the 
 **  case when there is no incoming sediment load, but there still may be
 **  erosion.
 **
 **    tBedErodeAParabolic1:  erorate ~ K fqs tau 
 **
 **  Simplifications -> not using a threshold term - ever.
 **  Assume that -> sediment transport capacity is already calculated.
 **              -> sediment routing is properly taken care of.
 **                  
  **  Created: FEB 2004 NG
 **
 */
/***************************************************************************/
class tBedErodeAParabolic1 : public tBedErode
{
public:
  tBedErodeAParabolic1( const tInputFile &infile );
  //Computes depth of potential erosion at node n over time interval dt
  double DetachCapacity( tLNode * n, double dt );
  //Computes rate of potential erosion of layer i at node n
  double DetachCapacity( tLNode * n, int i );
  //Computes rate of erosion at node n
  double DetachCapacity( tLNode * n );
  //Returns an estimate of maximum stable & accurate time step size
  double SetTimeStep( tLNode * n );

private:
  double kb;  // Erosion coefficient - get this from layer, but used for ts
  double mb;  // Exponent on total discharge
  double nb;  // Exponent on slope
  double beta;  // fraction of sediment which is bedload

};

/***************************************************************************\
 **  class tBedErodeGeneralFQS
 **
 **  E=K Qs/W (1- Qs/Qc ) A^m S^n
 **
 **  This function will never erode unless there is some source of sediment,
 **  other than fluvial, to start up erosion.  This sediment comes from
 **  diffusion and in order for proper accounting to take place, this function
 **  should be used along with detachErode2.
 **
 **  ASSUMPTIONS - incoming sediment load from a different process in upper 
 **                reaches - diffusion.
 **              - transport capacity is set elsewhere.
 **              - use hydraulic width => W = kb Q^b - calculated elsewhere.
 **
 **  created 01-2005 NMG
 ***************************************************************************/
class tBedErodeGeneralFQS : public tBedErode
{
public:
   tBedErodeGeneralFQS( const tInputFile &infile );
     //Computes depth of potential erosion at node n over time interval dt
   double DetachCapacity( tLNode * nd, double dt );
   //Computes rate of potential erosion of layer i at node n 
   double DetachCapacity( tLNode * nd, int i );
     //Computes rate of erosion at node n
   double DetachCapacity( tLNode * nd );
     //Does nothing but included in just in case it's needed in the future.
   double SetTimeStep( tLNode * nd );

  private:
   double m;
   double n;
   double K;
   double beta;
   
};

/***************************************************************************/
/**
 **  @class tErosion
 **
 **  Manages data and routines related to various aspects of erosion.
 **
 **  (class added by gt 3/98; routines were previously collected under
 **  tStreamNet).
 */
/***************************************************************************/
class tErosion
{
   tErosion(const tErosion&);
   tErosion& operator=(const tErosion&);
   tErosion();
public:
   tErosion( tMesh< tLNode > *, const tInputFile & );
   ~tErosion();
   void ErodeDetachLim( double dtg, tStreamNet *, tVegetation * );
   void ErodeDetachLim( double dtg, tStreamNet *, tUplift const * );
   void StreamErode( double dtg, tStreamNet * );
   void StreamErodeMulti( double dtg, tStreamNet *, double time);
   void DetachErode( double dtg, tStreamNet *, double time, tVegetation * pVegetation );
   void DetachErode2( double dtg, tStreamNet *, double time, tVegetation * pVegetation );
   void Diffuse( double dtg, bool detach );
   void UpdateExposureTime( double dtg);
   void DensifyMesh( double time );
   
private:
   tMesh<tLNode> *meshPtr;    // ptr to mesh
   tBedErode *bedErode;        // bed erosion object
   tSedTrans *sedTrans;        // sediment transport object
   double kd;                 // Hillslope transport (diffusion) coef
   double difThresh;          // Diffusion occurs only at areas < difThresh
   double mdMeshAdaptMaxFlux; // For dynamic point addition: max ero flux rate
   
};


#endif
