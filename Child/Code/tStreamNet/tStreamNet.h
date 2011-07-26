//-*-c++-*-

/**************************************************************************/
/**
**  @file tStreamNet.h
**  @brief Header file for class tStreamNet and related class tInlet.
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
**  TODO: decouple from landscape erosion functions; make this a purely
**  hydrologic entity. Currently tStreamNet is "hard wired" to operate
**  with a tLNode; it could be recoded to operate with a "hydrologic"
**  node or any descendents.
**
**  Modifications:
**   - added new class tParkerChannels to implement Parker-Paola
**     channel geometry model (GT 6/01)
**   - added option to tInlet to have sediment influx calculated for each
**     storm, based on a prescribed slope and bed grain-size distribution.
**     (GT 5/06)
**
**  $Id: tStreamNet.h,v 1.65 2006-11-12 23:39:46 childcvs Exp $
*/
/**************************************************************************/

#ifndef TSTREAMNET_H
#define TSTREAMNET_H

#include "../Definitions.h"
#include "../Classes.h"
#include "../tUplift/tUplift.h"
#include "../Erosion/erosion.h"
#include "../MeshElements/meshElements.h"
#include "../tMesh/tMesh.h"
#include "../tLNode/tLNode.h"
#include "../tInputFile/tInputFile.h"
#include "../tStorm/tStorm.h"
#include "../globalFns.h"
#include <vector>


/**************************************************************************/
/**
** @class tInlet
**
** A tInlet represents a point source of water and sediment, normally at
** one side of a computational mesh. Each tInlet maintains a pointer to
** the inlet node, a pointer to the mesh object, and values of drainage
** area, total sediment load, and sediment load by grain size.
**
*/
/**************************************************************************/
class tInlet
{
    friend class tStreamNet;
  tInlet( const tInlet&, tMesh<tLNode>* );
  tInlet& operator=(const tInlet&);
public:
    tInlet();
    tInlet( tMesh< tLNode > *, const tInputFile & );
    ~tInlet();
    void FindNewInlet();
    double getInSedLoad() const;
    double getInSedLoad( size_t );
    tArray< double > getInSedLoadm() const;
    void setInSedLoad( double );
    void setInSedLoad( size_t, double );
    double getInDrArea() const;
    void setInDrArea( double );
    tLNode *getInNodePtr();
    void setInNodePtr( tLNode * );

private:
    tLNode *innode;   // ptr to inlet node
    double inDrArea;  // drainage area at inlet
    double inSedLoad; // total sediment load
    tArray< double > inSedLoadm; // incoming sediment load if multi-sizes
    tMesh< tLNode > *meshPtr;  // ptr to mesh
    
    // IN PROGRESS: ALTERNATIVE HANDLING OF SED FEED, MAY 06
    bool optCalcSedFeed;  // Option to calculate sediment influx based on prescribed gradient & bed grain sizes
    double inletSlope;    // Fixed slope for inlet, if used
    tArray<double> inletSedSizeFraction;  // Proportion of each sediment class on bed at inlet
    
};


/**************************************************************************/
/**
**  @class tParkerChannels
**
**  A "Parker Channel" is one that obeys the Parker-Paola self-formed
**  channel hypothesis for gravel-bed streams. The hypothesis states
**  that a gravel channel with mobile bed and banks will tend to adjust
**  its width such that the ratio of bankfull bed shear stress to critical
**  shear stress for mobilizing mean-size bed sediment (D50) is equal to
**  a constant with a value between 1.2 - 1.4.
**
**  The tParkerChannels class implements this concept, and provides an
**  alternative to empirical "regime" channel geometry.
**
**  The class provides a constructor function that reads the necessary
**  parameters from the input file and sets up initial values. The
**  calculation of channel width and depth is performed in the
**  CalcChanGeom() member function.
**
**  In the present implementation, width is computed individually for
**  every storm. This approximation is convenient, but leaves something
**  to be desired -- it implies that channel form adjusts immediately
**  to each discharge event. An alternative would be to set bankfull
**  width according to the Parker-Paola hypothesis, and actual width
**  according to an empirical width-discharge function.
**
**  Created: June, 2001, GT
**
**  References:
**     Paola, C., Heller, P.L., and Angevine, C.L., 1992, The large-scale
**       dynamics of grain-size variation in alluvial basins, 1: theory:
**       Basin Research, v. 4, p. 73-90.
**    Parker, 1978a, Self-formed straight rivers with equilibrium banks
**      and mobile bed. Part 1. The sand-silt river: J. Fluid Mech, v. 89,
**      p.109-125
**    Parker, 1978b, Self-formed straight rivers with equilibrium banks and
**      mobile bed. Part 2. The gravel river: J. Fluid Mech, v. 89, p.127-148.
**    Parker, G., 1979, Hydraulic geometry of active gravel rivers: Journal
**      of Hydraulics Division, American Society of Civil Engineering,
**      v. 105, p. 1185-1201.
**
**  Modifications:
**
*/
/**************************************************************************/
class tParkerChannels
{
 public:
  tParkerChannels( const tParkerChannels & );
  tParkerChannels( const tInputFile &infile );
  void CalcChanGeom( tMesh<tLNode> *meshPtr );

 private:
  double mdPPfac;    // Multiplicative factor for chan. width
  double mdPPexp1;   // Exponent on slope in width equation
  double mdPPexp2;   // Exponent on D50 in width equation (used if multi)
  double mdRough;    // Roughness (eg, Manning's n), used for depth
  double mdDepthexp; // Exponent used in computing depth via Manning/Chezy
  int miNumGrainSizeClasses;  // Number of grain size classes
  tArray<double> mD50BySizeClass; // Array recording D50 of each size class

private:
  tParkerChannels();
};

inline tParkerChannels::tParkerChannels( const tParkerChannels &orig )
  : mdPPfac(orig.mdPPfac), mdPPexp1(orig.mdPPexp1), mdPPexp2(orig.mdPPexp2), 
    mdRough(orig.mdRough), mdDepthexp(orig.mdDepthexp), 
    miNumGrainSizeClasses(orig.miNumGrainSizeClasses), 
    mD50BySizeClass(orig.mD50BySizeClass) 
{}

/**************************************************************************/
/**
**  @class tStreamNet
**
**  The tStreamNet class handles routing of water across a landscape
**  surface, including partitioning of flow into surface and subsurface
**  components. Methods include setting flow directions, computing total
**  contributing areas and discharges, computing runoff, resolving
**  drainage for closed depressions ("lakes"), calculating hydraulic
**  geometry, and sorting nodes in upstream-to-downstream order.
**
**  tStreamNet also includes methods to compute edge slopes, though this
**  should be moved to mesh utilities (TODO). Some meander functionality
**  is here, and this should also be moved (TODO).
**
**  Modifications:
**   - 3/99 GT added data member bankfullevent: precip rate corresponding
**     to bankfull discharge.
**   - 6/99 GT removed unused parameter mndrchngprob and ref's to it
**   - 1/00 GT added RouteFlowKinWave and two data parameters,
**     mdKinWaveExp and mdKinWaveRough
**   - 2/00 GT added DensifyMeshDrArea plus two data members to support
**     dynamic node addition in areas of high drainage area
**   - 6/01 GT added miChannelType and mpParkerChannels parameters to
**     implement alternative "Parker Channels"
**     (see new class tParkerChannels)
**
*/
/**************************************************************************/
class tStreamNet
{
    tStreamNet& operator=(const tStreamNet&);
    tStreamNet();
public:
    typedef enum {
      kHortonian = 0,       // Option for uniform infilt-excess runoff
      kSaturatedFlow1 = 1,  // Option for sat-excess runoff w/ return flow
      kSaturatedFlow2 = 2,  // Option for sat-excess runoff w/o return flow
      kConstSoilStore = 3,  // Option for "bucket"-type flow generation
      k2DKinematicWave = 4, // Option for 2D steady kinematic wave multi-flow
      kHydrographPeakMethod = 5,  // Option for hydrograph peak method
      kSubSurf2DKinematicWave = 6 // Option for kinematic wave with Darcy's Law
    } kFlowGen_t;

    tStreamNet( tMesh< tLNode > &, tStorm &, const tInputFile & );
  tStreamNet( const tStreamNet&, tStorm *, tMesh<tLNode>* );
    ~tStreamNet();
    void ResetMesh( tMesh< tLNode > & );
    const tMesh< tLNode > *getMeshPtr() const;
  inline void setMeshPtr( tMesh<tLNode>* ptr ) {meshPtr = ptr;}
    tMesh< tLNode > *getMeshPtrNC();
    const tStorm *getStormPtr() const;
  inline void setStormPtr( tStorm* ptr ) {stormPtr = ptr;}
    tStorm *getStormPtrNC();
    kFlowGen_t getFlowGenOpt() const;
    bool getFillLakesOpt() const;
    double getRainRate() const;
    double getTransmissivity() const;
    double getInfilt() const;
    double getSoilStore() const;
    double getInDrArea() const;
    double getInSedLoad() const;
    tArray< double > getInSedLoadm() const;
    tLNode const *getInletNodePtr() const;
    tLNode *getInletNodePtrNC();
    double getInletSlope() const;   // Returns slope of inlet node if fixed slope option is used
    tArray< double > getInletSedSizeFraction() const;  // Returns sed size fractions on bed of inlet node if appl.
    void setFlowGenOpt( kFlowGen_t );
    void setFillLakesOpt( bool );
    void setRainRate( double );
    void setTransmissivity( double );
    void setInfilt( double );
    void setInDrArea( double );
    void setInSedLoad( double );
    void setInSedLoadm( size_t, double );
    void setInletNodePtr( tLNode * );
    void UpdateNet( double time );
    void UpdateNet( double time, tStorm & );
    void CheckNetConsistency();
    int CheckNetConsistencyFlowPath( tLNode **);
    void CheckMeander() const;
    void CalcSlopes();
    void InitFlowDirs();
    void ReInitFlowDirs();
    void FlowDirs();
    void DrainAreaVoronoi();
//   void DrainAreaVoronoiMFD();
    void FlowPathLength();
    void RouteFlowHydrographPeak();
    void MakeFlow( double tm );
    void FlowUniform();
    void FlowSaturated1();
    void FlowSaturated2();
    void FlowBucket();
    void RouteFlowKinWave( double );
    void FillLakes();
    bool FindLakeNodeOutlet( tLNode * ) const;
    void SortNodesByNetOrder( bool optMultiFlow=false );
    //find hydraulic and channel geometries, respectively;
    //FindHydrGeom is contingent upon current storm conditions
    //and storm variability;
    //FindChanGeom is based on the 1.5-yr storm event,
    //or the mean rainrate if no variability:
    void FindChanGeom();
    void FindHydrGeom();
    void DensifyMeshDrArea( double time=0.0 );  // Densifies mesh locally
    void ShowMeanderNeighbours(int) const;
  // find streamlines from points specified in input file:
  void FindStreamLines( const tInputFile &, tPtrList< tLNode > &, bool lvFEs = false );

protected:
    tLNode *BuildLakeList( tPtrList< tLNode > &, tLNode *);
    void FillLakesFlowDirs(tPtrListIter< tLNode > &, tLNode *) const;
    inline static void RouteFlowArea( tLNode *, double );
    inline static void RouteRunoff( tLNode *, double, double );
    static void RouteError( tLNode * ) ATTRIBUTE_NORETURN;
//   void RouteFlowAreaMultipleDirections( tLNode*, double );
	bool FlowDirBreaksMeanderChannel( tLNode *, tEdge * ) const;

    typedef enum {
      kRegimeChannels = 1,
      kParkerChannels = 2
    } kChannelType_t;
    static kChannelType_t IntToChannelType( int );

    tMesh< tLNode > * meshPtr;  // ptr to mesh
    tStorm *stormPtr; // ptr to storm object (for getting precip)
    kFlowGen_t miOptFlowgen; // option for runoff generation & routing method
    bool filllakes;   // option for filling lakes
    bool optrainvar;  //flag w/ 1=>varying storms=>hydraulic geom != chan. geom.
    double kwds, ewds, ewstn;//coefs & exps for dwnstrm & at-a-stn hydr. width
    double kdds, edds, edstn;//coefs & exps for dwnstrm & at-a-stn hydr. depth
    double knds, ends, enstn;//coefs & exps for ds & at-a-stn hydr. roughness
    double klambda, elambda; //coef & exp for downstrm bank roughness length
    double rainrate;      // current rainfall rate
    double bankfullevent; // rainfall rate corresponding to bankfull event
    double trans;         // soil transmissivity
    double infilt;        // soil infiltration capacity (also K_sat)
    double soilStore;     // soil water storage, depth equiv (m)
    tInlet inlet;         // inlet
    bool optSinVarInfilt;  // opt for sinusoidal variation in infilt cap
    kChannelType_t miChannelType; // code for type of channels: "regime", "parker"
    tParkerChannels *mpParkerChannels;  // -> tParkerChannels object
    double infilt_dev;    // max +/- variation from mean infilt cap
    double infilt0;    // mean infilt cap
    double twoPiLam;   // Parameter for sinusoidal variation: 2pi / period
    double mdKinWaveExp;  // Depth-disch exponent for kinematic wave routing
    double mdKinWaveRough; // Roughness coef, units m^(1/e)-2 yr, e=above exp
    double mdMeshAdaptMinArea;  // Dr area threshold for densifying mesh
    double mdMeshAdaptMaxVArea; // Max voronoi area for nodes above threshold
    double mdHydrgrphShapeFac;  // "Fhs" for hydrograph peak method
    double mdFlowVelocity;      // Runoff velocity for computing travel time
  bool optVariableTransmissivity; // option for soil depth-dependent transmissivity
//   bool optMultipleFlowDirections; // option for flow routing via MFD algorithm

  void DebugShowNbrs( tLNode * theNode ) const;  // debugging function shows neighbor nodes
};

#endif






