//-*-c++-*-

/************************************************************************/
/**
 **  @file tLNode.h
 **  @brief Header file for derived class tLNode and its member classes
 **
 **  A tLNode object encapsulates data and methods for a single landscape
 **  node. It inherits basic data and methods from the generic tNode
 **  class (see meshElements.h/.cpp). Embedded within tLNode are four
 **  small classes, tChannel, tSurface, tRegolith, and tRock, that
 **  contain data specific to channel geometry, surface characteristics,
 **  regolith properties, and rock properties (many of the properties in
 **  these embedded objects are now obsolete and need to be surgically
 **  removed). Embedded within tChannel is an additional class,
 **  tMeander, that contains node-data related to stream meandering.
 **
 **  Changes:
 **    - GT commented out inclusion of run timer, 1/99
 **    - removed obsolete class tErode, GT 4/99
 **    - added data mbr tau to track shear stress (or stream pwr, etc),
 **      inlined most tLayer fns. GT 1/00
 **    - GT made these changes as part of addition of vegetation module:
 **        - removed tSurface class
 **        - added tauc data mbr and retrieval functions
 **        - added embedded tVegCover object and retrieval fn
 **          (Jan 2000)
 **
 **  $Id: tLNode.h,v 1.98 2007/08/07 02:23:56 childcvs Exp $
 */
/************************************************************************/

#ifndef TLNODE_H
#define TLNODE_H

#include <iostream>
#include <string.h>

#include "../Definitions.h"
#include "../tArray/tArray.h"
#include "../MeshElements/meshElements.h"
#include "../tList/tList.h"
#include "../tInputFile/tInputFile.h"
#include "../globalFns.h"
#include "../tVegetation/tVegetation.h"

class tStratNode;
class tStratGrid;

#include "../compiler.h"

#define kVeryHigh 100000  // Used in FillLakes


/** @class tLayer
    Layer records */

class tLayer
{
public:
  // 0 = bedrock; 1 = some mixture of sediments so there
  // may be multi-sizes.  although multiple sizes are stored for
  // bedrock, this is the texture of what the bedrock will break up
  // into, not what is there on the bed.
  // Later it may be used as a flag for alluvium vs. regolith, etc.
  typedef enum {
    kBedRock = 0,
    kSed = 1
  } tSed_t;

  inline tLayer();
  inline tLayer( size_t );
  inline tLayer( const tLayer & );
  inline const tLayer &operator=( const tLayer & );
  inline void setCtime( double );
  inline double getCtime() const;
  inline void setRtime( double );
  inline double getRtime() const;
  inline void setEtime( double );
  inline void addEtime( double );
  inline double getEtime() const;
  void setDepth( double );
  // NOTE setDepth will also update the depths in dgrade so that the
  // total depth is consistent and the original texture is kept.
  // If texture needs to also be changed, call setDgrade.  This will
  // also update depth and the texture will change appropriately.
  inline double getDepth() const;
  inline void setErody( double );
  inline double getErody() const;
  inline void setSed( tSed_t );
  inline tSed_t getSed() const;
  inline void setDgradesize( size_t );
  inline size_t getDgradesize() const;
  void setDgrade( size_t, double );
  inline double getDgrade( size_t ) const;
  inline tArray< double > const & getDgrade() const;
  inline void addDgrade(size_t, double);
  inline void addPaleoCurrent(double, double);
  inline void setPaleoCurrent( double );
  inline double getPaleoCurrent() const;

protected:
  double ctime; // time of creation of layer
  double rtime; // most recent time (time steps) that there was erosion/depo.
  double etime; // exposure time, i.e. time material spent at surface
  double depth; // total depth of the layer
  double erody; // erodibility of layer (varies with material)
  tSed_t sed;   // type of sediment
  tArray< double > dgrade; // depth of each size if multi size
  double paleocurrent;
};


/*************************************************************************
 **
 **  Inlined functions for tLayer
 **
 *************************************************************************/

/*************************************************************************
 **  tLayer::tLayer : Constructor function for tLayer
 *************************************************************************/
inline tLayer::tLayer () :
  ctime(0.), rtime(0.), etime(0.),
  depth(0.), erody(0.), sed(kBedRock),
  dgrade(), paleocurrent(-1.)
{
  if (0) //DEBUG
    std::cout << "tLayer( )" << std::endl;
}

inline tLayer::tLayer ( size_t num ) :
  ctime(0.), rtime(0.), etime(0.),
  depth(0.), erody(0.), sed(kBedRock),
  dgrade( num ), paleocurrent(-1.)
{
  if (0) //DEBUG
    std::cout << "tLayer( num )" << std::endl;
}

//copy constructor
inline tLayer::tLayer( const tLayer &orig ) :                        //tLayer
  ctime(orig.ctime), rtime(orig.rtime), etime(orig.etime),
  depth(orig.depth), erody(orig.erody), sed(orig.sed),
  dgrade( orig.dgrade ), paleocurrent(orig.paleocurrent)
{
  if (0) //DEBUG
    std::cout << "tLayer copy const\n";
}

inline const tLayer &tLayer::operator=( const tLayer &right )     //tLayer
{
  if (0) //DEBUG
    std::cout << "tLayer op=\n";

  if( &right != this )
    {
      dgrade = right.dgrade;
      ctime=right.ctime;
      rtime=right.rtime;
      etime=right.etime;
      depth=right.depth;
      erody=right.erody;
      sed=right.sed;
      paleocurrent = right.paleocurrent;

    }
  return *this;
}


inline void tLayer::setCtime( double tt )
{
  ctime = tt;
}

inline double tLayer::getCtime() const
{
  return ctime;
}

inline void tLayer::setRtime( double tt )
{
  rtime = tt;
}

inline double tLayer::getRtime() const
{
  return rtime;
}

inline void tLayer::setEtime( double tt )
{
  etime = tt;
}

inline void tLayer::addEtime( double tt )
{
  etime += tt;
}

inline double tLayer::getEtime() const
{
  return etime;
}


inline double tLayer::getDepth() const
{
  return depth;
}

inline void tLayer::setErody( double ero)
{
  erody = ero;
}

inline double tLayer::getErody() const
{
  return erody;
}

inline void tLayer::setSed( tSed_t rg)
{
  sed = rg;
}

inline tLayer::tSed_t tLayer::getSed() const
{
  return sed;
}

inline void tLayer::setDgradesize( size_t i )
{
  dgrade.setSize(i);
}

inline size_t tLayer::getDgradesize( ) const
{
  return dgrade.getSize();
}

inline void tLayer::addDgrade( size_t i, double size )
{
   assert(i<dgrade.getSize());

   if(dgrade[i] < 0.0){
     std::cout <<"before dgrade "<< i <<" < 0., = " << dgrade[i] << " depth= " << depth << '\n';
   }
   assert( dgrade[i]>=0.0 );

   dgrade[i]+=size;

   if(dgrade[i] < 0.0){
     std::cout <<"after dgrade "<< i <<" < 0., = " << dgrade[i] << " depth= " << depth
          << " size=" << size << '\n';
   }
   assert( dgrade[i]>=0.0 );
   depth+=size;
}

inline void tLayer::addPaleoCurrent(double dh, double compass)
{

 if(paleocurrent >= 0 && paleocurrent <= 360 && dh > 0.0){
   if(compass == -1){
     // do nothing, do not change the existing paleocurrent info in this layer
   }
   else if(compass >= 0 && compass <= 360){
     paleocurrent = (paleocurrent*depth + compass*dh)/ (depth + dh);   //weighted for thickness/contribution
   }
 }
 else if(paleocurrent==-1.){
   paleocurrent = compass;
 }

}

inline void tLayer::setPaleoCurrent( double val ){
  paleocurrent = val;
}

inline double tLayer::getPaleoCurrent() const
{
 return paleocurrent;
}

inline double tLayer::getDgrade( size_t i) const
{
  assert( i<dgrade.getSize() );
  return dgrade[i];
}

inline tArray< double > const &
tLayer::getDgrade( ) const
{
  return dgrade;
}

/** class tErode ***********************************************************/
/*class tErode
  {
  friend class tChannel;
  friend class tLNode;
  public:
  tErode();
  tErode( const tErode & );
  tErode( int, int );
  ~tErode();
  const tErode &operator=( const tErode & );
  private:

  double sedinput;
  double dz;
  double zp;
  double qs;
  double qsp;
  double qsin;
  double qsinp;
  int nsmpts;
  tArray< double > smooth;
  double tau;
  };*/


/** @class tMeander
 */
class tMeander
{
  friend class tChannel;
  friend class tLNode;
public:
  tMeander();
  tMeander( const tMeander & );
  tMeander( bool, double, double );
  ~tMeander();
  const tMeander &operator=( const tMeander & );
private:
  double newx, newy;
  double deltax, deltay; /* Displacements in x and y from meandering*/
  double zoldright;	/* right bed elevation */
  double zoldleft;	/* left bed elevation*/
  double bankrough; //bank roughness (lambda in meander.f) w/ units of length
  tArray< double > xyzd;
  bool reachmember;  /* Flag indicating node has been included in a reach*/
  bool meander;      /* flag indicating if the point meanders */
};

/** @class tBedrock
 */
class tBedrock
{
  friend class tLNode;
public:
  tBedrock();
  tBedrock( const tBedrock & );
  ~tBedrock();
  const tBedrock &operator=( const tBedrock & );
private:
  double erodibility;
};

/** class tSurface **********************************************************/
/*class tSurface
  {
  friend class tLNode;
  public:
  tSurface();
  tSurface( const tSurface & );
  ~tSurface();
  const tSurface &operator=( const tSurface & );
  private:
  double veg;          // Percent vegetation cover
  double tauc;         // Threshold
  double vegerody;     //erodibility of vegetated surface (or channel bank)
  };*/

/** @class tRegolith
 */
class tRegolith
{
  friend class tLNode;
public:
  tRegolith();
  //tRegolith( tInputFile &infile ); /* Reads needed values from input file*/
  tRegolith( const tRegolith & );
  ~tRegolith();
  const tRegolith &operator=( const tRegolith & );
private:
  double thickness;  /* dynamic thickness of regolith */
  tArray< double > dgrade;/* depth of each sediment class in active layer [m] */
};

/** @class tChannel
 */
class tChannel
{
  friend class tLNode;
public:
  tChannel();
  tChannel( const tChannel & );
  ~tChannel();
  const tChannel &operator=( const tChannel & );
private:
  double drarea;       /* drainage area (2/97)*/
  double q;  /* discharge in m^3/yr */
  double mdFlowPathLength;  /* Longest flow path from divide (9/01) */
  double chanwidth;    /* Channel geometry: width*/
  double hydrwidth;    /* hydraulic geometry: width*/
  double channrough;       /* Channel roughness (Manning 'n')*/
  double hydrnrough;       /* Hydraulic roughness (Manning 'n')*/
  double chandepth;    /* Channel flow depth*/
  double hydrdepth;    /* Hydraulic flow depth*/
  double chanslope;
  double hydrslope;
  double diam;    	/* Grain diameter of bed material*/
  /*member objects:*/
  tMeander migration;
};

/** @class tLNode
 */
class tLNode : public tNode
{
public:
  typedef enum {
    kFlooded     = 1,  // Flooding (lake) codes: part of a lake...
    kNotFlooded  = 0,  // ...or not...
    kCurrentLake = 2,  // ...or part one that is currently being computed...
    kSink        = 3,  // ...or a dry sink (unfilled depression).
    kOutletFlag  = 4,  // Used as temporary flag in FillLakes.
    kOutletPreFlag= 5  // ditto
  } tFlood_t;

  static const char* FloodName( tFlood_t f ){
    switch(f){
    case kFlooded:      return "1-Flooded";
    case kNotFlooded:   return "0-NotFlooded";
    case kCurrentLake:  return "2-CurrentLake";
    case kSink:         return "3-Sink";
    case kOutletFlag:   return "4-OutletFlag";
    case kOutletPreFlag:return "5-OutletPreFlag";
    }
    ::abort();
  }

  tLNode();
  tLNode( const tInputFile &infile );
  tLNode( const tLNode & );
  //Syntax for calling copy constructor
  //tLNode *newnode = new tLNode( *oldtLNode );
  ~tLNode();
  const tLNode &operator=( const tLNode & );
  inline tVegCover &getVegCover();
  inline const tBedrock &getRock() const;
  inline const tRegolith &getReg() const;
  inline const tChannel &getChan() const;
  inline tFlood_t getFloodStatus() const;
  inline void setFloodStatus( tFlood_t status );
  inline tEdge * getFlowEdg();
  inline tEdge const * getFlowEdg() const;
  inline void setFlowEdg( tEdge * );
  void setFlowEdgToZero() { flowedge = 0; }
  inline void setDrArea( double );
  inline void setFlowPathLength( double );
  inline double getFlowPathLength() const;
  inline void AddDrArea( double );
  inline void AddDischarge( double );
  inline tLNode * getDownstrmNbr();
  inline tLNode const * getDownstrmNbr() const;
  tLNode* getUpstrmNbr();
  inline double getQ() const;  // Gets total discharge from embedded chan obj
  // fluvial discharge is in now in m^3/YR
  inline double calcSlope();    // Computes and returns slope in flow direction
  inline double getDSlopeDt();
  inline bool Meanders() const;
  inline void setMeanderStatus( bool );
  inline void setHydrWidth( double );
  inline void setChanWidth( double );
  inline double getHydrWidth() const;
  inline double getChanWidth() const;
  inline void setHydrDepth( double );
  inline void setChanDepth( double );
  inline double getHydrDepth() const;
  inline double getChanDepth() const;
  inline void setHydrRough( double );
  inline void setChanRough( double );
  inline double getHydrRough() const;
  inline double getChanRough() const;
  inline void setHydrSlope( double );
  inline void setChanSlope( double );
  inline double getHydrSlope() const;
  inline double getChanSlope() const;
  double getDiam() const;
  inline void setBankRough( double );
  inline double getBankRough() const;
  inline double getDrArea() const;
  double getTotalLayerDepth() const;
  inline tArray< double > getZOld() const;
  inline tArray< double > getNew2DCoords() const;   //for chan.migration.newx, newy
  inline void setNew2DCoords( double, double );      //        "
  inline tArray< double > getNew3DCoords() const;   //        "
  inline tArray< double > getLatDisplace() const;  //for chan.migration.deltax, deltay
  inline void setLatDisplace( double, double );      //        "
  inline void addLatDisplace( double, double );      //        "
  inline void setRock( const tBedrock & );
  inline void setVegCover( const tLNode * );
  inline void setReg( const tRegolith & );
  inline void setChan( const tChannel & );
  inline double getSubSurfaceDischarge() const;
  inline void setDischarge( double );
  inline void setSubSurfaceDischarge( double );
  inline void setZOld( double, double );
  inline void RevertToOldCoords();
  virtual inline void UpdateCoords();
  double DistNew( tLNode const *, tLNode const * ) const;
  void ActivateSortTracer();
  void DeactivateSortTracer();
  void MoveSortTracerDownstream();
  void FlagDownhillNodes();
  inline void AddTracer();
  bool NoMoreTracers() const;
  void EroDep( double dz );
  inline void setAlluvThickness( double );
  inline double getAlluvThickness() const;
  inline tArray< double > const & getAlluvThicknessm( ) const;
  inline void setBedErody( double );
  inline double getBedErody() const;
  inline void setReachMember( bool );
  inline bool getReachMember() const;
  // NOTE - For the get and set functions which involve arrays of size numg,
  // the arrays go from 0 to (numg-1) and must be indexed in this manner
   //Qs functions for transport rate
  inline void setQs( double );
  inline void setQs( size_t, double );
  inline double getQs() const;
  inline double getQs( size_t );
  inline tArray< double > const & getQsm( ) const;
  inline void addQs( double );
  inline void addQs( size_t, double );
  void addQs( tArray< double > const &);
   //Qsin functions for incoming fluvial sediment flux
  inline void setQsin( double );
  inline void setQsin( size_t, double );
  inline void setQsin( tArray< double > const & );
  void setQsinErrorHandler( size_t ) const ATTRIBUTE_NORETURN;
  inline void addQsin( double );
  inline void addQsin( size_t, double );
  void addQsin( tArray< double > const &);
  inline double getQsin() const;
  inline double getQsin( size_t );
  inline tArray< double > const & getQsinm( ) const;
   //Qsdin functions for incoming diffusive sediment flux
  inline void setQsdin( double );
  inline void setQsdin( size_t, double );
  inline void setQsdin( tArray< double > const & );
  void setQsdinErrorHandler( size_t ) const ATTRIBUTE_NORETURN;
  inline void addQsdin( double );
  inline void addQsdin( size_t, double );
  void addQsdin( tArray< double > const &);
  inline double getQsdin() const;
  inline double getQsdin( size_t );
  inline tArray< double > const & getQsdinm( ) const;

  inline void setGrade( size_t, double ) const;
  inline double getGrade( size_t ) const;
  inline tArray< double > const & getGrade( ) const;
  inline void setXYZD( tArray< double > const &);
  inline tArray< double > const & getXYZD() const;
  double DistFromOldXY() const;
  inline bool OnBedrock() const;
  inline double getDzDt() const;
  inline void setDzDt( double );
  inline void addDrDt(double);
  inline double getDrDt() const;
  inline void setDrDt( double );
  inline double getTau() const;
  inline void setTau( double );
  inline double getTauCrit() const;
  inline void setTauCrit( double );
  inline void setUplift( double );
  inline double getUplift() const;
  inline size_t getNumg() const;
  inline void setNumg( size_t ) const;
  inline double getMaxregdep() const;
  // NOTE for the get and set functions which involve the layerlist
  // the top layer is layer 0 and indexes go from 0 to (getNumLayer-1)
  // NOTE The set functions for layer depths (including grades)
  // do not obey the rules of maximum layer depths
  // and should be only used for inititalizing values
  // Other than that, use EroDep() for erosion or deposition.
  // This function takes the layer index because there may be
  // erosion of the first few layers if the surface layer is not deep enough
  // The addtoLayer() function is a helper to addtoSurfaceDgrade()
  inline double getLayerCtime(int) const;
  inline double getLayerRtime(int) const;
  inline double getLayerEtime(int) const;
  inline double getLayerDepth(int) const;
  inline double getLayerErody(int) const;
  inline tLayer::tSed_t getLayerSed(int) const;
  inline double getLayerDgrade(int, size_t) const;  // first int is layer index
  // second int is grade index - see note above for indexing directions
  inline int getNumLayer() const;
  inline void setLayerCtime(int, double);
  inline void setLayerRtime(int, double);
  inline void setLayerEtime(int, double);
  inline void addLayerEtime(int, double);
  inline void setLayerDepth(int, double);
  inline void setLayerErody(int, double);
  inline void setLayerSed(int, tLayer::tSed_t);
  inline void setLayerDgrade(int, size_t, double);
  void setStratNode(tStratNode *s_) {
    stratNode = s_;
  }
  tStratNode *getStratNode() {
    return stratNode;
  }
  void FindInitFlowDir(); // was tStreamNet::
  bool FindFlowDir(); // was tStreamNet::
  bool FindDynamicFlowDir();
  tArray<double> EroDep(int, tArray<double>, double);
  // returns the depth of of each size that was actually deposited or
  // eroded.  Important in case less can be eroded than planned.
  // Can be used for erosion of bedrock.
  // Algorithm assumes that the material being deposited is the
  // Same material as that in the layer you are depositing into.

  // New functions below, added 28-10-2003, Quintijn
  inline void ResetAccummulatedDh();                           // set to 0.0
  inline void IncrementAccummulatedDh(const tArray<double>&);  // increment dh's, called during tFloodplain functions
  inline const tArray<double>& getAccummulatedDh() const;      // returns array with current values
  inline double getAccCoarse() const;	             	       // returns gravel thickness, accummulated during an other function
  inline double getAccFine() const;                            // returns fines thickness, accummulated in a
  // end new functions added 28-10-2003

  tArray<double> addtoLayer(int, double);
  // Used if removing material from lower layers -
  // only called from EroDep
  // because appropriate checking needs to be done first.
  // array tells the composition of the material which was taken from layer
  void addtoLayer(int, int, double, double);
  // Used if depositing or eroding material to lower layer size by size
  // only called from EroDep because appropriate checking needs
  // to be done first - also used for erosion from the surface layer
  void makeNewLayerBelow(int, tLayer::tSed_t, double, tArray<double> const &,
			 double);
  void removeLayer(int);
  void InsertLayerBack( tLayer const & );
  void LayerInterpolation( tTriangle const *, double, double, double );
  virtual void WarnSpokeLeaving(tEdge *);
  virtual void InitializeNode();
  virtual tArray< double > FuturePosn();
  virtual bool isMobile() const { return Meanders();}
  inline virtual bool flowThrough( tEdge const *e) const;
  virtual tNode *splitFlowEdge();
  inline virtual void setDownstrmNbr( tNode* ); // overrides tNode

  virtual void PrepForAddition( tTriangle const *, double );
  virtual void PrepForMovement( tTriangle const *, double );

  void CopyLayerList( tLNode const * ); // Copy layerlist from another node (gt 12/99)

   virtual void ChangeZ( double ); // overrides tNode::ChangeZ
#ifndef NDEBUG
  void TellAll() const;
#endif

protected:
  double CalcSlopeMeander(); // specialisation of CalcSlope()
  tLNode *getDSlopeDtMeander( double &curlen );  // specialisation of getDSlopeDt()
protected:
  tVegCover vegCover;  // Vegetation cover properties (see tVegetation.h/.cpp)
  tBedrock rock;
  tRegolith reg;
  tChannel chan;
  tFlood_t flood;        /* flag: is the node part of a lake?*/
  tEdge *flowedge;
  int tracer;                       /* Used by network sorting algorithm*/
  double dzdt;                      /* Erosion rate */
  double drdt;                      /* Rock erosion rate */
  double tau;                       // Shear stress or equivalent (e.g., unit stream pwr)
  //070427 double tauc;                      // Critical (threshold) shear stress or equiv.
  double taucb;                      // Critical (threshold) shear stress or equiv. for bedrock 070427
  double taucr;                      // For regolith
  // NOTE - all sediment transport rates are volume per year
  double qs;                        /* Sediment transport rate*/
  tArray< double > qsm;             /* multi size; transport rate of each size fraction*/
  double qsin;                      /* Sediment influx rate*/
  tArray< double > qsinm;           /* multi size; influx rate of each size fraction*/
  double qsdin;                      /* Sediment influx rate from diffusion*/
  tArray< double > qsdinm;           /* multi size; influx rate of each size fraction*/
  double uplift;                    /* uplift rate*/
  tList< tLayer > layerlist;        /* list of the different layers */
  tStratNode *stratNode;            /* Pointer to rectangular grid node. */
  tArray< double > accumdh;         /* temp inremental record of dh-fine and dh-coarse, used for projection on the tStratGrid */
  static size_t numg;
  // number of grain sizes recognized NIC should be the same for all
  // nodes, maybe put this somewhere else when you figure out what is going on
  static tArray< double > grade;
  // size of each grain size class, NIC again, you may
  // want to put this somewhere else
  static double maxregdep;
  static double KRnew;
  double qsubsurf;   // Subsurface discharge
public:
  int public1; // a "public" member that can be used for various purpose
};

inline tLNode::tFlood_t tLNode::getFloodStatus() const { return flood; }

inline void tLNode::setFloodStatus( tFlood_t status )
{
  flood = status;
}

inline tEdge * tLNode::getFlowEdg()
{
  return flowedge;
}

inline tEdge const * tLNode::getFlowEdg() const
{
  return flowedge;
}

inline void tLNode::setFlowEdg( tEdge * newflowedge )
{
  assert( newflowedge != 0 );  // Fails when passed an invalid edge
  flowedge = newflowedge;
}

inline void tLNode::setDrArea( double val ) {chan.drarea = val;}
inline void tLNode::AddDrArea( double val ) {chan.drarea += val;}

inline void tLNode::AddDischarge( double val ) {chan.q += val;}

inline tLNode * tLNode::getDownstrmNbr()
{
  return ( flowedge != 0 ) ?
    static_cast<tLNode *>(flowedge->getDestinationPtrNC()):
    0;
}

inline tLNode const * tLNode::getDownstrmNbr() const
{
  return ( flowedge != 0 ) ?
    static_cast<tLNode const *>(flowedge->getDestinationPtr()):
    0;
}

// nb: if channel is integrated into node, change this
inline double tLNode::getQ() const
{
  return chan.q;
}

/************************************************************************\
 **  calcSlope: Computes and returns the slope of the node's flowedg, or
 **  zero if the slope is less than zero.
 **
 **  The name of this function ("getSlope") makes NG cringe.
 **  Change to CALCSLOPE!!!! Done AD 09/2003
 **
 **  Assumptions: edge lengths up to date and nonzero, flowedg's up to
 **    date.
\************************************************************************/
inline double tLNode::calcSlope()
{
  assert( flowedge->getLength()>0 ); // failure means lengths not init'd

  const double slp =
    Meanders() ?
    CalcSlopeMeander():
    (z - getDownstrmNbr()->z ) / flowedge->getLength();

  //if( timetrack >= kBugTime ) std::cout << "GS 4; " << std::endl;
  return ( slp>=0.0 ) ? slp : 0.0;
}

inline double tLNode::getDSlopeDt()
{
  assert( flowedge != 0 );
  assert( flowedge->getLength()>0 ); // failure means lengths not init'd
  double curlen, slp;
  tLNode *dn;
  if( Meanders() )
    {
      dn = getDSlopeDtMeander( curlen );
    }
  else
    {
      curlen = flowedge->getLength();
      assert( curlen > 0.0 );
      dn = getDownstrmNbr();
    }
  slp = ( dzdt - dn->dzdt + uplift - dn->uplift ) / curlen;
  return ( slp>=0.0 ) ? slp : 0.0;
}

inline bool tLNode::Meanders() const {return chan.migration.meander;}
inline void tLNode::setMeanderStatus( bool val )
{ chan.migration.meander = val;}

inline double tLNode::getHydrWidth() const {return chan.hydrwidth;}
inline double tLNode::getChanWidth() const {return chan.chanwidth;}
inline double tLNode::getHydrDepth() const {return chan.hydrdepth;}
inline double tLNode::getChanDepth() const {return chan.chandepth;}
inline double tLNode::getHydrRough() const {return chan.hydrnrough;}
inline double tLNode::getChanRough() const {return chan.channrough;}
inline double tLNode::getHydrSlope() const {return chan.hydrslope;}
inline double tLNode::getChanSlope() const {return chan.chanslope;}
inline double tLNode::getSubSurfaceDischarge() const {return qsubsurf;}

inline double tLNode::getBankRough() const {return chan.migration.bankrough;}

//TODO: suggest doing away with the zero test for performance enhancement
inline void tLNode::setHydrWidth( double val )  {chan.hydrwidth = ( val > 0 ) ? val : 0;}
inline void tLNode::setChanWidth( double val )  {chan.chanwidth = ( val > 0 ) ? val : 0;}
inline void tLNode::setHydrDepth( double val )  {chan.hydrdepth = ( val > 0 ) ? val : 0;}
inline void tLNode::setChanDepth( double val )  {chan.chandepth = ( val > 0 ) ? val : 0;}
inline void tLNode::setHydrRough( double val )  {chan.hydrnrough = ( val > 0 ) ? val : 0;}
inline void tLNode::setChanRough( double val )  {chan.channrough = ( val > 0 ) ? val : 0;}
inline void tLNode::setHydrSlope( double val )  {chan.hydrslope = ( val > 0 ) ? val : 0;}
inline void tLNode::setChanSlope( double val )  {chan.chanslope = ( val > 0 ) ? val : 0;}
inline void tLNode::setBankRough( double val )
{chan.migration.bankrough = ( val > 0 ) ? val : 0;}
inline void tLNode::setSubSurfaceDischarge( double val ) {qsubsurf = val;}

inline double tLNode::getDrArea() const {return chan.drarea;}

inline tArray< double >
tLNode::getZOld() const
{
  return
    tArray< double >( chan.migration.zoldright,
		      chan.migration.zoldleft );
}

inline void tLNode::setZOld( double right, double left )
{
  chan.migration.zoldright = right;
  chan.migration.zoldleft = left;
}

inline tArray< double >
tLNode::getNew2DCoords() const
{
  return Meanders() ?
    tArray< double >( chan.migration.newx, chan.migration.newy ):
    tArray< double >( x, y );
}

inline tArray< double >                                                   //tNode
tLNode::getNew3DCoords() const
{
  return Meanders() ?
    tArray< double >( chan.migration.newx, chan.migration.newy, z ):
    tArray< double >( x, y, z );
}

inline void tLNode::setNew2DCoords( double val1, double val2 )
{
  chan.migration.newx = val1;
  chan.migration.newy = val2;
}

inline tArray< double > tLNode::
getLatDisplace() const
{
  return
    tArray< double >( chan.migration.deltax,
		      chan.migration.deltay );
}

inline void tLNode::setLatDisplace( double dx, double dy )
{
  chan.migration.deltax = dx;
  chan.migration.deltay = dy;
}

inline void tLNode::addLatDisplace( double dx, double dy )
{
  chan.migration.deltax += dx;
  chan.migration.deltay += dy;
}

inline void tLNode::setDischarge( double val ) {chan.q = ( val > 0 ) ? val : 0;}

inline void tLNode::RevertToOldCoords()
{
  chan.migration.newx = x;
  chan.migration.newy = y;
}

inline void tLNode::UpdateCoords()
{
  if (Meanders()) {
    x = chan.migration.newx;
    y = chan.migration.newy;
  }
}

inline bool tLNode::flowThrough( tEdge const *e) const {
  return
    e == flowedge;
}

/***********************************************************************\
  tLNode::setDownstrmNbr
  Set flowedge to flow to node "ptr".
  09/2003 AD/QC
  SL, 10/03: replaces new flowTo; I've been using this function
             in my version of CHILD since ca. 1998.
\***********************************************************************/
void tLNode::setDownstrmNbr( tNode *ptr ){
  // New node flows to the next meander node downstream (ptr)
  setFlowEdg( EdgToNod( ptr ) );
}

inline void tLNode::setBedErody( double val )
{rock.erodibility = ( val >= 0.0 ) ? val : 0.0;}

inline double tLNode::getBedErody() const {return rock.erodibility;}

inline void tLNode::setReachMember( bool val )
{chan.migration.reachmember = val;}

inline bool tLNode::getReachMember() const {return chan.migration.reachmember;}

inline void tLNode::setQs( double val ) {qs = val;}

inline void tLNode::setQs( size_t i, double val )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  qsm[i] = val;
  qs += val;
}

inline double tLNode::getQs() const {return qs;}

inline double tLNode::getQs( size_t i )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  return qsm[i];
}

inline tArray< double > const &
tLNode::getQsm( ) const
{
  return qsm;
}

inline void tLNode::addQs( double val )
{
  qs += val;
}

inline void tLNode::addQs( size_t i, double val )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  qsm[i] += val;
  qs += val;

}

inline void tLNode::setQsin( double val ) {qsin = val;}

inline void tLNode::setQsin( size_t i, double val )
{
  if(unlikely( (i>=numg) || (i>=qsinm.getSize()) )) {
    setQsinErrorHandler( i );
  }
  qsinm[i]=val;
  double tot=0.;
  for(size_t j=0; j<numg; j++)
    tot+=qsinm[j];
  qsin=tot;
}

inline void tLNode::setQsin( tArray< double > const & q_ )
{
  assert( qsinm.getSize() == q_.getSize() );

  double tot = 0.;
  for(size_t j=0; j<numg; j++) {
    qsinm[j] = q_[j];
    tot += q_[j];
  }
  qsin = tot;
}

inline void tLNode::setQsdin( double val ) {qsdin = val;}

inline void tLNode::setQsdin( size_t i, double val )
{
  if(unlikely( (i>=numg) || (i>=qsdinm.getSize()) )) {
    setQsdinErrorHandler( i );
  }
  qsdinm[i]=val;
  double tot=0.;
  for(size_t j=0; j<numg; j++)
    tot+=qsdinm[j];
  qsdin=tot;
}

inline void tLNode::setQsdin( tArray< double > const & q_ )
{
  assert( qsdinm.getSize() == q_.getSize() );

  double tot = 0.;
  for(size_t j=0; j<numg; j++) {
    qsdinm[j] = q_[j];
    tot += q_[j];
  }
  qsdin = tot;
}

inline void tLNode::addQsin( double val )
{
  qsin += val;
}

inline void tLNode::addQsin( size_t i, double val )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  qsinm[i] += val;
  qsin += val;

}

inline void tLNode::addQsdin( double val )
{
  qsdin += val;
}

inline void tLNode::addQsdin( size_t i, double val )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  qsdinm[i] += val;
  qsdin += val;

}

inline double tLNode::getQsin() const {return qsin;}

inline double tLNode::getQsin( size_t i )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  return qsinm[i];
}

inline tArray< double > const &
tLNode::getQsinm( ) const
{
  return qsinm;
}

inline double tLNode::getQsdin() const {return qsdin;}

inline double tLNode::getQsdin( size_t i )
{
  if(unlikely(i>=numg))
    ReportFatalError( "Trying to index sediment sizes that don't exist ");
  return qsdinm[i];
}

inline tArray< double > const &
tLNode::getQsdinm( ) const
{
  return qsdinm;
}

inline void tLNode::setGrade( size_t i, double size ) const
{
  if(unlikely(i>=numg))
    ReportFatalError("Trying to set a grain size for an index which is too large");
  grade[i] = size;
}

double tLNode::getGrade( size_t i ) const
{
  return grade[i];
}

tArray< double > const &
tLNode::getGrade( ) const
{
  return grade;
}

inline void tLNode::setXYZD( tArray< double > const &arr )
{
  chan.migration.xyzd = ( arr.getSize() == 4 ) ? arr : tArray< double >(4);
}

inline tArray< double > const &
tLNode::getXYZD() const {return chan.migration.xyzd;}

// Tests whether bedrock is exposed at a node
inline bool tLNode::OnBedrock() const
{
  // For multi-size model, criterion might be active layer thickness less
  // than a nominal thickness; here, it's just an arbitrary alluvial depth
  return ( reg.thickness<0.1 );
}

inline void tLNode::setDzDt( double val ) {dzdt = val;}

inline double tLNode::getDzDt() const {return dzdt;}

inline void tLNode::setDrDt( double val ) {drdt = val;}

inline void tLNode::addDrDt( double val )
{
  drdt += val;
}

inline double tLNode::getDrDt() const {return drdt;}

inline double tLNode::getTau() const { return tau; }

inline void tLNode::setTau( double newtau )
{
  tau = newtau;
}

// Modified to return either tauc for bedrock or regolith, whichever is top layer, Apr 07 gt
inline double tLNode::getTauCrit() const 
{
	if( getLayerSed(0)==0 )
		return taucb;
    else
		return taucr;
}

// Modified to set BOTH bedrock and regolith taucrit -- beware! Apr 07, gt
inline void tLNode::setTauCrit( double newtauc )
{
  taucb = newtauc;
  taucr = newtauc;
}

void tLNode::setUplift( double val ) {uplift = val;}

double tLNode::getUplift() const {return uplift;}

inline size_t tLNode::getNumg() const
{
  return numg;
}

inline void tLNode::setNumg( size_t size ) const
{
  numg = size;
}

inline double tLNode::getMaxregdep() const
{
  return maxregdep;
}

inline double tLNode::getLayerCtime( int i ) const
{
  return layerlist.getIthDataRef(i).getCtime();
}

inline double tLNode::getLayerRtime( int i ) const
{
  return layerlist.getIthDataRef(i).getRtime();
}

inline double tLNode::getLayerEtime( int i ) const
{
  return layerlist.getIthDataRef(i).getEtime();
}

inline double tLNode::getLayerDepth( int i ) const
{
  if( unlikely(layerlist.isEmpty()) )
    {
      std::cout << "** WARNING lyr list empty\n"
	   << " NODE " << id << ":\n"
	   << "  x=" << x << " y=" << y << " z=" << z;
    }
  return layerlist.getIthDataRef(i).getDepth();
}

inline double tLNode::getLayerErody( int i ) const
{
  return layerlist.getIthDataRef(i).getErody();
}

inline tLayer::tSed_t tLNode::getLayerSed( int i ) const
{
  return layerlist.getIthDataRef(i).getSed();
}

inline double tLNode::getLayerDgrade( int i, size_t num ) const
{
   return layerlist.getIthDataRef(i).getDgrade(num);
}

inline int tLNode::getNumLayer() const
{
  return layerlist.getSize();
}

inline void tLNode::setLayerCtime( int i, double tt)
{
  layerlist.getIthDataPtrNC( i )->setCtime( tt );
}

inline void tLNode::setLayerRtime( int i, double tt)
{
  layerlist.getIthDataPtrNC( i )->setRtime( tt );
}

inline void tLNode::setLayerEtime( int i, double tt)
{
  layerlist.getIthDataPtrNC( i )->setEtime( tt );
}

inline void tLNode::addLayerEtime( int i, double tt)
{
  layerlist.getIthDataPtrNC( i )->addEtime( tt );
}

inline void tLNode::setLayerDepth( int i, double dep)
{
  assert( dep > 0.0 );
  layerlist.getIthDataPtrNC( i )->setDepth( dep );
}

inline void tLNode::setLayerErody( int i, double ero)
{
  layerlist.getIthDataPtrNC( i )->setErody( ero );
}

inline void tLNode::setLayerSed( int i, tLayer::tSed_t s)
{
  layerlist.getIthDataPtrNC( i )->setSed( s );
}

inline void tLNode::setLayerDgrade( int i, size_t g, double val)
{
  assert( val>=0.0 );
  layerlist.getIthDataPtrNC( i )->setDgrade(g, val );
}

inline tVegCover & tLNode::getVegCover()
{
  return vegCover;
}

inline void tLNode::setVegCover( const tLNode *node )
{
  vegCover = node->vegCover;
}

inline double tLNode::getFlowPathLength() const
{
  return chan.mdFlowPathLength;
}

inline void tLNode::setFlowPathLength( double fpl )
{
  chan.mdFlowPathLength = fpl;
}

inline const tBedrock &tLNode::getRock() const {return rock;}
//Xconst tSurface &tLNode::getSurf() const {return surf;}
inline const tRegolith &tLNode::getReg() const {return reg;}
inline const tChannel &tLNode::getChan() const {return chan;}

inline void tLNode::setRock( const tBedrock & val ) {rock = val;}
//Xvoid tLNode::setSurf( const tSurface & val ) {surf = val;}
inline void tLNode::setReg( const tRegolith & val ) {reg = val;}
inline void tLNode::setChan( const tChannel & val ) {chan = val;}

// Functions related to the conversion of dh from tLNode to tStratNode (e.g. on the tSTRATGRID)
inline void tLNode::ResetAccummulatedDh()
{
  accumdh[0] = 0.;
  accumdh[1] = 0.;
}

inline void tLNode::IncrementAccummulatedDh(const tArray<double> &dh)
{
  accumdh[0]+= dh[0];
  accumdh[1]+= dh[1];
}

inline const tArray<double>& tLNode::getAccummulatedDh() const
{
  return accumdh;
}

inline double tLNode::getAccCoarse() const
{
  return accumdh[0];
}

inline double tLNode::getAccFine() const
{
  return accumdh[1];
}
// end of 5 functions related to the conversion of dh from tLNode to tStratNode

/**************************************************************************\
 **
 **  Tracer-sorting routines:
 **
 **  These routines are utilities that are used in sorting the nodes
 **  according to their position within the drainage network. The main
 **  sorting algorithm is implemented in tStreamNet::SortNodesByNetOrder().
 **  The sorting method works by introducing a "tracer" at each point,
 **  then allowing the tracers to iteratively cascade downstream. At each
 **  step any nodes not containing tracers are moved to the back of the
 **  list. The result is a list sorted in upstream-to-downstream order.
 **
 **  These utilities do the following:
 **    ActivateSortTracer -- injects a single tracer at a node
 **    AddTracer -- adds a tracer to a node (ignored if node is a bdy)
 **    MoveSortTracerDownstream -- removes a tracer and sends it to the
 **                                downstream neighbor (unless the node is
 **                                a sink; then the tracer just vanishes)
 **    FlagDownhillNodes -- for multiple flow direction
 **                         routing: flag all downhill nodes
 **    NoMoreTracers -- reports whether there are any tracers left here
 **
 **  Created by GT 12/97.
 **
 **  Modifications:
 **   - added MoveSortTracersDownstrmMulti and moved all files to .h
 **     for inlining, 1/00, GT
 **
\**************************************************************************/

inline void tLNode::ActivateSortTracer()
{ tracer = 1; }

inline void tLNode::DeactivateSortTracer()
{ tracer = 0; }

inline void tLNode::MoveSortTracerDownstream()
{
  tracer--;
  if( flood!=kSink ) getDownstrmNbr()->AddTracer();
}

inline void tLNode::FlagDownhillNodes()
{
  tEdge *ce;

  // Flag all adjacent nodes that are lower than me
  ce = getEdg();
  do
    {
      if( ce->getDestinationPtr()->getZ() < z && ce->FlowAllowed() )
	(static_cast<tLNode *>(ce->getDestinationPtrNC()))->ActivateSortTracer();
      ce = ce->getCCWEdg();
    }
  while( ce!=getEdg() );

}

inline void tLNode::AddTracer()
{
  if( boundary == kNonBoundary ) tracer++;
}

inline void tLNode::setAlluvThickness( double val )
{
  //reg.thickness = ( val >= 0.0 ) ? val : 0.0;
  // gt changed for performance speedup (if stmt shouldn't ever be needed)
  assert( val>=0 );
  reg.thickness = val;
}

inline double tLNode::getAlluvThickness() const {return reg.thickness;}

inline tArray< double > const &
tLNode::getAlluvThicknessm( ) const
{
  return reg.dgrade;
}

inline bool tLNode::NoMoreTracers() const
{
  assert( tracer>=0 );
  return BOOL( tracer==0 );
}

// SL, 9/2003: Version to override tNode version. ChangeZ is only used
//  by tUplift functions, so this version will cause uplift to be applied
//  to "old" elevations so that, if a long time elapses between setting
//  the coordinates for node dropping and the dropping itself, the new
//  node will be at an appropriate elevation.
inline void tLNode::ChangeZ( double val )
{
   tNode::ChangeZ( val );
   if( Meanders() ){
      if( chan.migration.xyzd[3] != 0.0 ){
         chan.migration.xyzd[2] += val;
      }
   }
}

#endif
