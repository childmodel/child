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
**  $Id: tLNode.h,v 1.68 2003-08-05 13:14:09 childcvs Exp $
*/
/************************************************************************/

#ifndef TLNODE_H
#define TLNODE_H

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
# include <fstream>
using namespace std;
#else
# include <iostream.h>
# include <fstream.h>
#endif
#include <string.h>

#include "../tArray/tArray.h"
#include "../MeshElements/meshElements.h"
#include "../tList/tList.h"
#include "../tInputFile/tInputFile.h"
#include "../globalFns.h"
#include "../tVegetation/tVegetation.h"

#define kSink        3  // ...or a dry sink (unfilled depression).
#define kVeryHigh 100000  // Used in FillLakes


/** @class tLayer
   Layer records */
class tLayer
{
   friend class tListNode< tLayer >;

  public:
   tLayer();
   tLayer( int );
   tLayer( const tLayer & );
   const tLayer &operator=( const tLayer & );
   void setCtime( double );
   double getCtime() const;
   void setRtime( double );
   double getRtime() const;
   void setEtime( double );
   void addEtime( double );
   double getEtime() const;
   void setDepth( double );
   // NOTE setDepth will also update the depths in dgrade so that the
   // total depth is consistent and the original texture is kept.
   // If texture needs to also be changed, call setDgrade.  This will
   // also update depth and the texture will change appropriately.
   double getDepth() const;
   void setErody( double );
   double getErody() const;
   void setSed( int );
   int getSed() const;
   void setDgradesize( int );
   int getDgradesize() const;
   void setDgrade( int, double );
   double getDgrade( int );
   tArray< double > getDgrade() const;
   void addDgrade(int, double);
   
  protected:
   double ctime; // time of creation of layer
   double rtime; // most recent time (time steps) that there was erosion/depo.
   double etime; // exposure time, i.e. time material spent at surface
   double depth; // total depth of the layer
   double erody; // erodibility of layer (varies with material)
   int sed;  // 0 = bedrock; 1 = some mixture of sediments so there
   // may be multi-sizes.  although multiple sizes are stored for
   // bedrock, this is the texture of what the bedrock will break up
   // into, not what is there on the bed.
   // Later sed may be used as a flag for alluvium vs. regolith, etc.
   tArray< double > dgrade; // depth of each size if multi size
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
  depth(0.), erody(0.), sed(0),
  dgrade()
{
  if (0) //DEBUG
    cout << "tLayer( )" << endl;
}

inline tLayer::tLayer ( int num ) :
  ctime(0.), rtime(0.), etime(0.),
  depth(0.), erody(0.), sed(0),
  dgrade( num )
{
  if (0) //DEBUG
    cout << "tLayer( num )" << endl;
}

//copy constructor
inline tLayer::tLayer( const tLayer &orig ) :                        //tLayer
  ctime(orig.ctime), rtime(orig.rtime), etime(orig.etime),
  depth(orig.depth), erody(orig.erody), sed(orig.sed),
  dgrade( orig.dgrade )
{
  if (0) //DEBUG
    cout << "tLayer copy const\n";
}

inline const tLayer &tLayer::operator=( const tLayer &right )     //tLayer
{
   if (0) //DEBUG
     cout << "tLayer op=\n";
    
   if( &right != this )
   {
      dgrade = right.dgrade;
      ctime=right.ctime;
      rtime=right.rtime;
      etime=right.etime;
      depth=right.depth;
      erody=right.erody;
      sed=right.sed;
   
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

inline void tLayer::setSed( int rg)
{
   sed = rg;
}

inline int tLayer::getSed() const 
{
   return sed;
}

inline void tLayer::setDgradesize( int i )
{
   dgrade.setSize(i);
}

inline int tLayer::getDgradesize( ) const
{
   return dgrade.getSize();
}

inline void tLayer::addDgrade( int i, double size )
{
   assert(i<dgrade.getSize());
   dgrade[i]+=size;
   assert( dgrade[i]>=0.0 );
   depth+=size;
}

inline double tLayer::getDgrade( int i)
{
   assert( i<dgrade.getSize() );
   return dgrade[i];
}

inline tArray< double >
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
   bool meander;      /* flag indicating if the point meanders */
   double newx, newy;
   int head;         /* Flag indicating node is a reach head*/
   int reachmember;  /* Flag indicating node has been included in a reach*/
   double deltax, deltay; /* Displacements in x and y from meandering*/
   double zoldright;	/* right bed elevation */
   double zoldleft;	/* left bed elevation*/
   double bankrough; //bank roughness (lambda in meander.f) w/ units of length
   tArray< double > xyzd;
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
    tLNode();
    tLNode( tInputFile &infile );
    tLNode( const tLNode & );
   //Syntax for calling copy constructor
   //tLNode *newnode = new tLNode( *oldtLNode );
    ~tLNode();
    const tLNode &operator=( const tLNode & );
    inline tVegCover &getVegCover();
    const tBedrock &getRock() const;
    const tRegolith &getReg() const;
    const tChannel &getChan() const;
    inline int getFloodStatus() const;
    inline void setFloodStatus( int status );
    tEdge * getFlowEdg();
    void setFlowEdg( tEdge * );
    void setDrArea( double );
    inline void setFlowPathLength( double );
    inline double getFlowPathLength() const;
    void AddDrArea( double );
    void AddDischarge( double );
    inline tLNode * getDownstrmNbr();
    double getQ() const;  // Gets total discharge from embedded chan obj
   // fluvial discharge is in now in m^3/YR 
    double getSlope();    // Computes and returns slope in flow direction
    double getDSlopeDt();
    bool Meanders() const;
    void setMeanderStatus( bool );
    void setHydrWidth( double );
    void setChanWidth( double );
    double getHydrWidth() const;
    double getChanWidth() const;
    void setHydrDepth( double );
    void setChanDepth( double );
    double getHydrDepth() const;
    double getChanDepth() const;
    void setHydrRough( double );
    void setChanRough( double );
    double getHydrRough() const;
    double getChanRough() const;
    void setHydrSlope( double );
    void setChanSlope( double );
    double getHydrSlope() const;
    double getChanSlope() const;
    double getDiam() const;
    void setBankRough( double );
    double getBankRough() const;
    double getDrArea() const;
   double getTotalLayerDepth() const;
    tArray< double > getZOld() const;
    tArray< double > getNew2DCoords() const;   //for chan.migration.newx, newy
    void setNew2DCoords( double, double );      //        "
    tArray< double > getNew3DCoords() const;   //        "
    tArray< double > getLatDisplace() const;  //for chan.migration.deltax, deltay
    void setLatDisplace( double, double );      //        "
    void addLatDisplace( double, double );      //        "
    void setRock( const tBedrock & );
    inline void setVegCover( const tLNode * );
    void setReg( const tRegolith & );
    void setChan( const tChannel & );
    void setDischarge( double );
    void setZOld( double, double );
    void RevertToOldCoords();
    virtual void UpdateCoords();
    double DistNew( tLNode const *, tLNode const * ) const;
    void ActivateSortTracer();
    void DeactivateSortTracer();
    void MoveSortTracerDownstream();
    void FlagDownhillNodes();
    inline void AddTracer();
    int NoMoreTracers() const;
    void EroDep( double dz );
    void setAlluvThickness( double );
    double getAlluvThickness() const;
   tArray< double > getAlluvThicknessm( ) const;
    void setBedErody( double );
    double getBedErody() const;
    void setReachMember( int );
   int getReachMember() const;
// NOTE - For the get and set functions which involve arrays of size numg,
   // the arrays go from 0 to (numg-1) and must be indexed in this manner
   void setQs( double );
    void setQs( int, double );
    double getQs() const;
   double getQs( int );
    tArray< double > getQsm( ) const;
    void setQsin( double );
   void setQsin( int, double );
   void addQs( double );
   void addQs( int, double );
   void addQs( tArray< double > );
   void addQsin( double );
   void addQsin( int, double );
   void addQsin( tArray< double > );
    double getQsin() const;
   double getQsin( int );
   tArray< double > getQsinm( ) const;
   void setGrade( int, double ) const;
   double getGrade( int ) const;
   tArray< double > getGrade( ) const;
   void setXYZD( tArray< double > const &);
   tArray< double > getXYZD() const;
   double DistFromOldXY() const;
   int OnBedrock() const;
   double getDzDt() const;
   void setDzDt( double );
   void addDrDt(double);
   double getDrDt() const;
   void setDrDt( double );
   double getTau() const;
   inline void setTau( double );
   inline double getTauCrit() const;
   inline void setTauCrit( double );
   void setUplift( double );
   double getUplift() const;
   int getNumg() const;
   void setNumg( int ) const;
   double getMaxregdep() const;
   // NOTE for the get and set functions which involve the layerlist
   // the top layer is layer 0 and indexes go from 0 to (getNumLayer-1)
   // NOTE The set functions for layer depths (including grades)
   // do not obey the rules of maximum layer depths
   // and should be only used for inititalizing values
   // Other than that, use EroDep() for erosion or deposition.
   // This function takes the layer index because there may be
   // erosion of the first few layers if the surface layer is not deep enough
   // The addtoLayer() function is a helper to addtoSurfaceDgrade()
   double getLayerCtime(int) const;
   double getLayerRtime(int) const;
   double getLayerEtime(int) const;
   double getLayerDepth(int) const;
   double getLayerErody(int) const;
   int getLayerSed(int) const;
   double getLayerDgrade(int, int) const;  // first int is layer index
   // second int is grade index - see note above for indexing directions
   int getNumLayer() const;
   void setLayerCtime(int, double);
   void setLayerRtime(int, double);
   void setLayerEtime(int, double);
   void addLayerEtime(int, double);
   void setLayerDepth(int, double);
   void setLayerErody(int, double);
   void setLayerSed(int, int);
   void setLayerDgrade(int, int, double); 
   tArray<double> EroDep(int, tArray<double>, double);
   // returns the depth of of each size that was actually deposited or
   // eroded.  Important in case less can be eroded than planned.
   // Can be used for erosion of bedrock.
   // Algorithm assumes that the material being deposited is the
   // Same material as that in the layer you are depositing into.
   tArray<double> addtoLayer(int, double);
   // Used if removing material from lower layers -
   // only called from EroDep
   // because appropriate checking needs to be done first.
   // array tells the composition of the material which was taken from layer
   void addtoLayer(int, int, double, double);
   // Used if depositing or eroding material to lower layer size by size
   // only called from EroDep because appropriate checking needs
   // to be done first - also used for erosion from the surface layer
   void makeNewLayerBelow(int, int, double, tArray<double>, double);
   void removeLayer(int);
   void InsertLayerBack( tLayer const & );
   void LayerInterpolation( tTriangle const *, double, double, double );
   virtual void WarnSpokeLeaving(tEdge *);
   virtual void InitializeNode();
   virtual tArray< double > FuturePosn();
   virtual bool isMobile() const { return Meanders();}

   virtual void PrepForAddition( tTriangle const *, double );
   virtual void PrepForMovement( tTriangle const *, double );

   void CopyLayerList( tLNode const * ); // Copy layerlist from another node (gt 12/99)

#ifndef NDEBUG
   void TellAll();
#endif
   
protected:
   tVegCover vegCover;  // Vegetation cover properties (see tVegetation.h/.cpp)
   tBedrock rock;
   tRegolith reg;
   tChannel chan;
   int flood;        /* flag: is the node part of a lake?*/
   tEdge *flowedge;
   int tracer;       /* Used by network sorting algorithm*/
   double dzdt;      /* Erosion rate */
   double drdt;      /* Rock erosion rate */
   double tau;       // Shear stress or equivalent (e.g., unit stream pwr)
   double tauc;      // Critical (threshold) shear stress or equiv.
   // NOTE - all sediment transport rates are volume per year
   double qs;           /* Sediment transport rate*/
   tArray< double > qsm; /* multi size; transport rate of each size fraction*/
   double qsin;         /* Sediment influx rate*/
   tArray< double > qsinm; /* multi size; influx rate of each size fraction*/
   double uplift;  /*uplift rate*/
   tList< tLayer > layerlist; /* list of the different layers */
   static int numg;
   // number of grain sizes recognized NIC should be the same for all
   // nodes, maybe put this somewhere else when you figure out what is going on
   static tArray< double > grade;
   // size of each grain size class, NIC again, you may
   // want to put this somewhere else
   static double maxregdep;
   static double KRnew;
};

inline int tLNode::getFloodStatus() const { return flood; }

inline void tLNode::setFloodStatus( int status )
{
   flood = status;
}

inline tLNode * tLNode::getDownstrmNbr()
{
   return ( flowedge != 0 ) ?
     static_cast<tLNode *>(flowedge->getDestinationPtrNC()):
     0;
}

inline double tLNode::getTau() const { return tau; }

inline void tLNode::setTau( double newtau )
{
   tau = newtau;
}

inline double tLNode::getTauCrit() const { return tauc; }

inline void tLNode::setTauCrit( double newtauc )
{
   tauc = newtauc;
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
   if( !boundary ) tracer++;
}

inline int tLNode::NoMoreTracers() const 
{
   assert( tracer>=0 );
   return( tracer==0 );
}



#endif


