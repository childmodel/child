/************************************************************************\
**
**  tLNode.h
**
**  Header file for derived class tLNode and its member classes
**
**  $Id: tLNode.h,v 1.29 1998-07-15 22:27:36 nmgaspar Exp $
\************************************************************************/

#ifndef TLNODE_H
#define TLNODE_H
#include <iostream.h>
#include <fstream.h>

#include "../tArray/tArray.h"
#include "../GridElements/gridElements.h"
#include "../tList/tList.h"
#include "../tInputFile/tInputFile.h"
#include "../GlobalFns.h"
#include "../tRunTimer/tRunTimer.h"

/** class tLayer *********************************************************/
/* Layer records */
class tLayer
{
   friend class tListNode< tLayer >;

  public:
   tLayer();
   tLayer( int );
   tLayer( const tLayer & );
   ~tLayer();
   const tLayer &operator=( const tLayer & );
   void setCtime( double );
   double getCtime() const;
   void setRtime( double );
   double getRtime() const;
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
   void setFlag( int );
   int getFlag() const;   
   void setDgradesize( int );
   int getDgradesize();
   void setDgrade( int, double );
   double getDgrade( int );
   tArray< double > getDgrade() const;
   void addDgrade(int, double);
   
  protected:
   double ctime; // time of creation of layer
   double rtime; // most recent time (time steps) that there was erosion/depo.
   double depth; // total depth of the layer
   double erody; // erodibility of layer (varies with material)
   int sed;  // 0 = bedrock; 1 = some mixture of sediments so there
   // may be multi-sizes.  although multiple sizes are stored for
   // bedrock, this is the texture of what the bedrock will break up
   // into, not what is there on the bed.
   // Later sed may be used as a flag for alluvium vs. regolith, etc.
   tArray< double > dgrade; // depth of each size if multi size
   int flag; // put in to see if last change was erosion or deposition
   // 1=erosion, 2=deposition - testing purposes only
};

/** class tErode ***********************************************************/
class tErode
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
     /*int erodtype;*/
   double sedinput;     /* Sed. volume input (output if neg) during an iteration*/
   double dz;           /* Elevation change during an iteration*/
   double zp;           /* Predicted elevation (used in numerical scheme)*/
   double qs;           /* Sediment transport rate*/
   double qsp;          /* Predicted sed trans rate at new time step*/
   double qsin;         /* Sediment influx rate*/
   double qsinp;        /* Predicted sed influx at new time step*/
   int nsmpts;         /* # of points downstream over which to apply smoothing*/
   tArray< double > smooth; /*weights for erosion applied to downstrm nodes*/
   double tau;          /* Shear stress (or similar quantity)*/
};


/** class tMeander *************************************************************/
class tMeander
{
   friend class tChannel;
   friend class tLNode;
  public:
   tMeander();
   tMeander( const tMeander & );
   tMeander( int, double, double );
   ~tMeander();
   const tMeander &operator=( const tMeander & );
  private:
   int meander;      /* flag indicating if the point meanders */
   double newx, newy;
   int head;         /* Flag indicating node is a reach head*/
   int reachmember;  /* Flag indicating node has been included in a reach*/
   double deltax, deltay; /* Displacements in x and y from meandering*/
   double zoldright;	/* right bed elevation */
   double zoldleft;	/* left bed elevation*/
   double bankrough; //bank roughness (lambda in meander.f) w/ units of length
   tArray< double > xyzd;
};

/** class tBedrock *************************************************************/
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

/** class tSurface *************************************************************/
class tSurface
{
   friend class tLNode;
  public:
   tSurface();
   tSurface( const tSurface & );
   ~tSurface();
   const tSurface &operator=( const tSurface & );
  private:
   double veg;          /* Percent vegetation cover*/
   double tauc;         /* Threshold*/
   double vegerody;     //erodibility of vegetated surface (or channel bank)
};

/** class tRegolith ************************************************************/
class tRegolith
{
   friend class tLNode;
  public:
   tRegolith();
   tRegolith( tInputFile &infile ); /* Reads needed values from input file*/
   tRegolith( const tRegolith & );
   ~tRegolith();
   const tRegolith &operator=( const tRegolith & );
  private:
   double thickness;  /* dynamic thickness of regolith */
   tArray< double > dgrade;/* depth of each sediment class in active layer [m] */
};

/** class tChannel *************************************************************/
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
   double q;
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
   //tErode erosion;
   tMeander migration;
};

/** class tLNode *************************************************************/
class tLNode : public tNode
{
public:
    tLNode();
    tLNode( tInputFile &infile );
    tLNode( const tLNode & );
    ~tLNode();
    const tLNode &operator=( const tLNode & );   
    const tBedrock &getRock() const;
    const tSurface &getSurf() const;
    const tRegolith &getReg() const;
    const tChannel &getChan() const;
    int getFloodStatus();
    void setFloodStatus( int status );
    tEdge * getFlowEdg();
    void setFlowEdg( tEdge * );
    void setDrArea( double );
    void AddDrArea( double );
    void AddDischarge( double );
    tLNode * getDownstrmNbr();
    double getQ();        // Gets total discharge from embedded chan obj
   // fluvial discharge is in m^3/sec 
    double getSlope();    // Computes and returns slope in flow direction
    double getDSlopeDt();
    int Meanders() const;
    void setMeanderStatus( int );
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
    void setSurf( const tSurface & );
    void setReg( const tRegolith & );
    void setChan( const tChannel & );
    void setDischarge( double );
    void setDiam( double );
    void setZOld( double, double );
    void RevertToOldCoords();
    void UpdateCoords();
    double DistNew( tLNode *, tLNode * );
    void ActivateSortTracer();
    void MoveSortTracerDownstream();
    void AddTracer();
    int NoMoreTracers();
    void EroDep( double dz );
    void setAlluvThickness( double );
    double getAlluvThickness() const;
   tArray< double > getAlluvThicknessm( ) const;
    void setVegErody( double );
    double getVegErody() const;
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
   void setGrade( int, double );
   double getGrade( int );
   tArray< double > getGrade( ) const;
    void setXYZD( tArray< double > );
    tArray< double > getXYZD() const;
    double DistFromOldXY() const;
    int OnBedrock();
    double getDzDt();
    void setDzDt( double );
   void addDrDt(double);
    double getDrDt();
    void setDrDt( double );
    void setUplift( double );
    double getUplift() const;
   int getNumg() const;
   void setNumg( int );
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
   double getLayerDepth(int) const;
   double getLayerErody(int) const;
   int getLayerSed(int) const;
   int getLayerFlag(int) const;   
   double getLayerDgrade(int, int) const;  // first int is layer index
   // second int is grade index - see note above for indexing directions
   int getNumLayer() const;
   void setLayerCtime(int, double);
   void setLayerRtime(int, double);
   void setLayerDepth(int, double);
   void setLayerErody(int, double);
   void setLayerSed(int, int);
   void setLayerFlag(int, int);
   void setLayerDgrade(int, int, double); 
   tArray<double> EroDep(int, tArray<double>, double);
   // returns the depth of of each size that was actually deposited or
   // eroded.  Important in case less can be eroded than planned.
   // Can be used for erosion of bedrock.
   // Algorithm assumes that the material being deposited is the
   // Same material as that in the layer you are depositing into.
   tArray<double> addtoLayer(int, double, double);
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
   
#ifndef NDEBUG
   void TellAll();
#endif
   
protected:
   tBedrock rock;
   tSurface surf;
   tRegolith reg;
   tChannel chan;
   int flood;        /* flag: is the node part of a lake?*/
   tEdge *flowedge;
   int tracer;       /* Used by network sorting algorithm*/
   double dzdt;      /* Erosion rate */
   double drdt;      /* Rock erosion rate */
   // NOTE - all sediment transport rates are volume per year
   double qs;           /* Sediment transport rate*/
   tArray< double > qsm; /* multi size; transport rate of each size fraction*/
   double qsin;         /* Sediment influx rate*/
   tArray< double > qsinm; /* multi size; influx rate of each size fraction*/
   double uplift;  /*uplift rate*/
   tList< tLayer > layerlist; /* list of the different layers */
   //int numlay; /* number of layers in the list */
   // nic - not storing this for now - just use getSize() func from tList
   static int numg;
   // number of grain sizes recognized NIC should be the same for all
   // nodes, maybe put this somewhere else when you figure out what is going on
   static tArray< double > grade;
   // size of each grain size class, NIC again, you may
   // want to put this somewhere else
   static double maxregdep;
};

#endif


