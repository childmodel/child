/************************************************************************\
**
**  tLNode.h
**
**  Header file for derived class tLNode and its member classes
**
**  $Id: tLNode.h,v 1.25 1998-06-10 22:28:24 nmgaspar Exp $
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
   double getDepth() const;
   void setErody( double );
   double getErody() const;
   void setSed( int );
   int getSed() const;
   void setDgradesize( int );
   void setDgrade( int, double );
   double getDgrade( int );
   tArray< double > getDgrade() const;
   
  protected:
   double ctime; // time of creation of layer
   double rtime; // most recent time (time steps) that there was erosion/depo.
   double depth; // total depth of the layer
   double erody; // erodibility of layer (varies with material)
   int sed;  // 0 = bedrock; 1 = some mixture of sediments so there
   // may be multi-sizes.  although multiple sizes are stored for
   // bedrock, this is the texture of what the bedrock will break up
   // into, not what is there on the bed.
   // Later this may be used as a flag for alluvium vs. regolith, etc.
   tArray< double > dgrade; // depth of each size if multi size
   
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
    void setQs( double );
    void setQs( int, double );
    double getQs() const;
   double getQs( int );
    tArray< double > getQsm( ) const;
    void setQsin( double );
   void setQsin( int, double );
    void AddQsin( double );
   void AddQsin( int, double );
   void AddQsinm( tArray< double > );
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
    double getDrDt();
    void setDrDt( double );
    void setUplift( double );
    double getUplift() const;
   int getNumg() const;
   void setNumg( int );
   double getMaxregdepth() const;
   void setMaxregdepth( double );
   double getLayerCtime(int) const;
   double getLayerRtime(int) const;
   double getLayerDepth(int) const;
   double getLayerErody(int) const;
   int getLayerSed(int) const;
   double getLayerDgrade(int, int) const;
   int getNumLayer() const;
//    void setLayerCtime(int, double);
//    void setLayerRtime(int, double);
//    void setLayerDepth(int, double);
//    void setLayerErody(int, double);
//    void setLayerSed(int, int);
//    void setLayerDgrade(int, int, double);
    
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
    double qs;           /* Sediment transport rate*/
    tArray< double > qsm; /* multi size; transport rate of each size fraction*/
    double qsin;         /* Sediment influx rate*/
    tArray< double > qsinm; /* multi size; influx rate of each size fraction*/
   double uplift;  /*uplift rate*/
   tList< tLayer > layerlist; /* list of the different layers */
   //int numlay; /* number of layers in the list */
   // nic - not storing this for now - just use getSize() func from tList
   int numg; // number of grain sizes recognized NIC should be the same for all
   // nodes, maybe put this somewhere else when you figure out what is going on
   tArray< double > grade;  // size of each grain size class, again, you may
   // want to put this somewhere else
   double maxregdepth;  // maximum depth of a regolith layer
   // nic - for now the surface layer and max depth of layers below
   // will be the same
   
   
};

#endif
