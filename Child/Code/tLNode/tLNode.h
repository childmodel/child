/************************************************************************\
**
**  tLNode.h
**
**  Header file for derived class tLNode and its member classes
**
**  $Id: tLNode.h,v 1.18 1998-04-06 19:27:53 nmgaspar Exp $
\************************************************************************/

#ifndef TLNODE_H
#define TLNODE_H
#include <iostream.h>
#include <fstream.h>

#include "../tArray/tArray.h"
#include "../GridElements/gridElements.h"
#include "../tList/tList.h"
#include "../tInputFile/tInputFile.h"

/** class tDeposit *********************************************************/
/* Deposit records */
class tDeposit
{
   friend class tListNode< tDeposit >;

  public:
   tDeposit();
   tDeposit( int );
   tDeposit( const tDeposit & );
   ~tDeposit();
   const tDeposit &operator=( const tDeposit & );

  private:
   tArray< double > dgrade; /* depth of each sediment class in deposit [m]
			       dgrade[0]=total depth of deposit */
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
   double thickness;
   int numal;          /* total number of alluvium deposits below active layer
                          does NOT count the active layer*/
   int numg;           /* number of grain sizes */
   tArray< double > dgrade;/* depth of each sediment class in active layer [m]
			     dgrade[0]= total depth of active layer */
   static tArray< double > grade;/* representive grain diameter of 
			     each sediment class [m], grade[0]=0 */
   double dpth;         /* dynamic depth of active layer */
   double actdpth;      /* standard active layer depth */
   tList< tDeposit > depositList;
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

/** class tLNode ***************************************************************/
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
    int GetFloodStatus();
    void SetFloodStatus( int status );
    tEdge * GetFlowEdg();
    void SetFlowEdg( tEdge * );
    void SetDrArea( double );
    void AddDrArea( double );
    tLNode * GetDownstrmNbr();
    double GetQ();        // Gets total discharge from embedded chan obj
    double GetSlope();    // Computes and returns slope in flow direction
    double GetDSlopeDt();
    int Meanders() const;
    void SetMeanderStatus( int );
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
    void setVegErody( double );
    double getVegErody() const;
    void setBedErody( double );
    double getBedErody() const;
    void setReachMember( int );
    int getReachMember() const;
    void setQs( double );
    double getQs() const;
    void setQsin( double );
    void AddQsin( double );
    double getQsin() const;
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
};

#endif
