/************************************************************************\
**
**  tLNode.h
**
**  Header file for derived class tLNode and its member classes
**
**  $Id: tLNode.h,v 1.2 1998-01-15 19:38:52 gtucker Exp $
\************************************************************************/

#ifndef TLNODE_H
#define TLNODE_H


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
   float dpth;   /* depth of this deposit , NOTE, this ignores porosity */
   tArray< float > dgrade;/*( NUMG );depth of that sediment class in deposit [m]*/
};

/** class tErode ***************************************************************/
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
   float sedinput;     /* Sed. volume input (output if neg) during an iteration*/
   float dz;           /* Elevation change during an iteration*/
   tArray< float > newdz;  /* for each sediment class */
   float totdz;        /* total dz = sum of dz of all sizes*/
   float zp;           /* Predicted elevation (used in numerical scheme)*/
   float qs;           /* Sediment transport rate*/
   float qsp;          /* Predicted sed trans rate at new time step*/
   float qsin;         /* Sediment influx rate*/
   float qsinp;        /* Predicted sed influx at new time step*/
   int nsmpts;         /* # of points downstream over which to apply smoothing*/
   tArray< float > smooth; /*weights for erosion applied to downstrm nodes*/
   float tau;          /* Shear stress (or similar quantity)*/
};


/** class tMeander *************************************************************/
class tMeander
{
   friend class tChannel;
   friend class tLNode;
  public:
   tMeander();
   tMeander( const tMeander & );
   tMeander( int, float, float );
   ~tMeander();
   const tMeander &operator=( const tMeander & );
  private:
   int meander;      /* flag indicating if the point meanders */
   float newx, newy;
   int head;         /* Flag indicating node is a reach head*/
   int reachmember;  /* Flag indicating node has been included in a reach*/
   float deltax, deltay; /* Displacements in x and y from meandering*/
   float zoldright;	/* right bed elevation */
   float zoldleft;	/* left bed elevation*/
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
   float erodibility;
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
   float veg;          /* Percent vegetation cover*/
   float tauc;         /* Threshold*/
     /*float vegresistfactor;*/
};

/** class tRegolith ************************************************************/
class tRegolith
{
   friend class tLNode;
  public:
   tRegolith();
   tRegolith( const tRegolith & );
   tRegolith( int, float ); /*number of grain sizes and active layer thickness*/
   ~tRegolith();
   const tRegolith &operator=( const tRegolith & );
  private:
   float thickness;
   int numal;          /* total number of alluvium deposits below active layer
                          does NOT count the active layer*/
   tArray< float > dgrade;/* depth of each sediment class in active layer [m]*/
   float dpth;         /* depth of active layer, changes but always returns
                          to contant amount at end of an iteration*/
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
   float drarea;       /* drainage area (2/97)*/
   float q;
   float chanwidth;    /* Channel geometry: width*/
   float hydrwidth;    /* hydraulic geometry: width*/
   float nrough;       /* Hydraulic roughness (Manning 'n')*/
   float chandepth;    /* Channel flow depth*/
   float hydrdepth;    /* Hydraulic flow depth*/
   float diam;    	/* Grain diameter of bed material*/
/*member objects:*/
   tErode erosion;
   tMeander migration;
};

/** class tLNode ***************************************************************/
class tLNode : public tNode
{
  public:
   tLNode();
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
   void SetDrArea( float );
   void AddDrArea( float );
   tLNode * GetDownstrmNbr();
   float GetQ();        // Gets total discharge from embedded chan obj
   float GetSlope();    // Computes and returns slope in flow direction
   int Meanders() const;
   void SetMeanderStatus( int );
   float getHydrWidth() const;
   float getChanWidth() const;
   float getDrArea() const;
   tArray< float > getNew2DCoords() const;
   void setNew2DCoords( float, float );
   tArray< float > getNew3DCoords() const;
   void setRock( const tBedrock & );
   void setSurf( const tSurface & );
   void setReg( const tRegolith & );
   void setChan( const tChannel & );
   void setDischarge( float );
   void RevertToOldCoords();
   void UpdateCoords();
   float DistNew( tLNode *, tLNode * );
   void ActivateSortTracer();
   void MoveSortTracerDownstream();
   void AddTracer();
   int NoMoreTracers();
   void EroDep( float dz );
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
   int tracer;       // Used by network sorting algorithm
};

#endif
