/**************************************************************************\
**
**  gridElements.h: Header file for mesh (grid) elements tNode, tEdge,
**                  and tTriangle. Each of these grid elements is
**                  implemented as an object, as described below.
**
**  $Id: meshElements.h,v 1.9 1998-03-23 21:18:19 gtucker Exp $
**  (file consolidated from earlier separate tNode, tEdge, & tTriangle
**  files, 1/20/98 gt)
\**************************************************************************/

#ifndef GRIDELEMENTS_H
#define GRIDELEMENTS_H

#include <iostream.h>
#include "../Definitions.h"
#include "../tPtrList/tPtrList.h"
#include "../tArray/tArray.h"

class tEdge;

/**************************************************************************\
**  class tNode
**
**  tNode's are the points in a Delaunay triangulation, and their data
**  include x and y coordinates, a z value (which could be elevation or
**  some other variable), an ID, the point's Voronoi area and its reciprocal,
**  and a list of "spokes" (edges, of type tEdge) in the triangulation.  
**
**  Earlier modifications to tNode:
**  - GT changed getZ to type double, and defined the fn in .cpp,
**  11/16.
**  - added getX, getY functions
**  - added VoronoiArea function (was previously in tStreamNet)
**
\**************************************************************************/
class tNode
{
   friend ostream &operator<<( ostream &, tNode & );
   friend istream &operator>>( istream &, tNode & );
  public:
   tNode();
   tNode( const tNode & );
   ~tNode();
   const tNode &operator=( const tNode & );
   tArray< double > get3DCoords() const;
   tArray< double > get2DCoords() const;
   int getID() const;
   double getX() const;
   double getY() const;
   double getZ() const;
   double getVArea() const;
   double getVArea_Rcp() const;
   int getBoundaryFlag() const;
   const tPtrList< tEdge > &getSpokeList() const;
   tPtrList< tEdge > &getSpokeListNC();
   const tPtrListNode< tEdge > *getFirstSpokeNode() const;
   tPtrListNode< tEdge > *getFirstSpokeNodeNC();
   
   const tPtrListNode< tEdge > *
       getNextSpokeNode( const tPtrListNode< tEdge > * ) const;
   
   tPtrListNode< tEdge > *
       getNextSpokeNodeNC( tPtrListNode< tEdge > * ) const;
       
   void setID( int );
   void setX( double );
   void setY( double );
   void setZ( double );
   void ChangeZ( double );   // adds or subtracts from the current Z value
   void setVArea( double );
   void setVArea_Rcp( double );
   void set2DCoords( double, double );
   void set3DCoords( double, double, double );
   void insertFrontSpokeList( tEdge * );
   void insertBackSpokeList( tEdge * );
   void makeWheel();
   void setBoundaryFlag( int );
   tEdge * GetEdg();
   void SetEdg( tEdge * );
   double Dist( tNode *, tNode * );
   //X void CalcSpokeVEdgLengths();// TODO: delete; replaced by tEdge fn
   tEdge *EdgToNod( tNode * );
   double ComputeVoronoiArea();
   void makeCCWEdges();
   
  protected:
   int id;
   double x;
   double y;
   double z;
   double varea;        /* Voronoi area (2/97) */
   double varea_rcp;    /* Reciprocal of Voronoi area = 1/varea (for speed)*/
   int boundary;
   tEdge * edg;        // Ptr to one edge
   tPtrList< tEdge > spokeList; /* list of connected edges */
   const tEdge *NextSpoke( tPtrListNode< tEdge > * );
   const tEdge *NextSpokeNC( tPtrListNode< tEdge > * );
};



/***************************************************************************\
**
**  class tEdge
**
**  tEdge objects represent the directed edges in a Delaunay triangulation
**  of tNode objects. ("Directed" means that the edge has directionality
**  from one point to the other; one is the origin and the other the  
**  the destination). In addition to pointing to its origin and destination
**  nodes, each tEdge points to the tEdge that shares the same origin and
**  lies immediately counter-clockwise. This makes it possible to obtain,
**  given one tEdge, all of the tEdges connected to a given origin node.
**  Other data maintained by a tEdge includes its length, slope (if
**  applicable), a boundary flag, and the coordinates of the Voronoi
**  point associated with the left-side triangle.
**
**  Earlier modifications to tEdge:
**   - added data mbr slope and functions setSlope, getSlope,
**     getOrgZ, getDestZ, 11/16 gt
**   - definition of "rvpt" changed to "left-hand triangle" (name should
**     also be changed) gt 1/98
**
\***************************************************************************/

/** class tEdge ************************************************************/
class tEdge
{
   friend ostream &operator<<( ostream &, const tEdge & );
     /*friend istream &operator>>( istream &, tEdge & );*/
public:
    tEdge();
    tEdge( const tEdge & );
    ~tEdge();
    const tEdge &operator=( const tEdge & );
    int getID() const;
    int getBoundaryFlag() const;
    double getLength() const;
    double getSlope() const; // slope = "z" gradient from org to dest nodes 
    double getOrgZ();
    double getDestZ();
    const tNode *getOriginPtr() const;
    const tNode *getDestinationPtr() const;
    tNode *getOriginPtrNC();
    tNode *getDestinationPtrNC();
    int FlowAllowed();
    void setID( int );
    void setLength( double );
    void setSlope( double );
    void setOriginPtr( tNode * );
    void setDestinationPtr( tNode * );
    void setFlowAllowed( int );
    double CalcLength();
    double CalcSlope();
    tEdge * GetCCWEdg();
    void SetCCWEdg( tEdge * edg );
    tArray< double > getRVtx() const; // Get the Voronoi vertex for LH triangle
    void setRVtx( tArray< double > );
    double getVEdgLen() const;  // Get length of associated Voronoi cell edge
    void setVEdgLen( double );
    double CalcVEdgLen();   // Computes, sets & returns length of V cell edg
#ifndef NDEBUG
    void TellCoords();  // debug routine
#endif
    
  private:
   int id;
   int flowAllowed; // boundary flag, usu. false when org & dest = closed bds 
   double len;       // edge length
   double slope;
   tArray< double > rvtx; // (x,y) coords of Voronoi vertex in LH triangle
   double vedglen;        // length of Voronoi edge shared by org & dest
   tNode *org, *dest;    // ptrs to origin and destination tNodes
   tEdge *ccwedg;  // ptr to counter-clockwise (left-hand) edge w/ same origin 
};


/**************************************************************************\
**
**  class tTriangle
**
**  tTriangles are the Delaunay triangles that form the triangulated mesh.
**  Each tTriangle maintains pointers to its three nodes, three adjacent
**  triangles (if they exist), and to the three counter-clockwise directed
**  edges. For example, a tTriangle containing points a, b, and c (where
**  the order abc is counter-clockwise) would point to tNodes a,b,c and
**  tEdges a->b, b->c, and c->a, as well as to the tTriangles that share
**  sides ab, bc, and ac.
** 
**  Numbering convention: 
**   - points p0,p1,p2 are in counter-clockwise order
**   - adjacent triangle t0 lies opposite point p0
**   - directed edge e0 has points p0->p1
**
**  Modifications:
**   - added function FindCircumcenter(), 1/11/98 gt
**
\**************************************************************************/


/** class tTriangle *******************************************************/
class tTriangle
{
   friend ostream &operator<<( ostream &, const tTriangle & );
   friend istream &operator>>( istream &, tTriangle & );
public:
   tTriangle();
   tTriangle( const tTriangle & );
   ~tTriangle();
   /*overloaded operators*/
   const tTriangle &operator=( const tTriangle & );
   int getID() const;
   void setID( int );
   tNode *pPtr( int );
   tEdge *ePtr( int );
   tTriangle *tPtr( int );
   void setPPtr( int, tNode * );
   void setEPtr( int, tEdge * );
   void setTPtr( int, tTriangle * );
   int nVOp( tTriangle * );
   int nVtx( tNode * );
   tArray<double> FindCircumcenter();
#ifndef NDEBUG
   void TellAll();
#endif

private:
   int id;       /* Triangle ID # for testing*/
   tNode *p[3];
   tEdge *e[3];
   tTriangle *t[3];
};

#endif
