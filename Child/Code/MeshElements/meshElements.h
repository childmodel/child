/**************************************************************************\
**
**  gridElements.h: Header file for mesh (grid) elements tNode, tEdge,
**                  and tTriangle. Each of these grid elements is
**                  implemented as an object, as described below.
**
**  $Id: meshElements.h,v 1.1 1998-01-21 01:10:26 gtucker Exp $
**  (file consolidated from earlier separate tNode, tEdge, & tTriangle
**  files, 1/20/98 gt)
\**************************************************************************/

#ifndef GRIDELEMENTS_H
#define GRIDELEMENTS_H

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
**  - GT changed getZ to type float, and defined the fn in .cpp,
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
   tArray< float > get3DCoords() const;
   tArray< float > get2DCoords() const;
   int getID() const;
   float getX() const;
   float getY() const;
   float getZ() const;
   float getVArea() const;
   float getVArea_Rcp() const;
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
   void setX( float );
   void setY( float );
   void setZ( float );
   void setVArea( float );
   void setVArea_Rcp( float );
   void set2DCoords( float, float );
   void set3DCoords( float, float, float );
   void insertFrontSpokeList( tEdge * );
   void insertBackSpokeList( tEdge * );
   void makeWheel();
   void setBoundaryFlag( int );
   tEdge * GetEdg();
   void SetEdg( tEdge * );
   float Dist( tNode *, tNode * );
   void CalcSpokeVEdgLengths();
   tEdge *EdgToNod( tNode * );
   float ComputeVoronoiArea();
   
  protected:
   int id;
   float x;
   float y;
   float z;
   float varea;        /* Voronoi area (2/97) */
   float varea_rcp;    /* Reciprocal of Voronoi area = 1/varea (for speed)*/
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
   float getLength() const;
   float getSlope() const; // slope is the "z" gradient from org to dest nodes 
   float getOrgZ();
   float getDestZ();
   const tNode *getOriginPtr() const;
   const tNode *getDestinationPtr() const;
   tNode *getOriginPtrNC();
   tNode *getDestinationPtrNC();
   int FlowAllowed();
   void setID( int );
   void setLength( float );
   void setSlope( float );
   void setOriginPtr( tNode * );
   void setDestinationPtr( tNode * );
   void setFlowAllowed( int );
   float CalcLength();
   tEdge * GetCCWEdg();
   void SetCCWEdg( tEdge * edg );
   tArray< float > getRVtx() const; // Get the Voronoi vertex for LH triangle
   void setRVtx( tArray< float > );
   float getVEdgLen() const;  // Get length of associated Voronoi cell edge
   void setVEdgLen( float );

  private:
   int id;
   int flowAllowed; // boundary flag, usu. false when org & dest = closed bds 
   float len;       // edge length
   float slope;
   tArray< float > rvtx; // (x,y) coords of Voronoi vertex in LH triangle
   float vedglen;        // length of Voronoi edge shared by org & dest
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
   tArray<float> FindCircumcenter();

private:
   int id;       /* Triangle ID # for testing*/
   tNode *p[3];
   tEdge *e[3];
   tTriangle *t[3];
};

#endif
