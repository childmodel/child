/**************************************************************************\
**
**  gridElements.h: Header file for mesh (grid) elements tNode, tEdge,
**                  and tTriangle. Each of these grid elements is
**                  implemented as an object, as described below.
**
**  This file contains declarations of the three classes that collectively
**  make up the triangulated mesh. These classes are:
**   - tNode: nodes, or vertices
**   - tEdge: directed edges, described by a starting node and an
**            ending node
**   - tTriangle: triangles in the mesh, with each triangle maintaining
**                pointers to its 3 vertex nodes, its 3 neighboring
**                triangles, and its 3 clockwise-oriented edges
**
**  Lists of each of these 3 types of mesh element are maintained by the
**  tGrid class, which implements the mesh and its routines. Connectivity
**  between mesh elements is managed using pointers, as follows:
**   - Each tNode object points to one of its "spokes" (the tEdges that
**     originate at the node). In the current implementation, each
**     tNode also maintains a list of all its spokes, but in a future
**     version such lists will only be created when needed by mesh
**     modification routines.
**   - Each tEdge points to its origin and destination nodes, and to
**     the tEdge that lies counterclockwise relative to the origin node.
**     tEdge objects also contain the coordinates of the the Voronoi
**     vertex that lies on the righthand side of the tEdge. (A Voronoi
**     vertex is the intersection between 3 Voronoi cells, and in a
**     Delaunay triangulation is found at the circumcenter of triangle).
**   - Each tTriangle object points to its 3 vertex nodes, its 3
**     neighboring triangles (or null if no neighboring triangle exists
**     across a given face), and the 3 clockwise-oriented tEdges. The
**     data structure uses the "opposite" numbering scheme, so that
**     triangle node 1 represents the vertex that is opposite to
**     neighboring triangle 1, and so on. Node 1 is also the origin for
**     edge 1, etc.
**   
**  $Id: meshElements.h,v 1.16 1999-01-11 22:54:04 gtucker Exp $
**  (file consolidated from earlier separate tNode, tEdge, & tTriangle
**  files, 1/20/98 gt)
\**************************************************************************/

#ifndef GRIDELEMENTS_H
#define GRIDELEMENTS_H

#include <iostream.h>
#include "../Definitions.h"
#include "../tPtrList/tPtrList.h"
#include "../tArray/tArray.h"
#include "../tGridList/tGridList.h"

class tEdge;

/**************************************************************************\
**  Class tNode  ***********************************************************
**
**  tNodes are the points in a Delaunay triangulation, and their data
**  include x and y coordinates, a z value (which could be elevation or
**  some other variable), an ID, the point's Voronoi area and its reciprocal,
**  and a pointer to one of its "spokes" (edges that originate at the node).
**  Because each spoke points to its counter-clockwise neighbor, tNode
**  objects only need to point to one of their spokes. However, for
**  convenience, the current implementation of tNode also contains a list
**  of pointers to all spokes (of type tPtrList). In the future, to
**  conserve memory usage, these spoke lists will only be allocated when
**  needed by various mesh modification routines.
**
**  Other tNode variables include the area of the corresponding Voronoi
**  cell, an ID number, and a flag indicating the node's boundary status.
**  Possible boundary status codes are non-boundary (mesh interior), closed
**  (outer boundary, not counted as part of the solution domain), and
**  open (boundary node not part of the solution domain but represents
**  a valid exit point for mass or energy flows). Note that these boundary
**  codes are used by tGrid to segregate nodes according to whether they
**  are boundary or non-boundary points. Note also that while all hull
**  points must be boundaries, interior points do not necessarily have to
**  be flagged as non-boundaries (e.g., one could include an "island" of
**  boundary points in the interior of a mesh without difficulty).
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
   tEdge * getEdg();
   void setEdg( tEdge * );
   double Dist( tNode *, tNode * );
   tEdge *EdgToNod( tNode * );
   double ComputeVoronoiArea();
   void makeCCWEdges();
   void ConvertToClosedBoundary();
   void WarnSpokeLeaving( tEdge *);
   
#ifndef NDEBUG
   void TellAll();  // Debugging routine
#endif
   
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
**  Subsequent modifications:
**   - added FindComplement function, 4/98 GT
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
    tEdge * getCCWEdg();
    void setCCWEdg( tEdge * edg );
    tArray< double > getRVtx() const; // get the Voronoi vertex for LH triangle
    void setRVtx( tArray< double > );
    double getVEdgLen() const;  // get length of associated Voronoi cell edge
    void setVEdgLen( double );
    double CalcVEdgLen();   /* Computes, sets & returns length of V cell edg */
    tEdge * FindComplement();
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
