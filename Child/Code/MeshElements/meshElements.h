/**************************************************************************\
**
**  gridElements.h: Header file for mesh (grid) elements tNode, tEdge,
**                  and tTriangle. Each of these grid elements is
**                  implemented as an object, as described below.
**
**  This file contains declarations of the three classes that collectively
**  make up the triangulated mesh. These classes are:
**   - tNode: nodes (ie, the points in the triangulation)
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
**     tNode also maintains a list of all its spokes. In a future
**     version such lists will only be created when needed by mesh
**     modification routines (to reduce memory overhead).
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
**  $Id: meshElements.h,v 1.21 1999-02-04 21:56:51 gtucker Exp $
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
**  of pointers to all spokes (a "spoke list" of type tPtrList). In the
**  future, to conserve memory usage, these spoke lists will only be
**  allocated when needed by various mesh modification routines.
**
**  Other tNode variables include the area of the corresponding Voronoi
**  cell, an ID number, and a flag indicating the node's boundary status.
**  Possible boundary status codes are non-boundary (mesh interior), closed
**  (outer boundary, not counted as part of the solution domain), and
**  open (boundary node not part of the solution domain but representing
**  a valid exit point for mass or energy flows). Note that these boundary
**  codes are used by tGrid to segregate nodes according to whether they
**  are boundary or non-boundary points. Note also that while all hull
**  points must be boundaries, interior points do not necessarily have to
**  be flagged as non-boundaries (e.g., one could include an "island" of
**  boundary points in the interior of a mesh if needed).
**
**  Modifications:
**   - 2/99 GT added tNode::AttachNewSpoke and tEdge::WelcomeCCWNeighbor
**
\**************************************************************************/
class tNode
{
  friend ostream &operator<<( ostream &, tNode & );
  friend istream &operator>>( istream &, tNode & );

public:

  tNode();                                   // default constructor
  tNode( const tNode & );                    // copy constructor
  ~tNode();                                  // destructor

  const tNode &operator=( const tNode & );   // assignment operator
  tArray< double > get3DCoords() const;      // returns x,y,z
  tArray< double > get2DCoords() const;      // returns x,y
  int getID() const;                         // returns ID number
  double getX() const;                       // returns x coord
  double getY() const;                       // returns y coord
  double getZ() const;                       // returns z value
  double getVArea() const;                   // returns Voronoi area
  double getVArea_Rcp() const;               // returns 1/Voronoi area
  int getBoundaryFlag() const;               // returns boundary code
  tEdge * getEdg();                          // returns ptr to one spoke

  const tPtrList< tEdge > &getSpokeList() const;  // returns ref to spokelist
  tPtrList< tEdge > &getSpokeListNC();            // returns non-const ref "
  const tPtrListNode< tEdge > *getFirstSpokeNode() const; // returns 1st spoke
  tPtrListNode< tEdge > *getFirstSpokeNodeNC();   // returns non-const "
  void insertFrontSpokeList( tEdge * ); // adds edge ptr to front of spokelist
  void insertBackSpokeList( tEdge * );  // adds edge ptr to back of spokelist
  void makeWheel();                // makes spokelist circular
  void makeCCWEdges();             // sets up CCWEdg connectivity from spokelst
  
  // returns next spokelist element (as const and non-const, respectively)
  const tPtrListNode< tEdge > *
      getNextSpokeNode( const tPtrListNode< tEdge > * ) const;
  tPtrListNode< tEdge > *
      getNextSpokeNodeNC( tPtrListNode< tEdge > * ) const;
  
  void setID( int );              // sets ID number
  void setX( double );            // sets x coord
  void setY( double );            // sets y coord
  void setZ( double );            // sets z value
  void ChangeZ( double );         // adds or subtracts from the current z value
  void setVArea( double );        // sets Voronoi area
  void setVArea_Rcp( double );    // sets 1 / Voronoi area
  void set2DCoords( double, double );         // sets x and y values
  void set3DCoords( double, double, double ); // sets x, y, and z values
  void setBoundaryFlag( int );    // sets boundary status flag
  void setEdg( tEdge * );         // sets ptr to one spoke

  double Dist( tNode *, tNode * ); // distance from node to line (node1,node2)
  tEdge *EdgToNod( tNode * );      // finds spoke connected to given node
  double ComputeVoronoiArea();     // calculates node's Voronoi area
  void ConvertToClosedBoundary();  // makes node a closed bdy & updates edges
  void AttachFirstSpoke( tEdge * ); // welcomes first spoke (sets edg)
  virtual void WarnSpokeLeaving( tEdge *); // signals node that spoke is being deleted
   virtual void InitializeNode();  // used when new nodes are created, for now only has a purpose in inherited classes
   
#ifndef NDEBUG
   void TellAll();  // Debugging routine that outputs node data
#endif

   
protected:
  int id;           // ID number
  double x;         // x coordinate
  double y;         // y coordinate
  double z;         // z value (representing height or any other variable)
  double varea;     // Voronoi cell area
  double varea_rcp; // Reciprocal of Voronoi area = 1/varea (for speed)
  int boundary;     // Boundary status code
  tEdge * edg;      // Ptr to one edge
  tPtrList< tEdge > spokeList; // list of connected edges (spokes)

  const tEdge *NextSpoke( tPtrListNode< tEdge > * );  // returns -> next spoke
  const tEdge *NextSpokeNC( tPtrListNode< tEdge > * );// returns non-const "
};



/***************************************************************************\
**  class tEdge  ************************************************************
**
**  tEdge objects represent the directed edges in a Delaunay triangulation
**  of tNode objects. "Directed" means that the edge has directionality
**  from one point to the other; one is the origin and the other the  
**  the destination. In addition to pointing to its origin and destination
**  nodes, each tEdge points to the tEdge that shares the same origin and
**  lies immediately counter-clockwise. This makes it possible to obtain,
**  given one tEdge, all of the tEdges connected to a given origin node.
**  Other data maintained by a tEdge include its length, slope (if
**  applicable), a boundary flag, the coordinates of the Voronoi vertex
**  vertex associated with the right-hand triangle, and the length of
**  the corresponding Voronoi cell edge.
**
**  The boundary status of a tEdge object depends on the boundary status
**  of the two nodes to which it is connected: if either node is a closed
**  boundary, the edge's is a "no flow" (boundary) edge; otherwise it is a 
**  "flow allowed" (non-boundary) edge.
**
**  Note that an edge's slope is defined as the (Zo - Zd)/L, where Zo and
**  Zd are the z values of the origin and destination nodes, respectively,
**  and L is the edge's (projected) length.
**
**  Modifications:
**   - added FindComplement function, 4/98 GT
**
\***************************************************************************/
class tEdge
{
  friend ostream &operator<<( ostream &, const tEdge & );
     /*friend istream &operator>>( istream &, tEdge & );*/

public:

  tEdge();                // default constructor
  tEdge( const tEdge & ); // copy constructor
  ~tEdge();               // destructor

  const tEdge &operator=( const tEdge & );  // assignment operator
  int getID() const;            // returns ID number
  int getBoundaryFlag() const;  // returns boundary status (flow or no flow)
  double getLength() const;     // returns edge's length (projected)
  double getSlope() const;      // slope = "z" gradient from org to dest nodes 
  double getOrgZ();             // returns origin's z value
  double getDestZ();            // returns destination's z value
  const tNode *getOriginPtr() const;      // returns ptr to origin node (const)
  const tNode *getDestinationPtr() const; // returns ptr to dest node (const)
  tNode *getOriginPtrNC();      // returns ptr to origin node (non-const)
  tNode *getDestinationPtrNC(); // returns ptr to destination node (non-const)
  tEdge * getCCWEdg();          // returns ptr to counter-clockwise neighbor
  tArray< double > getRVtx() const;  // returns Voronoi vertex for RH triangle
  double getVEdgLen() const;    // returns length of assoc'd Voronoi cell edge
  int FlowAllowed();            // returns boundary status ("flow allowed")

  void setID( int );                 // sets ID number
  void setLength( double );          // sets edge length
  void setSlope( double );           // sets slope
  void setOriginPtr( tNode * );      // sets origin ptr
  void setDestinationPtr( tNode * ); // sets destination ptr
  void setFlowAllowed( int );        // sets boundary code
  double CalcLength();               // computes & sets length
  double CalcSlope();                // computes & sets slope
  void setCCWEdg( tEdge * edg );     // sets ptr to counter-clockwise neighbor
  void setRVtx( tArray< double > );  // sets coords of Voronoi vertex RH tri
  void setVEdgLen( double ); // sets length of corresponding Voronoi edge
  double CalcVEdgLen();      // computes, sets & returns length of V cell edg
  tEdge * FindComplement();  // returns ptr to edge's complement
  void WelcomeCCWNeighbor( tEdge * );  // Adds another edge ccw from this edge
  
#ifndef NDEBUG
  void TellCoords();  // debug routine that reports edge coordinates
#endif
  
private:
  int id;          // ID number
  int flowAllowed; // boundary flag, usu. false when org & dest = closed bds 
  double len;      // edge length
  double slope;    // edge slope
  tArray< double > rvtx; // (x,y) coords of Voronoi vertex in RH triangle
  double vedglen;        // length of Voronoi edge shared by org & dest cells
  tNode *org, *dest;     // ptrs to origin and destination nodes
  tEdge *ccwedg;   // ptr to counter-clockwise (left-hand) edge w/ same origin 
};


/**************************************************************************\
**  class tTriangle  *******************************************************
**
**  tTriangles are the Delaunay triangles that form the triangulated mesh.
**  Each tTriangle maintains pointers to its three nodes (vertices), three
**  adjacent triangles (if they exist), and to the three counter-clockwise
**  directed edges. For example, a tTriangle containing points a, b, and c
**  (where the order abc is counter-clockwise) would point to tNodes a,b,c
**  and tEdges a->b, b->c, and c->a, as well as to the tTriangles that share
**  sides ab, bc, and ac.
** 
**  Numbering convention: 
**   - points p0,p1,p2 are in counter-clockwise order
**   - adjacent triangle t0 lies opposite point p0
**   - directed edge e0 has points p0->p2
**
**  Modifications:
**   - added function FindCircumcenter(), 1/11/98 gt
**
\**************************************************************************/
class tTriangle
{
  friend ostream &operator<<( ostream &, const tTriangle & );
  friend istream &operator>>( istream &, tTriangle & );

public:
  tTriangle();                    // default constructor
  tTriangle( const tTriangle & ); // copy constructor
  ~tTriangle();                   // destructor

  const tTriangle &operator=( const tTriangle & ); // assignment operator
  int getID() const;                 // returns ID number
  tNode *pPtr( int );                // returns ptr to given vertex (0,1, or 2)
  tEdge *ePtr( int );                // returns ptr to given clockwise edge
  tTriangle *tPtr( int );            // returns ptr to given neighboring tri
  void setID( int );                 // sets ID number
  void setPPtr( int, tNode * );      // sets ptr to given vertex
  void setEPtr( int, tEdge * );      // sets ptr to given clockwise edge
  void setTPtr( int, tTriangle * );  // sets ptr to given neighboring tri
  int nVOp( tTriangle * );    // returns side # (0,1 or 2) of nbr triangle
  int nVtx( tNode * );        // returns vertex # (0,1 or 2) of given node
  tArray<double> FindCircumcenter(); // computes & returns tri's circumcenter

#ifndef NDEBUG
  void TellAll();  // debugging routine
#endif

private:
  int id;          // triangle ID number
  tNode *p[3];     // ptrs to 3 nodes (vertices)
  tEdge *e[3];     // ptrs to 3 clockwise-oriented edges
  tTriangle *t[3]; // ptrs to 3 neighboring triangles (or 0 if no nbr exists)
};

#endif
