//-*-c++-*-

/*************************************************************************/\
/**
**  @file meshElements.h
**  @brief Header file for mesh elements tNode, tEdge,
**         and tTriangle. Each of these mesh elements is
**         implemented as an object, as described below.
**         (formerly called gridElements)
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
**  tMesh class, which implements the mesh and its routines. Connectivity
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
**  Significant modifications:
**   - 2/2/00: GT transferred get/set, constructors, and other small
**     functions from .cpp file to inline them
**
**  $Id: meshElements.h,v 1.68 2004-02-18 16:51:49 childcvs Exp $
**  (file consolidated from earlier separate tNode, tEdge, & tTriangle
**  files, 1/20/98 gt)
*/
/**************************************************************************/

#ifndef MESHELEMENTS_H
#define MESHELEMENTS_H

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include <math.h>       // for sqrt() used in inlined fn below
#include "../Definitions.h"
#include "../tPtrList/tPtrList.h"
#include "../tArray/tArray.h"
#include "../tMeshList/tMeshList.h"
#include "../Geometry/geometry.h"   // for Point2D definitions & fns

class tEdge;

/**************************************************************************/
/**
**  @class tNode
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
**  codes are used by tMesh to segregate nodes according to whether they
**  are boundary or non-boundary points. Note also that while all hull
**  points must be boundaries, interior points do not necessarily have to
**  be flagged as non-boundaries (e.g., one could include an "island" of
**  boundary points in the interior of a mesh if needed).
**
**  Modifications:
**   - 2/99 GT added tNode::AttachNewSpoke and tEdge::WelcomeCCWNeighbor
**   - 2/00 GT added getVoronoiVertexList and getVoronoiVertexXYZList.
**          These fns can be used in adding new nodes at V. Vertices.
**
*/
/**************************************************************************/
class tNode
{
  friend ostream &operator<<( ostream &, tNode & );
  friend istream &operator>>( istream &, tNode & );

public:

  tNode();                                   // default constructor
  tNode( const tNode & );                    // copy constructor
  virtual ~tNode() { edg = 0; }

  const tNode &operator=( const tNode & );   // assignment operator
  tArray< double > get3DCoords() const;      // returns x,y,z
  tArray< double > get2DCoords() const;      // returns x,y
// SL, 7/2003: added less costly versions that pass references
  void get3DCoords( tArray< double >& ) const;
  void get2DCoords( tArray< double >& ) const;
  inline int getID() const;                         // returns ID number
  double getX() const;                       // returns x coord
  double getY() const;                       // returns y coord
  double getZ() const;                       // returns z value
  double getVArea() const;                   // returns Voronoi area
  double getVArea_Rcp() const;               // returns 1/Voronoi area
  tBoundary_t getBoundaryFlag() const;               // returns boundary code
  tEdge * getEdg();                          // returns ptr to one spoke
  tEdge const * getEdg() const;              // returns ptr to one spoke
  void getVoronoiVertexList( tList<Point2D> * );  // Returns list of V vertices
  void getVoronoiVertexXYZList( tList<Point3D> * ); // As above plus interp z

  void setID( int );              // sets ID number
  void setX( double );            // sets x coord
  void setY( double );            // sets y coord
  void setZ( double );            // sets z value
  virtual void ChangeZ( double );         // adds or subtracts from the current z value
  void setVArea( double );        // sets Voronoi area
  void setVArea_Rcp( double );    // sets 1 / Voronoi area
  void set2DCoords( double, double );         // sets x and y values
  void set3DCoords( double, double, double ); // sets x, y, and z values
  void setBoundaryFlag( tBoundary_t );    // sets boundary status flag
  void setEdg( tEdge * );         // sets ptr to one spoke

  double Dist( tNode const *, tNode const * ) const; // distance from node to line (node1,node2)
  tEdge *EdgToNod( tNode const * );// finds spoke connected to given node
  double ComputeVoronoiArea();     // calculates node's Voronoi area
  void ConvertToClosedBoundary();  // makes node a closed bdy & updates edges
  virtual void WarnSpokeLeaving( tEdge *); // signals node that spoke is being deleted
  virtual void InitializeNode();  // used when new nodes are created, 
                                  // for now only has a purpose in inherited classes
  virtual tArray< double > FuturePosn();
  virtual void UpdateCoords() {}
  virtual bool isMobile() const { return false;}
  virtual bool flowThrough( tEdge const * ) const { return false; }
  virtual tNode *splitFlowEdge() { return 0; }
  virtual void setDownstrmNbr( tNode * ) {}

  virtual void PrepForAddition( tTriangle const *, double ) {}
  virtual void PrepForMovement( tTriangle const *, double ) {}

#ifndef NDEBUG
   void TellAll() const;  // Debugging routine that outputs node data
#endif


protected:
  int id;           // ID number
  double x;         // x coordinate
  double y;         // y coordinate
  double z;         // z value (representing height or any other variable)
  double varea;     // Voronoi cell area
  double varea_rcp; // Reciprocal of Voronoi area = 1/varea (for speed)
  tBoundary_t boundary;     // Boundary status code
private:
  tEdge * edg;      // Ptr to one edge
};



/***************************************************************************\
**  @class tEdge
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
**   - added cwedg and tri pointers with corresponding get and set functions,
**     2/99 SL
**   - added compedg pointer and corresponding get and set functions, 4/00 SL
**
\***************************************************************************/
class tEdge
{
  friend ostream &operator<<( ostream &, const tEdge & );
     /*friend istream &operator>>( istream &, tEdge & );*/

public:

  tEdge();                // default constructor
  tEdge( const tEdge & ); // copy constructor
  tEdge( tNode*, tNode* ); // makes edge between two nodes
  tEdge( int, tNode*, tNode* ); // makes edge between two nodes w/ given id

  //~tEdge();               // destructor

  const tEdge &operator=( const tEdge & );  // assignment operator
  void InitializeEdge( tNode*, tNode*, tNode const *, bool useFuturePosn = false );
  inline int getID() const;            // returns ID number
  tBoundary_t getBoundaryFlag() const; // returns boundary status (flow or no flow)
  double getLength() const;     // returns edge's length (projected)
  double getSlope() const;      // slope = "z" gradient from org to dest nodes
  double getOrgZ() const;       // returns origin's z value
  double getDestZ() const;      // returns destination's z value
  const tNode *getOriginPtr() const;      // returns ptr to origin node (const)
  const tNode *getDestinationPtr() const; // returns ptr to dest node (const)
  tNode *getOriginPtrNC();      // returns ptr to origin node (non-const)
  inline tNode *getDestinationPtrNC(); // returns ptr to destination node (non-const)
  inline tEdge * getCCWEdg();          // returns ptr to counter-clockwise neighbor
  inline tEdge * getCWEdg();
  inline tEdge* getComplementEdge();
  inline tEdge const * getComplementEdge() const;
  inline void setComplementEdge( tEdge* );
  inline tArray< double > const & getRVtx() const;  // returns Voronoi vertex for RH triangle
  inline void getRVtx( tArray< double >& ) const; // less costly ref-passing version
  inline double getVEdgLen() const;    // returns length of assoc'd Voronoi cell edge
  inline tEdgeBoundary_t FlowAllowed() const; // returns boundary status ("flow allowed")

  inline void setID( int );                 // sets ID number
  inline void setLength( double );          // sets edge length
  inline void setSlope( double );           // sets slope
  inline void setOriginPtr( tNode * );      // sets origin ptr
  inline void setDestinationPtr( tNode * ); // sets destination ptr
  static tEdgeBoundary_t isFlowAllowed( const tNode*, const tNode* );
  void setFlowAllowed( tEdgeBoundary_t );        // sets boundary code
  inline void setFlowAllowed( const tNode*, const tNode* ); // sets boundary code
  double CalcLength();               // computes & sets length
  double CalcSlope();                // computes & sets slope
  void setCCWEdg( tEdge * edg );     // sets ptr to counter-clockwise neighbor
  void setCWEdg( tEdge * edg );
  void setRVtx( tArray< double > const &);  // sets coords of Voronoi vertex RH tri
  void setVEdgLen( double ); // sets length of corresponding Voronoi edge
  double CalcVEdgLen();      // computes, sets & returns length of V cell edg
  tEdge * FindComplement();  // returns ptr to edge's complement
  tTriangle* TriWithEdgePtr();
  void setTri( tTriangle* );

  inline bool isFlippable() const;

#ifndef NDEBUG
  void TellCoords();  // debug routine that reports edge coordinates
#endif

private:
  int id;          // ID number
  tEdgeBoundary_t flowAllowed; // boundary flag, usu. false when org & dest = closed bds
  double len;      // edge length
  double slope;    // edge slope
  tArray< double > rvtx; // (x,y) coords of Voronoi vertex in RH triangle
  double vedglen;        // length of Voronoi edge shared by org & dest cells
  tNode *org, *dest;     // ptrs to origin and destination nodes
  tEdge *ccwedg;  // ptr to counter-clockwise (left-hand) edge w/ same origin
  tEdge *cwedg;   // ptr to clockwise (right-hand) edge w/ same origin
  tEdge *compedg; // ptr to complement edge
  tTriangle *tri; // ptr to triangle (if any) that contains pointer to edge;
                  // it's set from tTriangle::setEPtr( tEdge* ptr )
};


/**************************************************************************/
/**
**  @class tTriangle
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
**   - added index_ and associated functions 08/2003
**
*/
/**************************************************************************/
class tTriangle
{
  friend ostream &operator<<( ostream &, const tTriangle & );
  friend istream &operator>>( istream &, tTriangle & );

public:
  tTriangle();                    // default constructor
  tTriangle( const tTriangle & ); // copy constructor
  tTriangle( int, tNode*, tNode*, tNode* );
  tTriangle( int, tNode*, tNode*, tNode*, tEdge*, tEdge*, tEdge* );
  ~tTriangle();                   // destructor

  const tTriangle &operator=( const tTriangle & ); // assignment operator
  void InitializeTriangle( tNode*, tNode*, tNode* );
  inline int getID() const;          // returns ID number
  inline tNode *pPtr( int ) const;          // returns ptr to given vertex (0,1, or 2)
  inline tEdge *ePtr( int ) const;          // returns ptr to given clockwise edge
  inline tTriangle *tPtr( int ) const;      // returns ptr to given neighboring tri
  inline void setID( int );                 // sets ID number
  inline void setPPtr( int, tNode * );      // sets ptr to given vertex
  inline void setEPtr( int, tEdge * );      // sets ptr to given clockwise edge
  inline void setTPtr( int, tTriangle * );  // sets ptr to given neighboring tri
  inline int nVOp( const tTriangle * ) const;// returns side # (0,1 or 2) of nbr triangle
  int nVtx( const tNode * ) const;  // returns vertex # (0,1 or 2) of given node
  tArray<double> FindCircumcenter() const; // computes & returns tri's circumcenter
  const unsigned char *index() const { return index_; }
  void SetIndexIDOrdered(); // build the ordering index array
  inline bool isIndexIDOrdered() const;
  bool containsPoint(double, double) const; // does "this" contains the point (x,y)
  tTriangle* NbrToward( double, double );

#ifndef NDEBUG
  void TellAll() const;  // debugging routine
#endif

private:
  tNode *p[3];     // ptrs to 3 nodes (vertices)
  tEdge *e[3];     // ptrs to 3 clockwise-oriented edges
  tTriangle *t[3]; // ptrs to 3 neighboring triangles (or 0 if no nbr exists)
  int id;          // triangle ID number
  unsigned char index_[3]; // index used for ordered output

  inline void SetIndex(); // build the ordering index array in simple order
};


/**************************************************************************/
/**
**  @class tSpkIter
**
*/
/**************************************************************************/
class tSpkIter
{
public:
   tSpkIter();
   tSpkIter( tNode* );
   ~tSpkIter();

   void Reset( tNode* );

   tNode* CurNode();
   tEdge* CurSpoke();
   int getNumSpokes();

   int First();
   int Last();
   int Next();
   int Prev();
   int Get( int );
   int Get( const tEdge* );
   int Where() const;

   tEdge* FirstP();
   tEdge* LastP();
   tEdge* NextP();
   tEdge* PrevP();
   tEdge* GetP( int );
   tEdge* GetP( tEdge const * );
   tEdge* ReportNextP();
   tEdge* ReportPrevP();

   inline int AtEnd();
   inline int isEmpty();

   int insertAtPrev( tEdge* );
   int insertAtNext( tEdge* );
   inline int insertAtFront( tEdge* );
   int insertAtBack( tEdge* );
   tEdge* removePrev();
   tEdge* removeNext();
   tEdge* removeFromFront();
   tEdge* removeFromBack();

private:
   tNode* curnode;
   tEdge* curedg;
   int counter;
};

/***************************************************************************\
\**  Inlined Functions for class tNode  ************************************/

/***********************************************************************\
**
**  Constructors & destructors:
**
**  Default:  initializes values to zero.
**  Copy:  copies all values and makes duplicate spoke list
**  Destructor:  no longer used
**
\***********************************************************************/

//default constructor
inline tNode::tNode() :
  id(0),
  x(0.), y(0.), z(0.),
  varea(0.), varea_rcp(0.),
  boundary(kNonBoundary), edg(0)
{}

//copy constructor
inline tNode::tNode( const tNode &original ) :
  id(original.id),
  x(original.x), y(original.y), z(original.z),
  varea(original.varea), varea_rcp(original.varea_rcp),
  boundary(original.boundary), edg(original.edg)
{}

/*X tNode::~tNode()
{
   if (0)//DEBUG
     cout << "    ~tNode()" << endl;
}*/


/***********************************************************************\
**
**  Overloaded operators:
**
**    assignment: copies all values (spokelist's assignment operator
**                creates duplicate copy of list)
**    right shift: takes input for x, y and z values from input stream
**    left shift: sends the following data to the output stream:
**                node ID, x, y, z values, and IDs of neighboring nodes,
**                which are obtained through the spokelist.
**
\***********************************************************************/

//assignment
inline const tNode &tNode::operator=( const tNode &right )
{
   if( &right != this )
   {
      tPtrListIter< tEdge > spokIter;
      id = right.id;
      x = right.x;
      y = right.y;
      z = right.z;
      boundary = right.boundary;
      varea = right.varea;
      varea_rcp = right.varea_rcp;
      edg = right.edg;
   }
   return *this;
}

//right shift
inline istream &operator>>( istream &input, tNode &node )
{
   cout << "x y z:" << endl;
   input >> node.x >> node.y >> node.z;
   return input;
}

//left shift
inline ostream &operator<<( ostream &output, tNode &node )
{
   output << node.id << ": " << node.x << " " << node.y << " "
          << node.z
	  << endl;
   return output;
}


/***********************************************************************\
**
**  tNode "get" functions:
**
**  get3DCoords - returns x, y, z as a 3-element array
**  get2DCoords - returns x & y coords as a 2-element array
**  getID - returns ID #
**  getX - returns node's x coord
**  getY - returns node's y coord
**  getZ - returns node's z value
**  getVArea - returns Voronoi area
**  getVArea_Rcp - returns 1 / Voronoi area
**  getBoundaryFlag - returns boundary code
**  getEdg - returns pointer to one spoke
**
\***********************************************************************/

inline tArray< double >
tNode::get3DCoords() const
{
  return tArray< double > (x, y, z);
}

inline void tNode::get3DCoords( tArray< double >& xyz ) const
{
   if( xyz.getSize() != 3 ) xyz.setSize(3);
   xyz[0] = x;
   xyz[1] = y;
   xyz[2] = z;
}

inline tArray< double >
tNode::get2DCoords() const
{
  return tArray< double > (x, y);
}

inline void tNode::get2DCoords( tArray< double >& xy ) const
{
   if( xy.getSize() != 2 ) xy.setSize(2);
   xy[0] = x;
   xy[1] = y;
}

inline int tNode::getID() const {return id;}
inline double tNode::getX() const {return x;}
inline double tNode::getY() const {return y;}
inline double tNode::getZ() const {return z;}
inline double tNode::getVArea() const {return varea;}
inline double tNode::getVArea_Rcp() const {return varea_rcp;}
inline tBoundary_t tNode::getBoundaryFlag() const {return boundary;}
inline tEdge * tNode::getEdg() {return edg;}
inline tEdge const * tNode::getEdg() const {return edg;}

/***********************************************************************\
**
**  tNode "set" functions:
**
**  setID - sets ID number to val
**  setX - sets x coord to val
**  setY - sets y coord to val
**  setZ - sets z value to val
**  setVArea - sets Voronoi area to val
**  setVArea_Rcp - sets 1/Voronoi area to val
**  setBoundaryFlag - returns boundary code
**  set3DCoords - sets x, y, z to val1, val2, val3
**  set2DCoords - sets x & y to val1 and val2
**  setEdg - sets edge ptr to theEdg
**
**  Note: unless otherwise noted, no runtime value checking is done
**        (aside from assert statements)
**
\***********************************************************************/

inline void tNode::setID( int val ) {id = val;}
inline void tNode::setX( double val ) {x = val;}
inline void tNode::setY( double val ) {y = val;}
inline void tNode::setZ( double val ) {z = val;}

inline void tNode::setVArea( double val )
{
  assert( val>=0.0 );
  varea = val;
  /*varea = ( val >= 0.0 ) ? val : 0.0;*/
}

inline void tNode::setVArea_Rcp( double val )
{
  assert( val>=0.0 );
  varea_rcp = val;
  /*varea_rcp =  ( val >= 0.0 ) ? val : 0.0;*/
}

inline void tNode::setBoundaryFlag( tBoundary_t val )
{
  boundary = val;
}

inline void tNode::set2DCoords( double val1, double val2 )
{
   setX( val1 );
   setY( val2 );
}

inline void tNode::set3DCoords( double val1, double val2, double val3 )
{
   setX( val1 );
   setY( val2 );
   setZ( val3 );
}

inline void tNode::setEdg( tEdge * theEdg )
{
   edg = theEdg;
   if (0)//DEBUG
     cout << "Assigning edge " << theEdg->getID()
	  << " to node " << getID() << endl;
}


/***********************************************************************\
**
**  tNode::ChangeZ:  Adds delz to current z value
**
\***********************************************************************/
inline void tNode::ChangeZ( double delz ) { z += delz; }

/*******************************************************************\
**
**  tNode::WarnSpokeLeaving( tEdge * edglvingptr )
**
**  This function is called when an edge is being removed from the edge list.
**  If edg (the edge pointer member of tNode) is pointing to the edge
**  which will be removed, this edg must be updated.
**
**  edglvingptr is as it says, a pointer to the edge which will be
**  removed.
**
**  Called from tMesh::ExtricateEdge
**
**  9/98 NG and GT
\*******************************************************************/
inline void tNode::WarnSpokeLeaving( tEdge * edglvingptr )
{
  assert(edg);
  if( edglvingptr == edg )
    edg = edg->getCCWEdg();
}

/**********************************************************************\
 **
 **  tNode::InitializeNode()
 **
 **  A virtual function.
 **  This functions doesn't do anything here, only in inherited classes.
 **  Used for initializing things in newly created nodes that are set up
 **  for the rest of the nodes when the mesh is created.
 **
 **  1/1999  NG
 \**********************************************************************/
inline void tNode::InitializeNode()
{
}

/*******************************************************************\
  tNode::FuturePosn() virtual function; here, just calls get2DCoords()

  11/98 SL
  5/2003 AD
\*******************************************************************/
inline tArray< double > tNode::FuturePosn() {return get2DCoords();}



/**************************************************************************\
\***  Functions for class tEdge  ******************************************/


/***********************************************************************\
**
**  Constructors & destructors:
**
**  Default:  initializes values to zero and makes rvtx a 2-elem array
**  Copy:  copies all values
**  Destructor:  no longer used
**
\***********************************************************************/

//default constructor
inline tEdge::tEdge() :
  id(0), flowAllowed(kFlowNotAllowed), len(0.), slope(0.),
  rvtx(2),
  vedglen(0.),
  org(0), dest(0), ccwedg(0), cwedg(0),
  compedg(0), tri(0)
{
   if (0)//DEBUG
     cout << "tEdge()" << endl;
}

//copy constructor
inline tEdge::tEdge( const tEdge &original ) :
  id(original.id), flowAllowed(original.flowAllowed),
  len(original.len), slope(original.slope),
  rvtx(original.rvtx),
  vedglen(original.vedglen),
  org(original.org), dest(original.dest), ccwedg(original.ccwedg),
  cwedg(original.cwedg), compedg(original.compedg), tri(original.tri)
{}

inline tEdge::tEdge(tNode* n1, tNode* n2) :
  id(0), flowAllowed(kFlowNotAllowed), len(0.), slope(0.),
  rvtx(2),
  vedglen(0.),
  org(0), dest(0), ccwedg(0), cwedg(0),
  compedg(0), tri(0)
{
  setOriginPtr( n1 );
  setDestinationPtr( n2 );
  setFlowAllowed( n1, n2 );
}

inline tEdge::tEdge(int id_, tNode* n1, tNode* n2) :
  id(id_), flowAllowed(kFlowNotAllowed), len(0.), slope(0.),
  rvtx(2),
  vedglen(0.),
  org(0), dest(0), ccwedg(0), cwedg(0),
  compedg(0), tri(0)
{
  setOriginPtr( n1 );
  setDestinationPtr( n2 );
  setFlowAllowed( n1, n2 );
}

//tEdge::~tEdge() {/*cout << "    ~tEdge()" << endl;*/}


/***********************************************************************\
**
**  Overloaded operators:
**
**    assignment: copies all values (spokelist's assignment operator
**                creates duplicate copy of list)
**    left shift: sends the following data to the output stream:
**                edge ID, length, slope, and origin and destination IDs
**
\***********************************************************************/
inline const tEdge &tEdge::operator=( const tEdge &original )
{
   if( &original != this )
   {
      id = original.id;
      len = original.len;
      slope = original.slope;
      rvtx = original.rvtx;
      vedglen = original.vedglen;
      org = original.org;
      dest = original.dest;
      ccwedg = original.ccwedg;
      cwedg = original.cwedg;
      flowAllowed = original.flowAllowed;
      compedg = original.compedg;
      tri = original.tri;
   }
   return *this;
}

//left shift
inline ostream &operator<<( ostream &output, const tEdge &edge )
{
   output << edge.id << " " << edge.len << " " << edge.slope << " "
          << edge.org->getID()
          << " " << edge.dest->getID() << endl;
   return output;
}

/***********************************************************************\
**
**  tEdge "get" functions:
**
**  getID - returns ID #
**  getBoundaryFlag - returns boundary code
**  getLength - returns projectd length
**  getSlope - returns slope
**  getOriginPtr - returns const ptr to origin node
**  getDestinationPtr - returns const ptr to destination node
**  getOriginPtrNC - returns non-const ptr to origin node
**  getDestinationPtrNC - returns non-const ptr to destination node
**  getOrgZ- returns z value of origin node
**  getDestZ - returns z value of destination node
**  getCCWEdg - returns ptr to counterclockwise neighboring edge
**  FlowAllowed - returns the boundary flag, which indicates whether
**                or not the edge is an active flow conduit (which is
**                true as long as neither endpoint is a closed bdy node)
**  getRVtx - returns coordinates of right-hand Voronoi vertex as a
**            2-element array
**  getVEdgLen - returns the length of the corresponding Voronoi edge
**
\***********************************************************************/

inline int tEdge::getID() const {return id;}

//return 0 if flow allowed to match kNonBoundary:
inline tBoundary_t tEdge::getBoundaryFlag() const
{return ( flowAllowed == kFlowAllowed )?kNonBoundary:kClosedBoundary; }

inline double tEdge::getLength() const {return len;}

inline double tEdge::getSlope() const {return slope;}

inline const tNode *tEdge::getOriginPtr() const {return org;}

inline const tNode *tEdge::getDestinationPtr() const {return dest;}

inline tNode *tEdge::getOriginPtrNC() {return org;}

inline tNode *tEdge::getDestinationPtrNC() {return dest;}

inline double tEdge::getOrgZ() const
{
   assert( org!=0 );
   return( org->getZ() );
}

inline double tEdge::getDestZ() const
{
   assert( dest!=0 );
   return( dest->getZ() );
}

inline tEdge * tEdge::getCCWEdg()
{
   return ccwedg;
}

inline tEdge * tEdge::getCWEdg()
{
  return cwedg;
}

inline tEdge* tEdge::getComplementEdge()
{
  return compedg;
}

inline tEdge const* tEdge::getComplementEdge() const
{
  return compedg;
}

inline void tEdge::setComplementEdge( tEdge* edg )
{
  compedg = edg;
}

inline tEdgeBoundary_t tEdge::FlowAllowed() const
{
   return flowAllowed;
}

inline tArray< double > const &
tEdge::getRVtx() const
{
   return rvtx;
}

inline void tEdge::getRVtx( tArray< double >& arr ) const
{
   if( arr.getSize() != 2 ) arr.setSize(2);
   arr[0] = rvtx[0];
   arr[1] = rvtx[1];
}

inline double tEdge::getVEdgLen() const {return vedglen;}

inline tTriangle* tEdge::TriWithEdgePtr() {return tri;}

/***********************************************************************\
**
**  tEdge "set" functions:
**
**  setID - sets ID # to val
**  setLength - sets length to val
**  setSlope - sets slope to slp
**  setOriginPtr - sets origin pointer to ptr (if ptr is nonzero)
**  setDestinationPtr - sets destination pointer to ptr (if nonzero)
**  setFlowAllowed - sets flowAllowed status to val
**  setCCWEdg - sets ptr to counter-clockwise neighbor to edg
**  setRVtx - sets the coordinates of the right-hand Voronoi vertex
**            (ie, the Voronoi vertex at the circumcenter of the RH
**            triangle) to the 1st two elements in arr, which is
**            assumed to be a 2-element array
**  setVEdgLen - sets vedglen to val (vedglen is the length of the
**               corresponding Voronoi cell edge)
**
**  Note: unless otherwise noted, no checking of range or validity is
**        performed in these routines (aside from assert statements)
**
\***********************************************************************/

inline void tEdge::setID( int val ) {
   assert( id>=0 );
   id = val;
}

inline void tEdge::setLength( double val )
{
   assert( val>=0.0 );
   len = val;
}

inline void tEdge::setSlope( double slp )
{ slope = slp; }

inline void tEdge::setOriginPtr( tNode * ptr ) {if( ptr != 0 ) org = ptr;}

inline void tEdge::setDestinationPtr( tNode * ptr )
{if( ptr != 0 ) dest = ptr;}

inline void tEdge::setFlowAllowed( tEdgeBoundary_t val )
{
   assert( val==kFlowAllowed || val==kFlowNotAllowed );
   flowAllowed = val;
}

inline tEdgeBoundary_t tEdge::isFlowAllowed( const tNode* n1, const tNode* n2 )
{
   assert( n1 && n2 );
   return ( n1->getBoundaryFlag() != kClosedBoundary
	    && n2->getBoundaryFlag() != kClosedBoundary
	    && !( n1->getBoundaryFlag()==kOpenBoundary
		  && n2->getBoundaryFlag()==kOpenBoundary ) ) ?
     kFlowAllowed : kFlowNotAllowed;
}

inline void tEdge::setFlowAllowed( const tNode* n1, const tNode* n2 )
{
   flowAllowed = tEdge::isFlowAllowed(n1, n2);
}

inline void tEdge::setCCWEdg( tEdge * edg )
{
   assert( edg != 0 );
   ccwedg = edg;
}

inline void tEdge::setCWEdg( tEdge * edg )
{
  cwedg = edg;
}

inline void tEdge::setRVtx( tArray< double > const & arr )
{
   assert( arr.getSize() == 2 );
   if (0)//DEBUG
     cout << "setRVtx for edge " << id
	  << " to x, y, " << arr[0] << ", " << arr[1] << endl;
   rvtx = arr;
}

inline void tEdge::setVEdgLen( double val )
{
   assert( val>=0.0 );
   vedglen = val;
}

inline void tEdge::setTri( tTriangle* tptr )
{
  tri = tptr;
}

inline bool tEdge::isFlippable() const
{
  // An edge is not flippable if
  // - both nodes are mobile
  // - the edge is a flow edge for the origin or the complement
  //   edge is a flow edge for the destination
  return !(
	   org->isMobile() && dest->isMobile() &&
	   (org->flowThrough(this) ||
	    dest->flowThrough(this->getComplementEdge())
	    )
	   );
}

/**************************************************************************\
**
**  tEdge::CalcSlope
**
**  Computes the slope of the edge as ( Zorg - Zdest ) / length.
**
**  Returns: the slope
**  Modifies: slope (data mbr)
**  Assumes: length >0; org and dest valid.
**
\**************************************************************************/
inline double tEdge::CalcSlope()
{
   assert( org!=0 );  // Failure = edge has no origin and/or destination node
   assert( dest!=0 );
   assert( len>0.0 );

   slope = ( org->getZ() - dest->getZ() ) / len;
   return slope;
}


/**************************************************************************\
**
**  tEdge::CalcVEdgLen
**
**  Calculates the length of the Voronoi cell edge associated with the
**  current triangle edge. The Voronoi cell edge length is equal to the
**  distance between the Voronoi vertex of the right-hand triangle and
**  the Voronoi vertex of the left-hand triangle. The vertex for the
**  right-hand triangle is stored with rvtx[] (and is assumed to be up to
**  date), and the vertex for the left-hand triangle is stored in the
**  edge's counter-clockwise (left-hand) neighbor (also assumed valid and
**  up to date).
**
**  Data mbrs modified:  vedglen
**  Returns:  the Voronoi edge length
**  Assumes:  ccwedg valid, rvtx[] up to date
**
\**************************************************************************/
inline double tEdge::CalcVEdgLen()
{
	assert( ccwedg!=0 );

	double dx, dy;

	dx = rvtx[0] - ccwedg->rvtx[0];
	dy = rvtx[1] - ccwedg->rvtx[1];
	vedglen = sqrt( dx*dx + dy*dy );
	return( vedglen );
}

/**************************************************************************\
**  Functions for class tTriangle.
\**************************************************************************/

/***********************************************************************\
**
**  Constructors & destructors:
**
**  Default:  initializes node, edge, and triangle ptrs to zero.
**  Copy:  copies all values
**  ID & Vertices: creates a triangle w/ pointers to three vertices,
**                 and sets up edge pointers as well. Does not set
**                 triangle pointers however (these are zero'd).
**  Destructor:  no longer used
**
**  Modifications:
**   - ID & vertices constructor added 1/2000, GT
**
\***********************************************************************/
inline
void tTriangle::SetIndex()
{
  index_[0] = 0;
  index_[1] = 1;
  index_[2] = 2;
}

//default
inline tTriangle::tTriangle() :
  id(-1)
{
   for( int i=0; i<3; i++ )
   {
      p[i] = 0;
      e[i] = 0;
      t[i] = 0;
   }
   SetIndex();
   if (0)//DEBUG
     cout << "tTriangle()" << endl;
}

//copy constructor
inline tTriangle::tTriangle( const tTriangle &init ) :
  id(init.id)
{
   for( int i=0; i<3; i++ )
     {
       p[i] = init.p[i];
       setEPtr( i, init.e[i] ); // sets edge's tri pointer as well!
       t[i] = init.t[i];
       index_[i] = init.index_[i];
     }
   if (0)//DEBUG
     cout << "tTriangle( orig )" << endl;
}

// construct with id and 3 vertices
inline tTriangle::tTriangle( int id_, tNode* n0, tNode* n1, tNode* n2 ) :
  id(id_)
{
  assert( n0 != 0 && n1 != 0 && n2 != 0 );
  p[0] = n0;
  p[1] = n1;
  p[2] = n2;
  setEPtr( 0, n0->EdgToNod( n2 ) );
  setEPtr( 1, n1->EdgToNod( n0 ) );
  setEPtr( 2, n2->EdgToNod( n1 ) );
  t[0] = t[1] = t[2] = 0;
  SetIndex();
}

inline tTriangle::tTriangle( int id_, tNode* n0, tNode* n1, tNode* n2,
                             tEdge* e0, tEdge* e1, tEdge* e2 ) :
  id(id_)
{
  assert( n0 != 0 && n1 != 0 && n2 != 0 && e0 != 0 && e1 != 0 && e2 != 0 );
  p[0] = n0;
  p[1] = n1;
  p[2] = n2;
  setEPtr( 0, e0 );
  setEPtr( 1, e1 );
  setEPtr( 2, e2 );
  t[0] = t[1] = t[2] = 0;
  SetIndex();
}
/***********************************************************************\
**
**  Overloaded operators:
**
**    assignment: copies all values
**    left shift: sends the following data to the output stream:
**                triangle ID and the IDs of its 3 nodes, clockwise
**                edges, ad neighboring triangles (or -1 if no
**                neighboring triangle exists across a given face)
**    right shift: reads triangle ID and 3 other unspecified IDs from
**                 the input stream (the latter are not currently used
**                 for anything)
**
\***********************************************************************/

//overloaded assignment operator
inline const tTriangle &tTriangle::operator=( const tTriangle &init )
{
   if( &init != this )
   {
      id = init.id;
      for( int i=0; i<3; i++ )
      {
         p[i] = init.p[i];
         setEPtr( i, init.e[i] ); // sets edge's tri pointer as well!
         t[i] = init.t[i];
	 index_[i] = init.index_[i];
      }
   }
   return *this;
}

//left shift
inline ostream &operator<<( ostream &output, const tTriangle &tri )
{
   int i;
   output << tri.id << ":";
   for( i=0; i<3; i++ )
       output << " " << tri.p[i]->getID();
   output << ";";
   for( i=0; i<3; i++ )
       output << " " << tri.e[i]->getID();
   output << ";";
   for( i=0; i<3; i++ )
   {
      if( tri.t[i] != 0 ) output << " " << tri.t[i]->getID();
      else  output << " -1";
   }
   output << endl;
   return output;
}

inline istream &operator>>( istream &input, tTriangle &tri )
{
   int id1, id2, id3;
   cout << "triangle id, origin id, dest id:";
   input >> tri.id >> id1 >> id2 >> id3; //temporarily assign id vals to ptrs
   return input;
}

/***********************************************************************\
**
**  tTriangle "get" functions:
**
**  getID - returns ID #
**  pPtr - returns ptr to one of the 3 vertex nodes, as specified by
**         _index_ (index is 0, 1, or 2)
**  ePtr - returns ptr to one of the 3 clockwise edges, as specified by
**         _index_ (index is 0, 1, or 2)
**  tPtr - returns ptr to one of the 3 adjacent triangles, as specified
**         by _index_ (index is 0, 1, or 2)
**
\***********************************************************************/

inline int tTriangle::getID() const {return id;}

inline tNode *tTriangle::pPtr( int index ) const
{
   assert( index >= 0 && index < 3 );
   return p[index];
}

inline tEdge *tTriangle::ePtr( int index ) const
{
   assert( index >= 0 && index < 3 );
   return e[index];
}

inline tTriangle *tTriangle::tPtr( int index ) const
{
   assert( index >= 0 && index < 3 );
   return t[index];
}


/***********************************************************************\
**
**  tTriangle "set" functions:
**
**  setID - sets ID #
**  setPPtr - sets pointer to one of the 3 vertex nodes, as specified
**            by _index_ (index is 0, 1, or 2)
**  setEPtr - sets pointer to one of the 3 clockwise edges, as specified
**            by _index_ (index is 0, 1, or 2)
**  setTPtr - sets pointer to one of the 3 adjacent triangles, as
**            specified by _index_ (index is 0, 1, or 2)
**
\***********************************************************************/
inline void tTriangle::setID( int val ) {id = ( val >= 0 ) ? val : 0;}

inline void tTriangle::setPPtr( int index, tNode * ndptr )
{
   assert( index >= 0 && index < 3 );
   p[index] = ndptr;
}

inline void tTriangle::setEPtr( int index, tEdge * egptr )
{
  assert( index >= 0 && index < 3 );
  if( egptr == 0 && e[index] != 0 &&
      e[index]->TriWithEdgePtr() == this) {
    e[index]->setTri(0);
  }
  e[index] = egptr;
  if( egptr != 0 )
    egptr->setTri( this );
}

inline void tTriangle::setTPtr( int index, tTriangle * trptr )
{
   assert( index >= 0 && index < 3 );
   t[index] = trptr;
}


/**************************************************************************\
**
**  tTriangle::nVOp
**
**  Returns the side number (0, 1, or 2) of the neighboring triangle ct.
**  Assumes that ct _is_ one of the neighboring triangles.
**
\**************************************************************************/
inline int tTriangle::nVOp( const tTriangle *ct ) const
{
  for( int i=0; i<3; ++i )
    if( t[i] == ct ) return i;
  assert( 0 ); /*NOTREACHED*/
  abort();
}


/**************************************************************************\
**
**  tTriangle::nVtx
**
**  Returns the vertex number (0, 1, or 2) associated with node cn.
**  (In other words, it says whether cn is vertex 0, 1, or 2 in the
**  triangle).
**  Assumes that cn _is_ one of the triangle's vertices.
**
\**************************************************************************/
inline int tTriangle::nVtx( const tNode *cn ) const
{
  for( int i=0; i<3; ++i )
    if( p[i] == cn ) return i;
  assert( 0 ); /*NOTREACHED*/
  abort();
}

/*****************************************************************************\
**
**  tTriangle::isIndexIDOrdered
**
**  Tell whether index_ has been ID ordered
**
\*****************************************************************************/
inline bool tTriangle::isIndexIDOrdered() const {
  return
    BOOL(
	 ( index_[1] == (index_[0]+1)%3 ) &&
	 ( index_[2] == (index_[1]+1)%3 ) &&
	 ( pPtr(index_[0])->getID() < pPtr(index_[1])->getID() ) &&
	 ( pPtr(index_[0])->getID() < pPtr(index_[2])->getID() )
	 );
}

/**************************************************************************\
**
**  tSpkIter
**
\**************************************************************************/
inline tSpkIter::tSpkIter() :
  curnode(0),
  curedg(0),
  counter(0)
{}

inline tSpkIter::tSpkIter( tNode* nPtr ) :
  curnode(0),
  curedg(0),
  counter(0)
{
   assert( nPtr != 0 );
   curnode = nPtr;
   curedg = curnode->getEdg();
}

inline tSpkIter::~tSpkIter()
{
   curnode = 0;
   curedg = 0;
}

inline tNode* tSpkIter::CurNode() {return curnode;}

inline tEdge* tSpkIter::CurSpoke() {return curedg;}

inline int tSpkIter::getNumSpokes()
{
   tEdge* fe = curnode->getEdg();
   tEdge* ce = fe;
   int ctr = 1;
   while( ( ce = ce->getCCWEdg() ) != fe )
       ++ctr;
   return ctr;
}

inline void tSpkIter::Reset( tNode* nPtr )
{
   assert( nPtr != 0 );
   curnode = nPtr;
   curedg = curnode->getEdg();
   counter = 0;
}

inline int tSpkIter::First()
{
   curedg = curnode->getEdg();
   counter = 0;
   if( curedg != 0 ) return 1;
   return 0;
}

inline tEdge* tSpkIter::FirstP()
{
   if( First() ) return curedg;
   return 0;
}


inline int tSpkIter::Last()
{
   curedg = curnode->getEdg()->getCWEdg();
   counter = -1;
   if( curedg != 0 ) return 1;
   return 0;
}

inline tEdge* tSpkIter::LastP()
{
   if( Last() ) return curedg;
   return 0;
}


inline int tSpkIter::Next()
{
   curedg = curedg->getCCWEdg();
   ++counter;
   if( curedg != 0 ) return 1;
   return 0;
}

inline tEdge* tSpkIter::NextP()
{
   if( Next() ) return curedg;
   return 0;
}

inline int tSpkIter::Prev()
{
   curedg = curedg->getCWEdg();
   --counter;
   if( curedg != 0 ) return 1;
   return 0;
}

inline tEdge* tSpkIter::PrevP()
{
   if( Prev() ) return curedg;
   return 0;
}

inline int tSpkIter::Get( int num )
{
   if( !First() ) return 0;
   while( curedg->getID() != num && !AtEnd() )
       curedg = curedg->getCCWEdg();
   if( !AtEnd() ) return 1;
   return 0;
}

inline tEdge* tSpkIter::GetP( int num )
{
   if( Get( num ) ) return curedg;
   return 0;
}

inline int tSpkIter::Get( const tEdge* ePtr )
{
   if( ePtr == 0 ) return 0;
   if( !First() ) return 0;
   while( curedg != ePtr && !AtEnd() )
       curedg = curedg->getCCWEdg();
   if( !AtEnd() ) return 1;
   return 0;
}

inline tEdge* tSpkIter::GetP( const tEdge* ePtr )
{
   if( Get( ePtr ) ) return curedg;
   return 0;
}

inline int tSpkIter::Where() const
{
   if( curedg == 0 ) return -1;
   return curedg->getID();
}

inline tEdge* tSpkIter::ReportNextP()
{
   return curedg->getCCWEdg();
}

inline tEdge* tSpkIter::ReportPrevP()
{
   return curedg->getCWEdg();
}

inline int tSpkIter::AtEnd()
{
   if( isEmpty() ) return 1;
   return ( curedg == curnode->getEdg() && counter != 0 );
}

inline int tSpkIter::isEmpty()
{
   return( curnode->getEdg() == 0 );
}

inline int tSpkIter::insertAtPrev( tEdge* ePtr )
{
   if( ePtr == 0 ) return 0;
   if( isEmpty() ) return insertAtFront( ePtr );
   assert( ePtr->getOriginPtr() == curnode );
   ePtr->setCCWEdg( curedg );
   ePtr->setCWEdg( curedg->getCWEdg() );
   curedg->getCWEdg()->setCCWEdg( ePtr );
   curedg->setCWEdg( ePtr );
   return 1;
}

inline int tSpkIter::insertAtNext( tEdge* ePtr )
{
   if( ePtr == 0 ) return 0;
   if( isEmpty() ) return insertAtFront( ePtr );
   assert( ePtr->getOriginPtr() == curnode );
   ePtr->setCWEdg( curedg );
   ePtr->setCCWEdg( curedg->getCCWEdg() );
   curedg->getCCWEdg()->setCWEdg( ePtr );
   curedg->setCCWEdg( ePtr );
   return 1;
}

inline int tSpkIter::insertAtFront( tEdge* ePtr )
{
   if( ePtr == 0 ) return 0;
   assert( ePtr->getOriginPtr() == curnode );
   if( isEmpty() )
   {
      curnode->setEdg( ePtr );
      ePtr->setCCWEdg( ePtr );
      ePtr->setCWEdg( ePtr );
      return 1;
   }
   tEdge* ofe = curnode->getEdg();
   ePtr->setCCWEdg( ofe );
   ePtr->setCWEdg( ofe->getCWEdg() );
   ofe->getCWEdg()->setCCWEdg( ePtr );
   ofe->setCWEdg( ePtr );
   curnode->setEdg( ePtr );
   return 1;
}

inline int tSpkIter::insertAtBack( tEdge* ePtr )
{
   if( ePtr == 0 ) return 0;
   if( isEmpty() ) return insertAtFront( ePtr );
   assert( ePtr->getOriginPtr() == curnode );
   tEdge* ofe = curnode->getEdg();
   ePtr->setCCWEdg( ofe );
   ePtr->setCWEdg( ofe->getCWEdg() );
   ofe->getCWEdg()->setCCWEdg( ePtr );
   ofe->setCWEdg( ePtr );
   return 1;
}

inline tEdge* tSpkIter::removePrev()
{
   if( isEmpty() ) return 0;
   if( curedg->getCWEdg() == curedg )
   {
      curnode->setEdg( 0 );
      return curedg;
   }
   tEdge* re = curedg->getCWEdg();
   re->getCWEdg()->setCCWEdg( re->getCCWEdg() );
   re->getCCWEdg()->setCWEdg( re->getCWEdg() );
   return re;
}

inline tEdge* tSpkIter::removeNext()
{
   if( isEmpty() ) return 0;
   if( curedg->getCCWEdg() == curedg )
   {
      curnode->setEdg( 0 );
      return curedg;
   }
   tEdge* re = curedg->getCCWEdg();
   re->getCWEdg()->setCCWEdg( re->getCCWEdg() );
   re->getCCWEdg()->setCWEdg( re->getCWEdg() );
   return re;
}

inline tEdge* tSpkIter::removeFromFront()
{
   if( isEmpty() ) return 0;
   tEdge* re = curnode->getEdg();
   if( re->getCCWEdg() == re )
   {
      curnode->setEdg( 0 );
      return re;
   }
   curnode->setEdg( re->getCCWEdg() );
   re->getCWEdg()->setCCWEdg( re->getCCWEdg() );
   re->getCCWEdg()->setCWEdg( re->getCWEdg() );
   if( curedg == re ) curedg = curnode->getEdg();
   return re;
}

inline tEdge* tSpkIter::removeFromBack()
{
   if( isEmpty() ) return 0;
   tEdge* re = curnode->getEdg()->getCWEdg();
   if( re->getCCWEdg() == re )
   {
      curnode->setEdg( 0 );
      return re;
   }
   re->getCWEdg()->setCCWEdg( re->getCCWEdg() );
   re->getCCWEdg()->setCWEdg( re->getCWEdg() );
   if( curedg == re ) curedg = curnode->getEdg()->getCWEdg();
   return re;
}

#endif
