/**************************************************************************\
**
**  meshElements.h: Header file for mesh elements tNode, tEdge,
**                  and tTriangle. Each of these mesh elements is
**                  implemented as an object, as described below.
**                  (formerly called gridElements)
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
**  $Id: meshElements.h,v 1.27 2002-04-23 10:02:10 arnaud Exp $
**  (file consolidated from earlier separate tNode, tEdge, & tTriangle
**  files, 1/20/98 gt)
\**************************************************************************/

#ifndef MESHELEMENTS_H
#define MESHELEMENTS_H

#include <iostream.h>
#include <math.h>       // for sqrt() used in inlined fn below
#include "../Definitions.h"
#include "../tPtrList/tPtrList.h"
#include "../tArray/tArray.h"
#include "../tMeshList/tMeshList.h"
#include "../Geometry/geometry.h"   // for Point2D definitions & fns

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
\**************************************************************************/
class tNode
{
  friend ostream &operator<<( ostream &, tNode & );
  friend istream &operator>>( istream &, tNode & );

public:

  tNode();                                   // default constructor
  tNode( const tNode & );                    // copy constructor
  virtual ~tNode() {}

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
  void getVoronoiVertexList( tList<Point2D> * );  // Returns list of V vertices
  void getVoronoiVertexXYZList( tList<Point3D> * ); // As above plus interp z
  
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
  //~tEdge();               // destructor

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
  tTriangle( int, tNode *, tNode *, tNode * );
  //~tTriangle();                   // destructor

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
  x(0.), y(0.), z(0.), varea(0.), varea_rcp(0.),
  boundary(0), edg(0)
{}

//copy constructor
inline tNode::tNode( const tNode &original )
{
   if( &original != 0 )
   {
      id = original.id;
      x = original.x;
      y = original.y;
      z = original.z;
      boundary = original.boundary;
      varea = original.varea;
      varea_rcp = original.varea_rcp;
      edg = original.edg;
      if( &(original.spokeList) != 0 )
      {
         for( int i=0; i<original.spokeList.getSize(); i++ )
         {
            insertBackSpokeList( original.spokeList.getIthPtrNC( i ) );
         }
      }
        //else spokeList = 0;
   }
     //cout << "tNode( original )" << endl;
}

/*X tNode::~tNode()                                                      //tNode
{
     //if( spokeList != 0 ) delete spokeList;
     //cout << "    ~tNode()" << endl;
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
      spokeList = right.spokeList;
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
   //tPtrListIter< tEdge > spokeIter( node.getSpokeListNC() );
   
   output << node.id << ": " << node.x << " " << node.y << " "
          << node.z << ";";
   output <<
       node.spokeList.getFirstNC()->getPtrNC()->getDestinationPtrNC()->getID()
   //output << spokeIter.DatPtr()->getDestinationPtrNC()->getID()
          << " ";
   output << endl;
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
**  getSpokeList - returns const reference to spoke list
**  getSpokeListNC - returns non-const reference to spoke list
**  getFirstSpokeNode - returns ptr to 1st spokelist item 
**                      (return type tPtrListNode *)
**  getFirstSpokeNodeNC - non-const version of the above
**  getNextSpokeNode - returns ptr to next spokelist item that follows
**                       prevedg
**  getNextSpokeNodeNC - non-const version of the above
**
\***********************************************************************/

inline tArray< double >
tNode::get3DCoords() const
{
   tArray< double > xyz(3);
   xyz[0] = x;
   xyz[1] = y;
   xyz[2] = z;
   return xyz;
}

inline tArray< double >
tNode::get2DCoords() const
{
   tArray< double > xy(2);
   xy[0] = x;
   xy[1] = y;
   return xy;
}

inline int tNode::getID() const {return id;}                    //tNode
inline double tNode::getX() const {return x;}
inline double tNode::getY() const {return y;}
inline double tNode::getZ() const {return z;}
inline double tNode::getVArea() const {return varea;}             //tNode
inline double tNode::getVArea_Rcp() const {return varea_rcp;}     //tNode
inline int tNode::getBoundaryFlag() const {return boundary;}      //tNode
inline tEdge * tNode::getEdg() {return edg;}

inline const tPtrList< tEdge > &                                     //tNode
tNode::getSpokeList() const {assert( &spokeList != 0 ); return spokeList;}

inline tPtrList< tEdge > &                                           //tNode
tNode::getSpokeListNC() {assert( &spokeList != 0 ); return spokeList;}

inline const tPtrListNode< tEdge > *                                 //tNode
tNode::getFirstSpokeNode() const {return spokeList.getFirst();}

inline tPtrListNode< tEdge > *                                       //tNode
tNode::getFirstSpokeNodeNC() {return spokeList.getFirstNC();}

inline const tPtrListNode< tEdge > *tNode::                          //tNode
getNextSpokeNode( const tPtrListNode< tEdge > *prevedg ) const
{return prevedg->getNext();}

inline tPtrListNode< tEdge > *tNode::                                //tNode
getNextSpokeNodeNC( tPtrListNode< tEdge > *prevedg ) const
{return prevedg->getNextNC();}


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

inline void tNode::setBoundaryFlag( int val )
{
  assert( val>=0 && val<=2 );
  boundary = val;
  /*boundary = (val >=0 && val <= 2) ? val : 0;*/
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
   assert( theEdg > 0 );
   edg = theEdg;
   //cout << "Assigning edge " << theEdg->getID() << " to node " << getID() << endl;
}


/***********************************************************************\
**
**  tNode::ChangeZ:  Adds delz to current z value
**
\***********************************************************************/
inline void tNode::ChangeZ( double delz ) { z += delz; }      //tNode

/***********************************************************************\
**
**  tNode::makeWheel:  makes the spoke list circular
**
\***********************************************************************/
inline void tNode::makeWheel() {spokeList.makeCircular();}

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
inline tEdge::tEdge()
        : rvtx(2)
{
   id = 0;
   len = 0;
   slope = 0;
   vedglen = 0;
   org = dest = 0;
   ccwedg = 0;
   flowAllowed = 0;
     //cout << "tEdge()" << endl;
}

//copy constructor
inline tEdge::tEdge( const tEdge &original )
{
   if( &original != 0 )
   {
      id = original.id;
      len = original.len;
      slope = original.slope;
      rvtx = original.rvtx;
      vedglen = original.vedglen;
      org = original.org;
      dest = original.dest;
      ccwedg = original.ccwedg;
      flowAllowed = original.flowAllowed;
   }
     //cout << "tEdge( orig )" << endl;
}

//tEdge::~tEdge() {/*cout << "    ~tEdge()" << endl;*/}      //tEdge


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
      flowAllowed = original.flowAllowed;
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

inline int tEdge::getID() const {return id;}                                //tEdge

//return 0 if flow allowed to match kNonBoundary:
inline int tEdge::getBoundaryFlag() const             
{return !( flowAllowed == kFlowAllowed );} 

inline double tEdge::getLength() const {return len;}  

inline double tEdge::getSlope() const {return slope;} 

inline const tNode *tEdge::getOriginPtr() const {return org;} 

inline const tNode *tEdge::getDestinationPtr() const {return dest;}

inline tNode *tEdge::getOriginPtrNC() {return org;}                

inline tNode *tEdge::getDestinationPtrNC() {return dest;}          

inline double tEdge::getOrgZ() 
{
   //Xconst tNode * org = getOriginPtr(); 5/99
   assert( org!=0 );
   return( org->getZ() );
}

inline double tEdge::getDestZ()
{
   //Xconst tNode * dest = getDestinationPtr(); 5/99
   assert( dest!=0 );
   return( dest->getZ() );
}

inline tEdge * tEdge::getCCWEdg() 
{
   return ccwedg;
}

inline int tEdge::FlowAllowed() 
{
   return flowAllowed;
}

inline tArray< double >
tEdge::getRVtx() const
{
   //tArray< double >  xy( rvtx );
   //return xy;
   //cout << "getRVtx: ";
   return rvtx;
}

inline double tEdge::getVEdgLen() const {return vedglen;}


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
   /*id = ( val >=0 ) ? val : 0;*/
}           //tEdge

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

inline void tEdge::setFlowAllowed( int val )
{
   assert( val==0 || val==1 );
   flowAllowed = val;
   /*flowAllowed = ( val == 0 || val == 1 ) ? val : 0;*/}

inline void tEdge::setCCWEdg( tEdge * edg )
{
   assert( edg > 0 );
     //assert( ccwedg > 0 );
   ccwedg = edg;
}

inline void tEdge::setRVtx( tArray< double > arr )
{
   assert( &arr != 0 );
   assert( arr.getSize() == 2 );
     //cout << "setRVtx for edge " << id
     //   << " to x, y, " << arr[0] << ", " << arr[1] << endl;
   rvtx = arr;
}

inline void tEdge::setVEdgLen( double val )
{
   assert( val>=0.0 );
   vedglen = val;
   /*vedglen = ( val > 0 ) ? val : 0;*/
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
   //Xconst tNode * org = getOriginPtr(); 5/99
   //Xconst tNode * dest = getDestinationPtr(); 5/99
   
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
**
**  tEdge::WelcomeCCWNeighbor
**
**  Welcomes a new spoke to the neighborhood! neighbor is a new edge to
**  be inserted counter-clockwise from this one. We point neighbor at
**  the edge we're currently pointing to, and then point ourself to
**  neighbor, thus maintaining the edge connectivity.
**
**  Data mbrs modified:  ccwedg, neighbor->ccwedg
**  Created: 2/4/99 GT
**
\**************************************************************************/
inline void tEdge::WelcomeCCWNeighbor( tEdge * neighbor )
{
   assert( neighbor!=0 );
   assert( neighbor->org == org );
   neighbor->ccwedg = ccwedg;
   ccwedg = neighbor;
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

//default
inline tTriangle::tTriangle()
{
   assert( p != 0 && e != 0 && t != 0 );
   for( int i=0; i<3; i++ )
   {
      p[i] = 0;
      e[i] = 0;
      t[i] = 0;
   }
     //cout << "tTriangle()" << endl;
}

//copy constructor
inline tTriangle::tTriangle( const tTriangle &init )
{
   assert( p != 0 && e != 0 && t != 0 );
   if( &init != 0 )
   {
      id = init.id;
      for( int i=0; i<3; i++ )
      {
         p[i] = init.p[i];
         e[i] = init.e[i];
         t[i] = init.t[i];
      }
   }
     //cout << "tTriangle( orig )" << endl;
}

// construct with id and 3 vertices
inline tTriangle::tTriangle( int num, tNode* n0, tNode* n1, tNode* n2 )
{
   id = num;
   assert( n0 > 0 && n1 > 0 && n2 > 0 );
   p[0] = n0;
   p[1] = n1;
   p[2] = n2;
   setEPtr( 0, n0->EdgToNod( n2 ) );
   setEPtr( 1, n1->EdgToNod( n0 ) );
   setEPtr( 2, n2->EdgToNod( n1 ) );
   t[0] = t[1] = t[2] = 0;
}

//destructor
/*tTriangle::~tTriangle()                                          //tTriangle
{
     //cout << "    ~tTriangle()" << endl;
}*/


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
         e[i] = init.e[i];
         t[i] = init.t[i];
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
     //tri.setPPtr( tMesh::h.getList().
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

inline tNode *tTriangle::pPtr( int index )
{
   assert( index >= 0 && index <= 3 );
   return p[index];
}

inline tEdge *tTriangle::ePtr( int index )
{
   assert( index >= 0 && index <= 3 );
   return e[index];
}

inline tTriangle *tTriangle::tPtr( int index )
{
   assert( index >= 0 && index <= 3 );
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
   assert( index >= 0 && index <= 3 );
   p[index] = ndptr;
}

inline void tTriangle::setEPtr( int index, tEdge * egptr )
{
   assert( index >= 0 && index <= 3 );
   e[index] = egptr;
}

inline void tTriangle::setTPtr( int index, tTriangle * trptr )
{
   assert( index >= 0 && index <= 3 );
   t[index] = trptr;
}


/**************************************************************************\
**
**  tTriangle::nVOp
**
**  Returns the side number (0, 1, or 2) of the neighboring triangle ct.
**  Assumes that ct _is_ one of the neighboring triangles.
**
** NOTE: for runtime error checking, may want to take out the assert
** and instead generate a runtime error when i>2.
\**************************************************************************/
inline int tTriangle::nVOp( tTriangle *ct )
{
   int i;

   for( i=0; i<4; i++ )
   {
      assert( i<3 );
      if( t[i] == ct ) return i;
   }
   return i;
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
** NOTE: for runtime error checking, may want to take out the assert
** and instead generate a runtime error when i>2.
\**************************************************************************/
inline int tTriangle::nVtx( tNode *cn )
{
   int i;
   for( i=0; i<4; i++ )
   {
      assert( i<3 );
      if( p[i] == cn ) return i;
   }
   return i;
}

#endif
