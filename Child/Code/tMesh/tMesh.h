//-*-c++-*-

/***************************************************************************/
/**
**  @file tMesh.h
**  @brief Header file for class tMesh (formerly tGrid)
**
**  tMesh is the master class that handles the implementation of the
**  triangulated mesh. The class includes lists of the mesh elements
**  (nodes, triangles, and edges; see meshElements.h/.cpp), and provides
**  functionality to:
**    - read in or create meshes, either from scratch, from a list of
**      points, from a pre-existing set of triangulation files (e.g., a
**      previous run), or an Arc/Info grid DEM
**    - move, add, and/or delete nodes
**    - update Delaunay and Voronoi geometry
**
**  Summary of recent changes:
**    - added data mbr mSearchOriginTriPtr to implement triangle search
**      starting from a given location, GT, 1/2000
**    - added default argument "interpFlag" to MoveNodes() in order
**      to have nodes moved w/o interpolation (eg, for tectonic movement)
**      (GT, 4/00)
**
**  $Id: tMesh.h,v 1.56 2003-05-23 11:57:56 childcvs Exp $
*/
/***************************************************************************/

#ifndef TMESH_H
#define TMESH_H

#include "../tAssert.h"
#include <math.h>
#include <stdlib.h>
#include "../Classes.h"
#include "../Definitions.h"
#include "../tArray/tArray.h"
#include "../tMatrix/tMatrix.h"
#include "../Mathutil/mathutil.h"
#include "../errors/errors.h"
#include "../tPtrList/tPtrList.h"
#include "../tList/tList.h"
#include "../tMeshList/tMeshList.h"
#include "../MeshElements/meshElements.h"
#include "../tLNode/tLNode.h"
#include "../tListInputData/tListInputData.h"
//#include "../tListOutputData/tListOutputData.h"
#include "../globalFns.h"
#include "../Predicates/predicates.h"


/***************************/
/**** Defined constants ****/
/***************************/

/* codes passed to DeleteNode (why don't default args work?) */
#define kRepairMesh true
#define kNoRepair false


/****************************/
/**** Class Declarations ****/
/****************************/

class ParamMMFS_t;

/** @class tMesh
 */
template< class tSubNode >
class tMesh
{
   friend class tListOutputData< tSubNode >;
   tMesh(const tMesh&);
   tMesh& operator=(const tMesh&);

   void MakePointBoundary( const ParamMMFS_t &, tInputFile &, tPtrList< tSubNode > & );
   void MakePointInterior( const ParamMMFS_t &, tInputFile &,
			   bool makeMesh);
   void BuildDelaunayMeshTipper();
   void MeshDensification( tInputFile & );
   void MakeDelaunay( tPtrList< tTriangle > &);
public:
   tMesh();
   tMesh( tInputFile & );
   tMesh( tMesh const * );
   ~tMesh();
   void BatchAddNodes(); // quickly adds many nodes when starting w/ dense mesh
   void MakeMeshFromScratch( tInputFile & );   // creates a new mesh
   void MakeMeshFromScratchTipper( tInputFile & );   // creates a new mesh
   void MakeMeshFromInputData( tInputFile & ); // reads in an existing mesh
   void MakeMeshFromPoints( tInputFile & );    // creates mesh from list of pts
   void MakeMeshFromPointsTipper( tInputFile & ); // creates mesh from list of pts
   void MakeRandomPointsFromArcGrid( tInputFile & ); // mesh from arc (rand)
   void MakeHexMeshFromArcGrid( tInputFile &infile );// mesh from arc (hex)
   void MakeLayersFromInputData( tInputFile & );
   void Print() const;
   void setVoronoiVertices();
   void CalcVoronoiEdgeLengths();
   void CalcVAreas();
   tTriangle *LocateTriangle( double, double );
   tTriangle *LocateNewTriangle( double, double );
   /*returns ptr to triangle which points to edge, or zero if none:*/
   tTriangle *TriWithEdgePtr( tEdge * ) const;
   /*only routine needed to delete node; calls ExNode, RepairMesh:*/
   int DeleteNode( tListNode< tSubNode > *, bool repairFlag=true );
   int DeleteNode( tSubNode const *, bool repairFlag=true );
   /*deletes spokes, *calls ExEdge, makes nbr ptr list:*/
   int ExtricateNode( tSubNode *, tPtrList< tSubNode > & );
   int DeleteEdge( tEdge * );
   /*calls ExTriangle; deals w/edge and compliment, deletes nbr tri(s):*/
   int ExtricateEdge( tEdge * );
   int ClearEdge( tEdge* ) const; // calls ClearTriangle
   int DeleteTriangle( tTriangle const * );
   /*"un-points" nbr triangles:*/
   int ExtricateTriangle( tTriangle const * );
   int ClearTriangle( tTriangle const* ) const;
   /*complicated; fills in (any) hole defined by circular node ptr list:*/
   int RepairMesh( tPtrList< tSubNode > & );
   int AddEdgeAndMakeTriangle( tPtrList< tSubNode > &,
                               tPtrListIter< tSubNode > & );
   int AddEdgeAndMakeTriangle( tSubNode*, tSubNode*, tSubNode* );
   int MakeTriangle( tPtrList< tSubNode > const&,
                     tPtrListIter< tSubNode > & );
   int MakeTriangle( tSubNode*, tSubNode*, tSubNode* );
   int AddEdge( tSubNode *, tSubNode *, tSubNode const * );
   void CheckTrianglesAt( tSubNode * );
   tSubNode *AddToList( tSubNode const& );
   //add a node with referenced value/properties, update mesh connectivity
   tSubNode *AddNode( tSubNode &, bool updatemesh = false, double time = 0.0 );
   tSubNode* InsertNode( tSubNode*, double );
   //add a generic node at the referenced coordinates
   tSubNode *AddNodeAt( tArray< double > &, double time = 0.0 );
   tSubNode* AttachNode( tSubNode*, tTriangle* );
   tMeshList<tEdge> * getEdgeList();
   tMeshList<tSubNode> * getNodeList();
   tList< tTriangle > * getTriList();
   tEdge *getEdgeComplement( tEdge * ) const;
   /* Tests consistency of a user-defined mesh */
   void CheckMeshConsistency( bool boundaryCheckFlag=true );
   /* Updates mesh by comp'ing edg lengths & slopes & node Voronoi areas */
   void UpdateMesh();
   /* computes edge slopes as (Zorg-Zdest)/Length */
   //void CalcSlopes(); /* WHY is this commented out? */
   /*routines used to move points; MoveNodes is "master" function*/
   bool CheckForFlip( tTriangle *, int, bool, bool useFuturePosn=true );
   void FlipEdge( tTriangle *, tTriangle *, int, int );
   tEdge * IntersectsAnyEdge( tEdge * );
   //CheckTriEdgeIntersect and CheckLocallyDelaunay are the essential functions
   //for point moving and are called from MoveNodes; in general, they should be
   //after the three above functions and in the order listed below.
   void CheckTriEdgeIntersect();
   void CheckLocallyDelaunay();
   //once 'newx' and 'newy' are set, this is the only function one needs to call
   //to execute the move; maybe should have a separate function for doing things
   //peculiar to channels, but now this is the only one.
   void MoveNodes( double time = 0.0, bool interpFlag=true );
   /*end moving routines*/

   void AddNodesAround( tSubNode *, double time=0.0 );  // Local mesh densify

   void ResetNodeID(); // reset node IDs in list order
   void ResetEdgeID(); // reset edge IDs in list order
   void ResetTriangleID(); // reset triangle IDs in list order
   void SetmiNextNodeID(int);
   void SetmiNextEdgID(int);
   void SetmiNextTriID(int);

#ifndef NDEBUG
   /*'dump' routines for debugging*/
   void DumpEdges();
   void DumpSpokes( tSubNode * ) const;
   void DumpTriangles();
   void DumpNodes();
#endif

protected:
   int nnodes, nedges, ntri;       // # of nodes, edges, and tri's (obsolete?)
   tMeshList< tSubNode > nodeList; // list of nodes
   tMeshList< tEdge > edgeList;    // list of directed edges
   tList< tTriangle > triList;  // list of ptrs to triangles
   int miNextNodeID;                   // next ID for added node
   int miNextEdgID;                    // next ID for added edge
   int miNextTriID;                    // next ID for added triangle
   long seed;                      // random seed
   bool layerflag;                 // flag indicating whether nodes have layers
   tTriangle* mSearchOriginTriPtr; // ptr to tri. from which to start searches

};

/** @class tIdArray
    @brief Lookup table per Id for a tList
*/
template< class T >
class tIdArray
{
  tArray< T* > e_;
public:
  tIdArray(tList< T >& List);
  T* operator[]( int subscript ) const {
    // return a value and not a reference, hence the "const".
    return e_[subscript];
  }
};

template< class T >
tIdArray< T >::tIdArray(tList< T >& List) :
  e_(List.getSize())
{
  tListIter< T > Iter( List );
  T *c;
  for( c=Iter.FirstP(); !(Iter.AtEnd()); c=Iter.NextP() )
    e_[c->getID()] = c;
}

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tMesh.cpp"
#endif

#endif
