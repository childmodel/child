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
**  $Id: tMesh.h,v 1.82 2008/07/07 16:18:58 childcvs Exp $
*/
/***************************************************************************/

#ifndef TMESH_H
#define TMESH_H

#include <assert.h>
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

/** @class tIdArray
    @brief Lookup table per Id for a tList
*/
template< class T, class LN  = tListNodeBasic<T> > class tIdArray;

template< class T, class LN >
class tIdArray
{
  tIdArray();
  tArray< T* > e_;
public:
  tIdArray(tList< T, LN >& List);
  T* operator[]( int subscript ) const {
    // return a value and not a reference, hence the "const".
    return e_[subscript];
  }
};

template< class T, class LN >
tIdArray< T, LN >::tIdArray(tList< T, LN >& List) :
  e_(List.getSize())
{
  tListIter< T, LN > Iter( List );
  T *c;
  for( c=Iter.FirstP(); !(Iter.AtEnd()); c=Iter.NextP() )
    e_[c->getID()] = c;
}

/***************************/
/****** Defined types ******/
/***************************/

typedef enum {
  FLIP_NOT_NEEDED,
  FLIP_NEEDED,
  FLIP_ERROR,
  FLIP_DONE,
  FLIP_NOT_ALLOWED
} flipStatus_t;

/* codes passed to DeleteNode/AddNode */
typedef enum {
  kNoRepairMesh,
  kRepairMesh
} kRepairMesh_t;
typedef enum {
  kNoUpdateMesh,
  kUpdateMesh
} kUpdateMesh_t;
typedef enum {
  kNoFlip,
  kFlip
} kFlip_t;


/****************************/
/**** Class Declarations ****/
/****************************/

class ParamMMFS_t;

/** @class tMesh
 */
template< class tSubNode >
class tMesh
{
   tMesh(const tMesh&);
   tMesh& operator=(const tMesh&);
   tMesh();

   void MakePointBoundary( const ParamMMFS_t &, const tInputFile &,
			   tPtrList< tSubNode > &, tRand &);
   void MakePointInterior( const ParamMMFS_t &, const tInputFile &,
			   bool makeMesh, tRand &);
   void BuildDelaunayMeshTipper();
   void MeshDensification( const tInputFile & );
   void MakeDelaunay( tPtrList< tTriangle > &, double time );
   void SplitNonFlippableEdge( tPtrList< tEdge > &, double time );
public:
   // list types
   typedef tMeshList< tSubNode, tListNodeListable< tSubNode > > nodeList_t;
   typedef tMeshList< tEdge, tListNodeListable< tEdge > > edgeList_t;
   typedef tList< tTriangle, tListNodeListable< tTriangle > > triList_t;

   // iterators types
   typedef tMeshListIter< tSubNode, tListNodeListable< tSubNode > > nodeListIter_t;
   typedef tMeshListIter< tEdge, tListNodeListable< tEdge > > edgeListIter_t;
   typedef tListIter< tTriangle, tListNodeListable< tTriangle > > triListIter_t;

   // tListNode types
   typedef tListNodeListable< tSubNode > nodeListNode_t;
   typedef tListNodeListable< tEdge > edgeListNode_t;
   typedef tListNodeListable< tTriangle > triListNode_t;

   // tIdArray
   typedef tIdArray< tSubNode, tListNodeListable< tSubNode > > tIdArrayNode_t;
   typedef tIdArray< tEdge, tListNodeListable< tEdge > > tIdArrayEdge_t;
   typedef tIdArray< tTriangle, tListNodeListable< tTriangle > > tIdArrayTri_t;

   tMesh( const tInputFile &, bool checkMeshConsistency );
   tMesh( tMesh const * );
   ~tMesh();
   void BatchAddNodes(); // quickly adds many nodes when starting w/ dense mesh
   void MakeMeshFromScratch( const tInputFile &, tRand & ); // creates a new mesh
   void MakeMeshFromScratchTipper( const tInputFile &, tRand & );   // creates a new mesh
   void MakeMeshFromInputData( const tInputFile & ); // reads in an existing mesh
   void MakeMeshFromPoints( const tInputFile & );    // creates mesh from list of pts
   void MakeMeshFromPointsTipper( const tInputFile & ); // creates mesh from list of pts
   void MakeRandomPointsFromArcGrid( const tInputFile & ); // mesh from arc (rand)
   void MakeHexMeshFromArcGrid( const tInputFile & );// mesh from arc (hex)
   void MakeLayersFromInputData( const tInputFile & );
   void Print() const;
   void setVoronoiVertices();
   void CalcVoronoiEdgeLengths();
   void CalcVAreas();
   tTriangle *LocateTriangle( double, double, bool useFuturePosn=false );
   tTriangle *LocateNewTriangle( double, double );
   /*returns ptr to triangle which points to edge, or zero if none:*/
   tTriangle *TriWithEdgePtr( tEdge * ) const;
   /*only routine needed to delete node; calls ExNode, RepairMesh:*/
   int DeleteNode( nodeListNode_t *, kRepairMesh_t repairFlag=kRepairMesh,
		   kUpdateMesh_t updateFlag=kUpdateMesh,
		   bool allowMobileDeletion = false );
   int DeleteNode( tSubNode const *, kRepairMesh_t repairFlag=kRepairMesh,
		   kUpdateMesh_t updateFlag=kUpdateMesh,
		   bool allowMobileDeletion = false );
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
   void CheckTrianglesAt( tSubNode *, double time  );
   tSubNode *AddToList( tSubNode const& );
   void RemoveFromList( tSubNode * ); // removes the last node added
   //add a node with referenced value/properties, update mesh connectivity
   tSubNode *AddNode( tSubNode &, kUpdateMesh_t updatemesh = kNoUpdateMesh,
		      double time = 0.0, kFlip_t flip = kFlip );
   tSubNode* InsertNode( tSubNode*, double );
   //add a generic node at the referenced coordinates
   tSubNode *AddNodeAt( tArray< double > &, double time = 0.0 );
   tSubNode* AttachNode( tSubNode*, tTriangle* );
   edgeList_t * getEdgeList() { return &edgeList; }
   nodeList_t * getNodeList() { return &nodeList; }
   triList_t * getTriList() { return &triList; }
   tEdge *getEdgeComplement( tEdge * ) const;
   /* Tests consistency of a user-defined mesh */
   void CheckMeshConsistency( bool boundaryCheckFlag=true );
   /* Updates mesh by comp'ing edg lengths & slopes & node Voronoi areas */
   void UpdateMesh( bool checkMeshConsistency = true );
   /* computes edge slopes as (Zorg-Zdest)/Length */
   //void CalcSlopes(); /* WHY is this commented out? */
   /*routines used to move points; MoveNodes is "master" function*/
   flipStatus_t CheckForFlip( tTriangle *, int, bool, bool useFuturePosn=true );
   bool FlipEdge( tTriangle *, tTriangle *, int, int, bool useFuturePosn=false );
   tEdge * IntersectsAnyEdge( tEdge * );
   //CheckTriEdgeIntersect and CheckLocallyDelaunay are the essential functions
   //for point moving and are called from MoveNodes; in general, they should be
   //after the three above functions and in the order listed below.
   void CheckTriEdgeIntersect();
   void CheckLocallyDelaunay( double time );
   //once 'newx' and 'newy' are set, this is the only function one needs to call
   //to execute the move; maybe should have a separate function for doing things
   //peculiar to channels, but now this is the only one.
   void MoveNodes( double time = 0.0, bool interpFlag=true );
   /*end moving routines*/

   void AddNodesAround( tSubNode *, double time=0.0 );  // Local mesh densify

   static bool IDTooLarge(int, int);
   void ResetNodeIDIfNecessary();
   void ResetEdgeIDIfNecessary();
   void ResetTriangleIDIfNecessary();
   void ResetNodeID(); // reset node IDs in list order
   void ResetEdgeID(); // reset edge IDs in list order
   void ResetTriangleID(); // reset triangle IDs in list order
   void RenumberIDCanonically(); // reset IDs in canonical order
   void SetmiNextNodeID(int);
   void SetmiNextEdgID(int);
   void SetmiNextTriID(int);
   // find triangles between one node and the next, not connected by an edge
   tPtrList< tTriangle > InterveningTriangles( tNode*, tNode* ) const;
   void ForceFlow( tSubNode*, tSubNode*, double );

#ifndef NDEBUG
   /*'dump' routines for debugging*/
   void DumpEdges();
   void DumpSpokes( tSubNode * ) const;
   void DumpTriangles();
   void DumpNodes();
#endif


private:
   static int orderRNode(const void*, const void*);
   static int orderREdge(const void*, const void*);
   static int orderRTriangle(const void*, const void*);

protected:
   nodeList_t nodeList; // list of nodes
   edgeList_t edgeList;    // list of directed edges
   triList_t triList;  // list of ptrs to triangles

   tTriangle* mSearchOriginTriPtr; // ptr to tri. from which to start searches
   int nnodes, nedges, ntri;       // # of nodes, edges, and tri's (obsolete?)
   int miNextNodeID;                   // next ID for added node
   int miNextPermNodeID;               // next Permanent Node ID
   int miNextEdgID;                    // next ID for added edge
   int miNextTriID;                    // next ID for added triangle
   bool layerflag;                 // flag indicating whether nodes have layers
   bool runCheckMeshConsistency;    // shall we run the tests ?
};

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tMesh.cpp"
#endif

#endif
