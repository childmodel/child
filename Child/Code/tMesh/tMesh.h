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
**  $Id: tMesh.h,v 1.41 2003-04-07 17:01:03 childcvs Exp $
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
#define kRepairMesh 1
#define kNoRepair 0


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
public:
   tMesh();
   tMesh( tInputFile & );
   tMesh( tMesh * );
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
   void Print();
   /*makes edg, ccwedg structure from spokelists*/
   void MakeCCWEdges();
   void setVoronoiVertices();
   void CalcVoronoiEdgeLengths();
   void CalcVAreas();
   tTriangle *LocateTriangle( double, double );
   tTriangle *LocateNewTriangle( double, double );
   /*returns ptr to triangle which points to edge, or zero if none:*/ 
   tTriangle *TriWithEdgePtr( tEdge * );
   /*only routine needed to delete node; calls ExNode, RepairMesh:*/
   int DeleteNode( tListNode< tSubNode > *, int repairFlag=1 );
   int DeleteNode( tSubNode *, int repairFlag=1 );
   /*deletes spokes, *calls ExEdge, makes nbr ptr list:*/
   int ExtricateNode( tSubNode *, tPtrList< tSubNode > & );
   int DeleteEdge( tEdge * );
   /*calls ExTriangle; deals w/edge and compliment, deletes nbr tri(s):*/
   int ExtricateEdge( tEdge * );
   int DeleteTriangle( tTriangle * );
   /*"un-points" nbr triangles:*/
   int ExtricateTriangle( tTriangle * );
   /*complicated; fills in (any) hole defined by circular node ptr list:*/
   int RepairMesh( tPtrList< tSubNode > & );
   int AddEdgeAndMakeTriangle( tPtrList< tSubNode > &,
                               tPtrListIter< tSubNode > & );
   int AddEdgeAndMakeTriangle( tSubNode*, tSubNode*, tSubNode* );
   int MakeTriangle( tPtrList< tSubNode > &,
                     tPtrListIter< tSubNode > & );
   int AddEdge( tSubNode *, tSubNode *, tSubNode * );
   void CheckTrianglesAt( tSubNode * );
   void AddToList( tSubNode& );
   //add a node with referenced value/properties, update mesh connectivity
   tSubNode *AddNode( tSubNode &, int updatemesh = 0, double time = 0.0 );
   //add a generic node at the referenced coordinates
   tSubNode *AddNodeAt( tArray< double > &, double time = 0.0 );
   tSubNode* AttachNode( tSubNode*, tTriangle* );
   tMeshList<tEdge> * getEdgeList();
   tMeshList<tSubNode> * getNodeList();
   tList< tTriangle > * getTriList();
   tEdge *getEdgeComplement( tEdge * );
   /* Tests consistency of a user-defined mesh */
   void CheckMeshConsistency( int boundaryCheckFlag=1 );
   /* Updates mesh by comp'ing edg lengths & slopes & node Voronoi areas */ 
   void UpdateMesh();
   /* computes edge slopes as (Zorg-Zdest)/Length */
   //void CalcSlopes(); /* WHY is this commented out? */
   /*routines used to move points; MoveNodes is "master" function*/
   int CheckForFlip( tTriangle *, int, int );
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
   void MoveNodes( double time = 0.0, int interpFlag=1 );
   /*end moving routines*/

   void AddNodesAround( tSubNode *, double time=0.0 );  // Local mesh densify

#ifndef NDEBUG
   /*'dump' routines for debugging*/
   void DumpEdges();
   void DumpSpokes( tSubNode * );
   void DumpTriangles();
   void DumpNodes();
#endif
   
protected:
   int nnodes, nedges, ntri;       // # of nodes, edges, and tri's (obsolete?)
   tMeshList< tSubNode > nodeList; // list of nodes
   tMeshList< tEdge > edgeList;    // list of directed edges
   tList< tTriangle > triList;  // list of ptrs to triangles
   int miNextNodeID;                    // next ID for added triangle
   int miNextEdgID;                    // next ID for added triangle
   int miNextTriID;                    // next ID for added triangle
   long seed;                      // random seed
   int layerflag;                 // flag indicating whether nodes have layers
   tTriangle* mSearchOriginTriPtr; // ptr to tri. from which to start searches

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
