/***************************************************************************\
**
**  tGrid.h: Header file for class tGrid
**
**  tGrid is the master class that handles the implementation of the
**  triangulated mesh. The class includes lists of the mesh elements
**  (nodes, triangles, and edges; see gridElements.h/.cpp), and provides
**  functionality to:
**    - read in or create meshes, either from scratch, from a list of
**      points, from a pre-existing set of triangulation files (e.g., a
**      previous run), or an Arc/Info grid DEM
**    - move, add, and/or delete nodes
**    - update Delaunay and Voronoi geometry
**
**  $Id: tMesh.h,v 1.20 1999-04-01 16:02:49 gtucker Exp $
\***************************************************************************/

#ifndef TGRID_H
#define TGRID_H

#include <assert.h>
#include <math.h>
#include "../Classes.h"
#include "../Definitions.h"
#include "../tArray/tArray.h"
#include "../tMatrix/tMatrix.h"
#include "../Mathutil/mathutil.h"
#include "../errors/errors.h"
#include "../tPtrList/tPtrList.h"
#include "../tList/tList.h"
#include "../tGridList/tGridList.h"
#include "../GridElements/gridElements.h"
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

/** class tGrid ************************************************************/
template< class tSubNode >
class tGrid
{
   friend class tListOutputData< tSubNode >;
public:
   tGrid();
   tGrid( tInputFile & );
   ~tGrid();
   void BatchAddNodes(); // quickly adds many nodes when starting w/ dense mesh
   void MakeGridFromScratch( tInputFile & );   // creates a new mesh
   void MakeGridFromInputData( tInputFile & ); // reads in an existing mesh
   void MakeGridFromPoints( tInputFile & );    // creates mesh from list of pts
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
   int MakeTriangle( tPtrList< tSubNode > &,
                     tPtrListIter< tSubNode > & );
   int AddEdge( tSubNode *, tSubNode *, tSubNode * );
   //add a node with referenced value/properties, update mesh connectivity
   tSubNode *AddNode( tSubNode &, int updatemesh = 1, double time = 0.0 );
   //add a generic node at the referenced coordinates
   tSubNode *AddNodeAt( tArray< double > &, double time = 0.0 );
   tGridList<tEdge> * getEdgeList();
   tGridList<tSubNode> * getNodeList();
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
   void MoveNodes( double time = 0.0 );
   /*end moving routines*/

#ifndef NDEBUG
   /*'dump' routines for debugging*/
   void DumpEdges();
   void DumpSpokes( tSubNode * );
   void DumpTriangles();
   void DumpNodes();
#endif
   
protected:
   int nnodes, nedges, ntri;       // # of nodes, edges, and tri's (obsolete?)
   tGridList< tSubNode > nodeList; // list of nodes
   tGridList< tEdge > edgeList;    // list of directed edges
   tList< tTriangle > triList;     // list of triangles
   long seed;                      // random seed
   bool layerflag;                 // flag indicating whether nodes have layers

};

#endif
