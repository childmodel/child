/***************************************************************************\
**
**  tGrid.h: Header file for class tGrid
**
**  $Id: tMesh.h,v 1.18 1999-02-04 19:26:28 nmgaspar Exp $
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
#include "../tListOutputData/tListOutputData.h"
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
   //tGrid( tListInputData< tSubNode > & );
   //tGrid( tListInputData< tSubNode > &, int );
   tGrid( tInputFile & );
   ~tGrid();
   void BatchAddNodes();
   void MakeGridFromScratch( tInputFile & );
   void MakeGridFromInputData( tInputFile & );
   void MakeGridFromPoints( tInputFile & );
	 void MakeRandomPointsFromArcGrid( tInputFile & );
   void MakeHexMeshFromArcGrid( tInputFile &infile );
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
   int nnodes, nedges, ntri;
   tGridList< tSubNode > nodeList;
   tGridList< tEdge > edgeList;
   tList< tTriangle > triList;
   long seed;
   bool layerflag;

};

#endif
