/***************************************************************************\
**
**  tGrid.h: Header file for class tGrid
**
**  $Id: tMesh.h,v 1.7 1998-03-21 21:50:58 gtucker Exp $
\***************************************************************************/

#ifndef TGRID_H
#define TGRID_H

#include <assert.h>
#include <math.h>
#include "../Classes.h"
#include "../Definitions.h"
#include "../tLNode/tLNode.h"
#include "../Mathutil/mathutil.h"
#include "../errors/errors.h"
#include "../tList/tList.h"
#include "../tGridList/tGridList.h"
#include "../GridElements/gridElements.h"
#include "../tPtrList/tPtrList.h"
#include "../tListInputData/tListInputData.h"
#include "../tListOutputData/tListOutputData.h"

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
    void MakeGridFromScratch( tInputFile & );
    void MakeGridFromInputData( tInputFile & );
    void Print();
    /*makes edg, ccwedg structure from spokelists*/
    void MakeCCWEdges();
    void SetVoronoiVertices();
    void CalcVoronoiEdgeLengths();
    void CalcVAreas();
    tTriangle *LocateTriangle( double, double );
    tTriangle *LocateNewTriangle( double, double );
     /*returns ptr to triangle which points to edge, or zero if none:*/ 
    tTriangle *TriWithEdgePtr( tEdge * );
     /*only routine needed to delete node; calls ExNode, RepairMesh:*/
    int DeleteNode( tListNode< tSubNode > * );
    int DeleteNode( tSubNode * );
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
    int AddNode( tSubNode & );
     //add a generic node at the referenced coordinates
    tSubNode *AddNodeAt( tArray< double > & );
    tGridList<tEdge> * GetEdgeList();
    tGridList<tSubNode> * GetNodeList();
    tList< tTriangle > * GetTriList();
    tEdge *getEdgeCompliment( tEdge * );
    void CheckMeshConsistency();  // Tests consistency of a user-defined mesh
     /*Updates mesh by computing edge lengths & slopes & node Voronoi areas */ 
    void UpdateMesh();
     /* computes edge slopes as (Zorg-Zdest)/Length */
   //void CalcSlopes();
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
    void MoveNodes();
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
};

#endif
