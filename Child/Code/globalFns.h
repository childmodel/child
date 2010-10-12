//-*-c++-*- 

/**************************************************************************/
/**
**  @file globalFns.h
**  @brief Header file for global functions used by tGrid and other
**         modules of CHILD.
**
**  This file contains geometric functions used for computing and
**  updating triangulations, a function for interpolating values on
**  a regular grid, and several functions for plane- and line-fitting
**  which can be used by interpolation procedures (such as the layer
**  interpolation routines in CHILD).
**
**  $Id: globalFns.h,v 1.19 2004-06-16 13:37:24 childcvs Exp $
*/
/**************************************************************************/

#ifndef GLOBALFNS_H
#define GLOBALFNS_H

#include <math.h>
#include "tArray/tArray.h"
#include "tArray/tArray2.h"
#include "MeshElements/meshElements.h"
#include "tLNode/tLNode.h"
#include "tPtrList/tPtrList.h"
#include "Predicates/predicates.h"
#include "tMatrix/tMatrix.h"
#include <vector>
using namespace std;
class tLNode;

// template< class tSubNode >
// class tMesh;

// Macros for max, min
/*#undef max
#undef min
#define max(a, b)               ((a) < (b) ? (b) : (a))
#define min(a, b)               ((a) > (b) ? (b) : (a))*/ // commented out Oct 09 to avoid macro conflict with <vector>


extern Predicates predicate; // object should be declared elsewhere, e.g. main

/******** Global Function Declarations **************************************/
tArray< double > UnitVector( tEdge const * );

inline tArray< double > UnitVector( double dx, double dy )
{
   const double mag = sqrt( dx * dx + dy * dy );
   return tArray< double >( dx/mag, dy/mag );
}

// tArray< double > UnitVector( double, double );

tArray< double > NewUnitVector( tEdge * );

double FindCosineAngle0_2_1( tArray< double > const &, tArray< double > const &,
                             tArray< double > const & );

int TriPasses( tArray< double > const &, tArray< double > const &,
               tArray< double > const &, tArray< double > const & );

bool PointsCCW( tArray< double > const &, tArray< double > const &, tArray< double > const & );

bool PointsCCW( tArray2< double > const &, tArray2< double > const &, tArray2< double > const & );

int Orientation( tArray< double > const &,
                 tArray< double > const &,
                 tArray< double > const & );

int NewTriCCW( tTriangle const * );

int InNewTri( tArray< double > const &, tTriangle const * );

template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &, tPtrListIter< tSubNode > & );

int Intersect( tEdge *, tEdge * );

tEdge* IntersectsAnyEdgeInList( tEdge*, tPtrList< tEdge >& );

double InterpSquareGrid( double, double, tMatrix< double > const &, int );

bool InBoundsOnMaskedGrid( double, double, tMatrix<double> const &, int );

tArray2< double > FindIntersectionCoords( tArray2< double > const&,
					  tArray2< double > const&,
					  tArray2< double > const&,
					  tArray2< double > const&);

//double timetrack;

//Returns z value at location x,y on the plane defined by the
//x-y coordinates in p0, p1, and p2, and their zvals in the zs array
double PlaneFit(double x, double y, tArray<double> const &p0,
                tArray<double>const & p1, tArray<double> const &p2, tArray<double> const &zs);

double PlaneFit(double tx, double ty, tTriangle* tri );

double LineFit(double x1, double y1, double x2, double y2, double nx);

double DistanceBW2Points(double x1, double y1, double x2, double y2 );

double DistanceToLine( double x2, double y2, double a, double b, double c );

double DistanceToLine( double x2, double y2, tNode const *p0, tNode const *p1 );

// global function to be used with veg. growth stuff:
// added by SL, 8/10
double Richards_Chapman_equ( const double t, const double Smax, 
			     const double decay, const double shape );

//
// New computational geometry global functions by SL, 9/2010
//

// DotProduct2D: returns dot product of edge from p0 to p1 (presumed to exist;
// doesn't actually use edge) and (unit) vector passed by reference:
// -STL, 7/2010
double DotProduct2D( tLNode *p0, tLNode *p1, tArray<double> &uv );

// function to find nodes along a line as long as they have the same public flag
void ScanLineNodesByPublicFlag( tLNode*, tPtrList<tLNode>& );

// function to find nodes along a line and within a certain distance
int ScanLineNodesByDistance( tLNode*, double, tPtrList<tLNode>& );

// similar to above, but only finds neighbors
void FindNbrsOnLine( tLNode* atPtr, tArray<double> &abc, 
		     vector<tLNode*> &BVec );

// find equation of line perpendicular to node's flow edge
void LinePerpendicularToFlowEdge( tLNode* atPtr, tArray<double> &abc );

// find equation of line perpendicular to vector and through (x,y)
inline void LinePerpendicularToVector( const double x, const double y, 
				       tArray<double> &dxy, 
				       tArray<double> &abc )
{
   abc[0] = dxy[0];
   abc[1] = dxy[1];
   abc[2] = -dxy[1] * y - dxy[0] * x;
}

tArray< double > RightUnitVector( tEdge* ePtr );
tArray< double > LeftUnitVector( tEdge* ePtr );

tLNode* ScanLineNextNode( tLNode* prvPtr, tLNode* atPtr, tArray<double> &abc );

void FillCrossSection( int iStream, double areaFrontal, double widthFrontal,
		       vector<tLNode*> &scanLineNode, 
		       tArray<double> &leftXYZ, tArray<double> &rightXYZ );





#endif
