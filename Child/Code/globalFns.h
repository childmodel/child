/**************************************************************************\
**
**  globalFns.h: Header file for global functions used by tGrid and other
**               modules of CHILD.
**
**  This file contains geometric functions used for computing and
**  updating triangulations, a function for interpolating values on
**  a regular grid, and several functions for plane- and line-fitting
**  which can be used by interpolation procedures (such as the layer
**  interpolation routines in CHILD).
**
**  $Id: globalFns.h,v 1.3 1999-05-04 17:16:49 gtucker Exp $
\**************************************************************************/

#ifndef GLOBALFNS_H
#define GLOBALFNS_H

#include <iostream.h>
#include <math.h>
#ifdef __DECCXX
#include <macros.h> // for max(a,b), min(a,b), etc.
#else
#define max(a,b) ( (a>b) ? a : b )
#define min(a,b) ( (a<b) ? a : b )
#endif
#include "tArray/tArray.h"
#include "MeshElements/meshElements.h"
#include "tLNode/tLNode.h"
#include "tPtrList/tPtrList.h"
#include "Predicates/predicates.h"
#include "tMatrix/tMatrix.h"

extern Predicates predicate; // object should be declared elsewhere, e.g. main

/******** Global Function Declarations **************************************/
double ran3( long * );  // Random number generator from Numerical Recipes in C

tArray< double > UnitVector( tEdge* );

double FindCosineAngle0_2_1( tArray< double > &, tArray< double > &,
                             tArray< double > & );

int TriPasses( tArray< double > &, tArray< double > &,
               tArray< double > &, tArray< double > & );

int PointsCCW( tArray< double > &, tArray< double > &, tArray< double > & );

int NewTriCCW( tTriangle * );

int InNewTri( tArray< double > &, tTriangle * );

template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &, tPtrListIter< tSubNode > & );

template< class tSubNode >
int PointAndNext2Delaunay( tSubNode &, tPtrList< tSubNode > &,
                           tPtrListIter< tSubNode > & );

int Intersect( tEdge *, tEdge * );

tEdge* IntersectsAnyEdgeInList( tEdge*, tPtrList< tEdge >& );

double InterpSquareGrid( double, double, tMatrix< double >&, int );

tArray< double > FindIntersectionCoords( tArray< double >,
                                         tArray< double >,
                                         tArray< double >,
                                         tArray< double > );

template< class T > 
ostream &operator<<( ostream &, const tArray< T > & );

double timetrack;

//Returns z value at location x,y on the plane defined by the
//x-y coordinates in p0, p1, and p2, and their zvals in the zs array
double PlaneFit(double x, double y, tArray<double> p0,
                tArray<double> p1, tArray<double> p2, tArray<double> zs);

double LineFit(double x1, double y1, double x2, double y2, double nx);

double DistanceBW2Points(double x1, double y1, double x2, double y2 );


#endif
