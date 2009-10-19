/**************************************************************************/
/**
**  @file globalFns.cpp
**  @brief Global functions used by tGrid and other modules of
**         CHILD (see globalFns.h).
**
**  $Id: globalFns.cpp,v 1.20 2004-06-16 13:37:24 childcvs Exp $
*/
/**************************************************************************/

#include "globalFns.h"
#include <iostream>
#include <algorithm>  // added Oct 09 to replace max and min macros
  using std::max;  // ditto
  using std::min;  // ditto

/**************************************************************************/
/**
**  UnitVector
**
**  global function to return an array with the magnitudes of the x and y
**  components of the unit vector parallel to the pointed-to edge.
**  Created: SL 10/98
**
*/
/**************************************************************************/
tArray< double > UnitVector( tEdge const * ePtr )
{
   assert( ePtr != 0 );
   tArray< double > oxy( ePtr->getOriginPtr()->get2DCoords() );
   tArray< double > dxy( ePtr->getDestinationPtr()->get2DCoords() );
   const double dx = dxy[0] - oxy[0];
   const double dy = dxy[1] - oxy[1];
   const double mag = sqrt( dx * dx + dy * dy );
   return
     tArray< double >( dx/mag, dy/mag );
}

tArray< double > UnitVector( double dx, double dy )
{
   const double mag = sqrt( dx * dx + dy * dy );
   return tArray< double >( dx/mag, dy/mag );
}

/**************************************************************************/
/**
**  NewUnitVector
**
**  global function to return an array with the magnitudes of the x and y
**  components of the unit vector parallel to the pointed-to edge (with
**  endpoints at "new" coordinates).
**  Modified from UnitVector.
**  Created: SL 8/2003
**
*/
/**************************************************************************/
tArray< double > NewUnitVector( tEdge * ePtr )
{
   assert( ePtr != 0 );
   tArray< double > oxy( ePtr->getOriginPtrNC()->FuturePosn() );
   tArray< double > dxy( ePtr->getDestinationPtrNC()->FuturePosn() );
   const double dx = dxy[0] - oxy[0];
   const double dy = dxy[1] - oxy[1];
   const double mag = sqrt( dx * dx + dy * dy );
   return
     tArray< double >( dx/mag, dy/mag );
}

/**************************************************************************/
/**
**  FindCosineAngle0_2_1
**
*/
/**************************************************************************/
double FindCosineAngle0_2_1( tArray< double > const &p0,
                             tArray< double > const &p1,
                             tArray< double > const &p2 )
{
   assert( (&p0 != 0) && (&p1 != 0) && (&p1 != 0) );
   double dx0, dx1, dy0, dy1;
   double dotp, magp;
   dx0 = p0[0] - p2[0];
   dx1 = p1[0] - p2[0];
   dy0 = p0[1] - p2[1];
   dy1 = p1[1] - p2[1];
   dotp = dx0 * dx1 + dy0 * dy1;
   magp = sqrt( dx0 * dx0 + dy0 * dy0 ) * sqrt( dx1 * dx1 + dy1 * dy1 );
   return dotp / magp;
}


/***************************************************************************/
/**
**  TriPasses
**
**  Determines whether test point violates "Delaunay-ness"
**  of other three; i.e., does the triangle 'pass the test' against the
**  other?
**
*/
/***************************************************************************/
int TriPasses( tArray< double > const &ptest,
               tArray< double > const &p0,
               tArray< double > const &p1,
               tArray< double > const &p2 )
{
   if (0) //DEBUG
     std::cout << "TriPasses? ";
#if 1
   double ans = predicate.incircle( p0.getArrayPtr(), p1.getArrayPtr(),
				    p2.getArrayPtr(), ptest.getArrayPtr() );
   return ans > 0. ? 0 : 1;
#else
   double dx0, dx1, dy0, dy1;
   double crossp,      // cross-product
       dotp,           // dot-product
       angle0_2_1,     // angle p0-p2-p1
       angle0_test_1;  // angle p0-ptest-p1

   // Find angle p0-p2-p1
   dx0 = p0[0] - p2[0];
   dx1 = p1[0] - p2[0];
   dy0 = p0[1] - p2[1];
   dy1 = p1[1] - p2[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   dotp = dx0 * dx1 + dy0 * dy1;
   angle0_2_1 = atan2( crossp, dotp );

   // Find angle p0-ptest-p1
   dx0 = p0[0] - ptest[0];
   dx1 = p1[0] - ptest[0];
   dy0 = p0[1] - ptest[1];
   dy1 = p1[1] - ptest[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   dotp = dx0 * dx1 + dy0 * dy1;
   angle0_test_1 = atan2( crossp, dotp );

   // Compare and return the result
   if( angle0_2_1 < angle0_test_1 )
   {
        //std::cout << "Yes" << std::endl;
      return 0;
   }
   else
   {
        //std::cout << "No" << std::endl;
      return 1;
   }
#endif
}


/***************************************************************************/
/**
**  PointsCCW
**
**  Determines whether the points p0, p1, and p2 are in counter-clockwise
**  order.
**
**  Inputs: p0, p1, p2 -- 2-element arrays containing x and y coords
**  Returns: 1 if CCW, 0 if not
**
*/
/***************************************************************************/
bool PointsCCW( tArray< double > const &p0,
		tArray< double > const &p1,
		tArray< double > const &p2 )
{
   if (0) //DEBUG
     std::cout << "PointsCCW? ";

   if( p0 == p1 || p0 == p2 || p1 == p2 )
       return false;

   const double* a0 = p0.getArrayPtr();
   const double* a1 = p1.getArrayPtr();
   const double* a2 = p2.getArrayPtr();
   // call exact arithmetic predicate:
   return ( predicate.orient2d( a0, a1, a2 ) > 0 );
}


bool PointsCCW( tArray2< double > const &p0,
		tArray2< double > const &p1,
		tArray2< double > const &p2 )
{
   if (0) //DEBUG
     std::cout << "PointsCCW? ";

   if( p0 == p1 || p0 == p2 || p1 == p2 )
     return false;

   const double* a0 = p0.getArrayPtr();
   const double* a1 = p1.getArrayPtr();
   const double* a2 = p2.getArrayPtr();
   // call exact arithmetic predicate:
   return ( predicate.orient2d( a0, a1, a2 ) > 0 );
}
/***************************************************************************/
/**
**  Orientation
**
**  global function; determines whether points are counter-clockwise (1),
**  clockwise (-1) or colinear (0):
**
**  Inputs: p0, p1, p2 -- 2-element arrays containing x and y coords
**
*/
/***************************************************************************/
int Orientation( tArray< double > const &p0,
                 tArray< double > const &p1,
                 tArray< double > const &p2 )
{
   const double* a0 = p0.getArrayPtr();
   const double* a1 = p1.getArrayPtr();
   const double* a2 = p2.getArrayPtr();
   //NEW: call exact arithmetic predicate:
   const double ans = predicate.orient2d( a0, a1, a2 );
   if( ans > 0. ) return 1;
   if( ans < 0. ) return -1;
   return 0;
}

/***************************************************************************/
/**
**  NewTriCCW
**
**  Determines whether the "newx" and "newy" coordinates of the 3 points
**  in the triangle are counter-clockwise.
**
**  Inputs: ct -- pointer to the triangle to be tested
**  Returns: 1 if CCW, 0 if not
**
**  TODO: use general "moving" indicator, or dispense w/ the mdr test
*/
/***************************************************************************/
int NewTriCCW( tTriangle const *ct )
{
   assert( ct != 0 );
   tLNode *cn;

   cn = static_cast<tLNode *>(ct->pPtr(0));
   tArray< double > p0( cn->get2DCoords() );

   if( cn->Meanders() ) p0 = cn->getNew2DCoords();
   cn = static_cast<tLNode *>(ct->pPtr(1));
   tArray< double > p1( cn->get2DCoords() );
   if( cn->Meanders() ) p1 = cn->getNew2DCoords();
   cn = static_cast<tLNode *>(ct->pPtr(2));
   tArray< double > p2( cn->get2DCoords() );
   if( cn->Meanders() ) p2 = cn->getNew2DCoords();
   if( PointsCCW( p0, p1, p2 ) ) return 1;
   else return 0;
}


/***************************************************************************/
/**
**  InNewTri
**
**  Determines whether the point xy lies within the triangle formed by
**  the "newx" and "newy" coordinates of the 3 points
**  in the triangle.
**
**  Inputs: xy -- reference to 2-element array w/ x and y coords of point
**          ct -- pointer to the triangle to be tested
**  Returns: 1 if point is in the "new" triangle, 0 if not
**
**  TODO: use general "moving" indicator, or dispense w/ the mdr test
*/
/***************************************************************************/
int InNewTri( tArray< double > const &xy, tTriangle const *ct )
{
   int j;
   tLNode *vtx;
   tArray< double > xy1, xy2;
   for( j=0; j<3; j++ )
   {
      vtx = static_cast<tLNode *>(ct->pPtr(j));
      if( vtx->Meanders() ) xy1 = vtx->getNew2DCoords();
      else xy1 = vtx->get2DCoords();
      vtx = static_cast<tLNode *>(ct->pPtr( (j+1)%3 ));
      if( vtx->Meanders() ) xy2 = vtx->getNew2DCoords();
      else xy2 = vtx->get2DCoords();
      if ( ( (xy1[1] - xy[1]) * (xy2[0] - xy[0]) ) >
           ( (xy1[0] - xy[0]) * (xy2[1] - xy[1])) )
          break;
   }
   if( j == 3) return 1;
   else return 0;
}


/*****************************************************************************/
/**
**      Intersect
**
**      Tests for the intersection of two edges, using newx, newy rather than
**      x, y of the endpoints.
**
**      Inputs: ae, be -- ptrs to the two edges to be tested
**      Returns: 1 if they intersect, 0 if they don't
**      Called by: 
**      Calls:  
**
*/
/*****************************************************************************/
int Intersect( tEdge * ae, tEdge * be )
{
   if (0) //DEBUG
     std::cout << "Intersect(...)..." << std::endl;
   tLNode * lnode;
   
   if( !ae || !be )
   {
      std::cout<<"Intersect: Warning: invalid edge(s)"<<std::endl;
      return( 0 );
   }
   if( !ae->getOriginPtr() || !ae->getDestinationPtr() ||
       !be->getOriginPtr() || !be->getOriginPtr() )
   {
      std::cout<<"Intersect: Warning: invalid org or dest"<<std::endl;
      return( 0 );
   }
   lnode = static_cast<tLNode *>(ae->getOriginPtrNC());
   tArray< double > A( lnode->get2DCoords() );
   if( lnode->Meanders() ) A = lnode->getNew2DCoords();
   lnode = static_cast<tLNode *>(ae->getDestinationPtrNC());
   tArray< double > B( lnode->get2DCoords() );
   if( lnode->Meanders() ) B = lnode->getNew2DCoords();
   lnode = static_cast<tLNode *>(be->getOriginPtrNC());
   tArray< double > C( lnode->get2DCoords() );
   if( lnode->Meanders() ) C = lnode->getNew2DCoords();
   lnode = static_cast<tLNode *>(be->getDestinationPtrNC());
   tArray< double > D( lnode->get2DCoords() );
   if( lnode->Meanders() ) D = lnode->getNew2DCoords();
   /* algorithm from comp.graphics.algorithms FAQ; cited were: O'Rourke, pp.249-51
      and Gems III, pp. 199-202, "Faster Line Segment Intersection":

      For points in 2d real space, A, B, C, D, forming line segments AB, CD,

      AB = A+r(B-A), r in [0,1]; CD = C+s(D-C), s in [0,1];

      if AB and CD intersect, then A+r(B-A) = C+s(D-C), or

      XA+r(XB-XA) = XC+s(XD-XC) and YA+r(YB-YA) = YC+s(YD-YC);

      solving for r and s yields:

          (YA-YC)(XD-XC)-(XA-XC)(YD-YC)
      r = -----------------------------
          (XB-XA)(YD-YC)-(YB-YA)(XD-XC)

          (YA-YC)(XB-XA)-(XA-XC)(YB-YA)
      s = -----------------------------
          (XB-XA)(YD-YC)-(YB-YA)(XD-XC)

      for 0<=r<=1 and 0<=s<=1, the segments intersect;
      for r<0 or r>1 or s<0 or s>1, the segments do not intersect;
      if the denomenator of r = 0, then the segments are parallel;
      if the numerator of r is also = 0, then the segments are coincident.
      Intersection at the endpoints is r = 0 or 1 and s = 0 or 1; since
      we want to allow intersecting endpoints, we'll allow these equivalences.*/

   double rnumsign;
   double snumsign;
   
   double denomsign =
       predicate.DifferenceOfProductsOfDifferences( B[0], A[0],
                                                      D[1], C[1],
                                                      B[1], A[1],
                                                      D[0], C[0] );
   if( denomsign == 0.0 ) // segments parallel
   {
      rnumsign =
          predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         D[0], C[0],
                                                         A[0], C[0],
                                                         D[1], C[1] );
      if( rnumsign != 0.0 ) return 0;
      else
      {
         if( A[0] != B[0] ) //segments not vertical
         {
            if( max( A[0], B[0] ) <= min( C[0], D[0] ) ) return 0; // segments
            if( max( C[0], D[0] ) <= min( A[0], B[0] ) ) return 0; // do not overlap
            return 1; // segments overlap
         }
         else // vertical
         {
            if( max( A[1], B[1] ) <= min( C[1], D[1] ) ) return 0; // segments
            if( max( C[1], D[1] ) <= min( A[1], B[1] ) ) return 0; // do not overlap
            return 1; // segments overlap
         }
      }
   }
   else if( denomsign > 0.0 ) // not parallel
   {
      if( A[0] == C[0] && A[1] == C[1] ) return 0; // segments
      if( A[0] == D[0] && A[1] == D[1] ) return 0; // intersect 
      if( B[0] == C[0] && B[1] == C[1] ) return 0; // at
      if( B[0] == D[0] && B[1] == D[1] ) return 0; // endpoints
      
      rnumsign =
          predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         D[0], C[0],
                                                         A[0], C[0],
                                                         D[1], C[1] );
      if( rnumsign < 0.0 ) return 0; // r < 0
      snumsign =
          predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         B[0], A[0],
                                                         A[0], C[0],
                                                         B[1], A[1] );
      if( snumsign < 0.0 ) return 0; // s < 0
   }
   else
   {
      if( A[0] == C[0] && A[1] == C[1] ) return 0; // segments
      if( A[0] == D[0] && A[1] == D[1] ) return 0; // intersect 
      if( B[0] == C[0] && B[1] == C[1] ) return 0; // at
      if( B[0] == D[0] && B[1] == D[1] ) return 0; // endpoints
      
      rnumsign =
          predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         D[0], C[0],
                                                         A[0], C[0],
                                                         D[1], C[1] );
      if( rnumsign > 0.0 ) return 0; // r < 0
      snumsign =
          predicate.DifferenceOfProductsOfDifferences( A[1], C[1],
                                                         B[0], A[0],
                                                         A[0], C[0],
                                                         B[1], A[1] );
      if( snumsign > 0.0 ) return 0; // s < 0
   }
   // reverse directions of segments so we can still detect
   // intersection by the signs of the numerators and denomenators:
   // A -> B, B->A, C->D, D->C
   denomsign = 
       predicate.DifferenceOfProductsOfDifferences( A[0], B[0],
                                                      C[1], D[1],
                                                      A[1], B[1],
                                                      C[0], D[0] );
   if( denomsign > 0.0 )
   {
      rnumsign =
          predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
                                                         C[0], D[0],
                                                         B[0], D[0],
                                                         C[1], D[1] );
      if( rnumsign < 0.0 ) return 0; // r > 1
      snumsign =
          predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
                                                         A[0], B[0],
                                                         B[0], D[0],
                                                         A[1], B[1] );
      if( snumsign < 0.0 ) return 0; // s > 1
   }
   else if( denomsign < 0.0 )
   {
      rnumsign =
          predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
                                                         C[0], D[0],
                                                         B[0], D[0],
                                                         C[1], D[1] );
      if( rnumsign > 0.0 ) return 0; // r > 1
      snumsign =
          predicate.DifferenceOfProductsOfDifferences( B[1], D[1],
                                                         A[0], B[0],
                                                         B[0], D[0],
                                                         A[1], B[1] );
      if( snumsign > 0.0 ) return 0; // s > 1
   }
   //else segments must intersect other than at endpoints:
   return 1;
   
      
}


/*****************************************************************************/
/**
**      IntersectsAnyEdgeInList
**
**         Returns the first edge in the list (passed by
**         pointer) which intersects "edge" or NULL if "edge" intersects no
**         other edges
**
**      Data members updated: 
**      Called by: tMesh::BatchAddNodes()
**      Calls: Intersect
**      Created: SL 10/98
**      Based on tMesh::IntersectsAnyEdge()
**
*/
/*****************************************************************************/
tEdge* IntersectsAnyEdgeInList( tEdge* edge, tPtrList< tEdge >& edglistRef )
{
   if (0) //DEBUG
     std::cout << "IntersectsAnyEdge( tEdge * edge )..." << std::endl;
   tEdge * ce;
   tPtrListIter< tEdge > edgIter( edglistRef );
   if( !edge )
   {
      std::cout<<"IntersectsAnyEdge: Warning: invalid edge"<<std::endl;
      return 0;
   }
   
   // check all edges in list; implication is that the list contains
   // only edges that need to be checked:
   for( ce = edgIter.FirstP(); !(edgIter.AtEnd()); ce = edgIter.NextP() )
   {
      if( edge != ce ) // shouldn't be the same, but just to make sure...
      {
         if( Intersect( edge, ce ) ) return( ce );
      }
      
   }
   assert( edgIter.AtEnd() );
   return 0;
}


/**************************************************************************/
/**
**  InterpSquareGrid
**
**   Was part of tMesh::MakeRandomPointsFromArcGrid();
**   pulled it out so it can be used generally.
**		Use Tetzlaff & Harbaugh, '88, grid interpolation:
**			z(z1, z2, z3, z4, X, Y) = z1 * XY + z2 * (1-X)Y + z3 * X(1-Y) 
**												+ z4 * (1-X)(1-Y),  0<=X,Y<=1
**
**	 Assumes: passed coordinates are normalized to grid spacing = 1.
**   Created: 10/98 SL
**   Modified: to reflect difference between y (positive "up") and
**      j (pos. "down")
**
**
*/
/**************************************************************************/
double InterpSquareGrid( double xgen, double ygen, tMatrix< double > const & elev,
                         int nodata )
{
   int nodatacount = 0;
   int nx = elev.getNumCols();
   int ny = elev.getNumRows();
   
   int ix = static_cast<int>(xgen);  //                      z3----z4
   int iy = static_cast<int>(ygen);  //                      |     |
   double xrem = xgen - static_cast<double>(ix); //          |     |
   double yrem = ygen - static_cast<double>(iy); //          z1----z2
   double yrows = elev.getNumRows();
   int i = ix;
   int j = static_cast<int>(yrows - 1.0 - ygen);
   
   double z3 = elev( j, i );
   double z4 = ( i < nx - 1 ) ? elev( j, i + 1 ) : nodata;
   double z1 = ( j < ny - 1 ) ? elev( j + 1, i ) : nodata;
   double z2 = ( i < nx - 1 && j < ny - 1 ) ? elev( j + 1, i + 1 ) : nodata;

   if( z1 == nodata ) nodatacount++;
   if( z2 == nodata ) nodatacount += 2;
   if( z3 == nodata ) nodatacount += 4;
   if( z4 == nodata ) nodatacount += 8;
   if( nodatacount == 15 ) return nodata;					// all 4 corners have no data
   if( nodatacount == 0 )					// all 4 corners have data
       return z2 + (z4 - z2) * xrem + (z1 - z2) * yrem 
           - (z1 + z4 - z2 - z3) * xrem * yrem;
   if( nodatacount == 9 &&					// z1 and z4 have no data
       ( ( xrem <= 0.5 && yrem >= 0.5 ) || ( xrem >= 0.5 && yrem <= 0.5 ) ) )
       return z2 + 0.5 * ( z3 - z2 ) * ( xrem + yrem ); //z1 = z4 = (z2 + z3) / 2.0;
   if( nodatacount == 6 &&					// z2 and z3 have no data
       ( ( xrem <= 0.5 && yrem <= 0.5 ) || ( xrem >= 0.5 && yrem >= 0.5 ) ) )
       return 0.5 * ( z1 + z4 + ( z4 - z1 ) * ( xrem - yrem ) ); // z2 = z3 = (z1 + z4) / 2.0;
   if( nodatacount == 1 &&					// z1 has no data
       !( xrem < 0.5 && yrem < 0.5 ) )
       return z2 + ( z4 - z2 ) * xrem + 0.5 * ( z3 - z2 ) * yrem
           - ( z4 - 0.5 * ( z2 + z3 ) ) * xrem * yrem; // z1 = (z2 + z3) / 2.0;
   if( nodatacount == 2 &&					// z2 has no data
       !( xrem > 0.5 && yrem < 0.5 ) )
       return 0.5 * ( z1 + z4 + ( z4 - z1 ) * ( xrem - yrem ) )
           - ( 0.5 * ( z1 + z4 ) - z3 ) * xrem * yrem; //z2 = (z1 + z4) / 2.0;
   if( nodatacount == 4 &&					// z3 has no data
       !( xrem < 0.5 && yrem > 0.5 ) )
       return z2 + ( z4 - z2 ) * xrem + ( z1 - z2 ) * yrem
           - ( 0.5 * ( z1 + z4 ) - z2 ) * xrem * yrem; //z3 = (z1 + z4) / 2.0;
   if( nodatacount == 8 &&				  // z4 has no data
       !( xrem > 0.5 && yrem > 0.5 ) )
       return z2 + 0.5 * ( z3 - z2 ) * xrem + ( z1 - z2 ) * yrem 
           - ( z1 - 0.5 * ( z2 + z3 ) ) * xrem * yrem; //z4 = (z2 + z3) / 2.0;
   if( nodatacount == 3 && yrem >= 0.5 )        // z1 & z2 have no data
       return z4 + (z3 - z4) * yrem;              // => z1 = z3; z2 = z4;
   if( nodatacount == 10 && xrem <= 0.5 )       // z2 & z4 have no data
       return z1 + (z3 - z1) * xrem;              // => z2 = z1;  z4 = z3;
   if( nodatacount == 12 && yrem <= 0.5 )       // z3 & z4 have no data
       return z2 + (z1 - z2) * yrem;              // => z3 = z1; z4 = z2;
   if( nodatacount == 5 && xrem >= 0.5 )        // z1 & z3 have no data
       return z2 + (z4 - z2) * xrem;              // => z1 = z2; z3 = z4;
   if( nodatacount == 14 && xrem <= 0.5 && yrem <= 0.5 ) // z1 has data
       return z1;
   if( nodatacount == 13 && xrem >= 0.5 && yrem <= 0.5 ) // z2 has data
       return z2;
   if( nodatacount == 11 && xrem <= 0.5 && yrem >= 0.5 ) // z3 has data
       return z3;
   if( nodatacount == 7 && xrem >= 0.5 && yrem >= 0.5 )  // z4 has data
       return z4;
   return nodata;
}


/**********************************************************************/
/** 
**  PlaneFit 
** 
**  A plane is fit given the x,y,z coordinates of three points.
**  Returns the z value on this plane at the x and y location
**  that is sent to it.
**
**  tx, ty are the location where you want to know z
**  p0, p1, p2 are arrays containing x and y values (p0[0]=x; p0[1]=y)
**  zs contain the z values at p0, p1, and p2, respectively in the array
**
**  Created 2/1999 ng
**
*/
/**********************************************************************/
double PlaneFit(double x, double y, tArray<double> const &p0,
                tArray<double> const &p1, tArray<double> const &p2, tArray<double> const &zs)
{
   double a, b, c;
   double y0, y1, y2, x0, x1, x2, z0, z1, z2;
   y0=p0[1];
   y1=p1[1];
   y2=p2[1];
   x0=p0[0];
   x1=p1[0];
   x2=p2[0];
   z0=zs[0];
   z1=zs[1];
   z2=zs[2];

   //std::cout<<"PlaneFit"<<std::endl;
   //std::cout<<"x0 "<<x0<<" x1 "<<x1<<" x2 "<<x2<<std::endl;
   //std::cout<<"y0 "<<y0<<" y1 "<<y1<<" y2 "<<y2<<std::endl;
   //std::cout<<"z0 "<<z0<<" z1 "<<z1<<" z2 "<<z2<<std::endl;

   a=(-y1*z2+z2*y0+z1*y2-y2*z0+z0*y1-y0*z1)/(y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
   b=-(x2*z1-z1*x0-z2*x1+x1*z0-z0*x2+x0*z2)/(y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
   c=(y2*x1*z0-z0*x2*y1+z2*y1*x0-y2*z1*x0+y0*x2*z1-z2*x1*y0)/(y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);

   //std::cout<<"interpolated z is "<<a*x+b*y+c<<std::endl;
   //std::cout<<"at x = "<<x<<" y = "<<y<<std::endl;
   return(a*x+b*y+c);
   
}

double LineFit(double x1, double y1, double x2, double y2, double nx)
{
   const double slope = (y2-y1)/(x2-x1);
   return slope*(nx-x1)+y1;
}                  

double DistanceBW2Points(double x1, double y1, double x2, double y2)
{
   return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
