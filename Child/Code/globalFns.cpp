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
   return UnitVector( dx, dy );
//    const double mag = sqrt( dx * dx + dy * dy );
//    return
//      tArray< double >( dx/mag, dy/mag );
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
**	z(z1, z2, z3, z4, X, Y) = z1 * XY + z2 * (1-X)Y + z3 * X(1-Y) 
**				  + z4 * (1-X)(1-Y),  0<=X,Y<=1
**
**	 Assumes: passed coordinates are normalized to grid spacing = 1.
**   Created: 10/98 SL
**   Modified: to reflect difference between y (positive "up") and
**      j (pos. "down")
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

/**************************************************************************/
/**
**  InBoundsOnMaskedGrid
**
**   Adapted from InterpSquareGrid. Instead of interpolating gridded 
**   values, simply returns a boolean indicating whether coordinates
**   (xgen, ygen) are within the bounds, i.e., point is among grid points
**   with data, or out of bounds, i.e., among grid points with "nodata."
**
**	 Assumes: passed coordinates are normalized to grid spacing = 1
**     and lower-left corner = (0,0), and reflects difference between y
**     (positive "up") and j (pos. "down")
**
**   Created: 9/2010 SL
*/
/**************************************************************************/
bool InBoundsOnMaskedGrid( double xgen, double ygen, 
			   tMatrix<double> const &elev, int mask )
{
  int maskcount = 0;
  const int nx = elev.getNumCols();
  const int ny = elev.getNumRows();
   
  const int ix = static_cast<int>(xgen);  //                      z3----z4
  const int iy = static_cast<int>(ygen);  //                      |     |
  const double xrem = xgen - static_cast<double>(ix); //          |     |
  const double yrem = ygen - static_cast<double>(iy); //          z1----z2
  const double yrows = elev.getNumRows();
  const int i = ix;
  const int j = static_cast<int>(yrows - 1.0 - ygen);
  if( i < 0 || j < 0 || i > nx - 1 || j > ny - 1 ) return false;

  const double z3 = elev( j, i );
  const double z4 = ( i < nx - 1 ) ? elev( j, i + 1 ) : mask;
  const double z1 = ( j < ny - 1 ) ? elev( j + 1, i ) : mask;
  const double z2 = ( i < nx - 1 && j < ny - 1 ) ? elev( j + 1, i + 1 ) : mask;

  if( z1 == mask ) maskcount++;
  if( z2 == mask ) maskcount += 2;
  if( z3 == mask ) maskcount += 4;
  if( z4 == mask ) maskcount += 8;
  if( maskcount == 15 ) return false;                      // all corners are masked
  if( maskcount == 0 ) return true;                     // all corners are in bounds
  if( ( maskcount == 9 &&				    // z1 and z4 are masked
	( ( xrem <= 0.5 && yrem >= 0.5 ) || ( xrem >= 0.5 && yrem <= 0.5 ) ) )
      || ( maskcount == 6 &&				   // z2 and z3 are masked
	   ( ( xrem <= 0.5 && yrem <= 0.5 ) || ( xrem >= 0.5 && yrem >= 0.5 ) ) ) 
      || ( maskcount == 1 &&				   // z1 is masked
	   !( xrem < 0.5 && yrem < 0.5 ) ) 
      || ( maskcount == 2 &&				   // z2 is masked
	   !( xrem > 0.5 && yrem < 0.5 ) )
      || ( maskcount == 4 &&				   // z3 is masked
	   !( xrem < 0.5 && yrem > 0.5 ) )
      || ( maskcount == 8 &&				   // z4 is masked
	   !( xrem > 0.5 && yrem > 0.5 ) )
      || ( maskcount == 3 && yrem >= 0.5 )                 // z1 & z2 are masked
      || ( maskcount == 10 && xrem <= 0.5 )                // z2 & z4 are masked
      || ( maskcount == 12 && yrem <= 0.5 )                // z3 & z4 are masked
      || ( maskcount == 5 && xrem >= 0.5 )                 // z1 & z3 are masked
      || ( maskcount == 14 && xrem <= 0.5 && yrem <= 0.5 ) // z1 is in bounds
      || ( maskcount == 13 && xrem >= 0.5 && yrem <= 0.5 ) // z2 is in bounds
      || ( maskcount == 11 && xrem <= 0.5 && yrem >= 0.5 ) // z3 is in bounds
      || ( maskcount == 7 && xrem >= 0.5 && yrem >= 0.5 ) )// z4 is in bounds
    return true;
  else
    return false;
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

/***********************************************************************\
  PlaneFit
  Overloaded version that finds elevation as z.
  Fits a plane to three points and returns the "z" value on this plane
  at the tx and ty values passed to this function.
  nds = Array of size three pointers to nodes.  Use the x and y
        values from these three nodes for the plane fit.
  zs = Array of size three of the "z" values used for the fit.

  Modified 9/2010 SL
\***********************************************************************/
double PlaneFit(double tx, double ty, tTriangle* tri )
{
   double y0 = tri->pPtr(0)->getY();
   double y1 = tri->pPtr(1)->getY();
   double y2 = tri->pPtr(2)->getY();
   double x0 = tri->pPtr(0)->getX();
   double x1 = tri->pPtr(1)->getX();
   double x2 = tri->pPtr(2)->getX();
   double z0 = tri->pPtr(0)->getZ();
   double z1 = tri->pPtr(1)->getZ();
   double z2 = tri->pPtr(2)->getZ();

   double a = (-y1*z2+z2*y0+z1*y2-y2*z0+z0*y1-y0*z1)
       / (y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
   double b = -(x2*z1-z1*x0-z2*x1+x1*z0-z0*x2+x0*z2)
       / (y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
   double c = (y2*x1*z0-z0*x2*y1+z2*y1*x0-y2*z1*x0+y0*x2*z1-z2*x1*y0)
       / (y2*x1-x1*y0-x2*y1+y1*x0-x0*y2+y0*x2);
   return(a*tx+b*ty+c);  
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

/*****************************************************************************\
**
**      DistanceToLine: given x,y coords, finds distance to the line
**              defined by given a, b, and c (ax + by + c = 0)
**      Global function.
**
**      Data members updated:
**      Called by:
**      Calls:
**
\*****************************************************************************/
double DistanceToLine( double x2, double y2, double a, double b, double c )
{
   double x, y;
   double f = -b;
   double g = a;
   double h = b * x2 - a * y2;
   if( fabs(b) > 0 && fabs(a) > 0 )
   {
      x = (g * c / b - h) / (f - g * a / b);
      y = (-c - a * x) / b;
   }
   else
   {
      if( fabs(a) == 0.0 )
      {
         y = -c / b;
         x = -h / f;
      }
      else
      {
         y = -h / g;
         x = -c / a;
      }
   }
   double d = sqrt( (x2 - x) * (x2 - x) + (y2 - y) * (y2 - y) );
   return d;
}


/*****************************************************************************\
**
**      DistanceToLine: given x,y coords, finds distance to the line
**              formed by points  p0->(x, y) and p1->(x, y)
**      Global function.
**
**      Data members updated:
**      Called by:
**      Calls:
**
\*****************************************************************************/
double DistanceToLine( double x2, double y2, tNode const *p0, tNode const *p1 )
{
   double a, b, c, x0, y0, x1, y1, d;

   x0 = p0->getX();
   y0 = p0->getY();
   x1 = p1->getX();
   y1 = p1->getY();
   a = y1 - y0;
   b = x0 - x1;
   c = -( a * x0 + b * y0 );

   d = DistanceToLine( x2, y2, a, b, c );
   return d;
}


// global function to be used with veg. growth stuff:
  /*
   * this equation was called Richard's Chapman equation (in the Ecological
   * world ?). Here I used the form as introduced by Mark Harmon.
   * -Duan
   */
// added by SL, 8/10
double Richards_Chapman_equ( const double t, const double Smax, 
			     const double decay, const double shape )
{return Smax * pow( ( 1.0 - exp( -decay * t ) ), shape );}

//
// New computational geometry global functions by SL, 7/2010
//
// DotProduct2D: returns dot product of edge from p0 to p1 (presumed to exist;
// doesn't actually use edge) and (unit) vector passed by reference:
// -STL, 7/2010
double DotProduct2D( tLNode *p0, tLNode *p1, tArray<double> &uv ) 
{	    
  return 
    ( p1->getX() - p0->getX() ) * uv[0] 
    + ( p1->getY() - p0->getY() ) * uv[1];
}


/***************************************************************************\
 **  void ScanLineNodesByPublicFlag( tLNode*, vector<tLNode*> & )
 **
 **  Global function to find a node's neighbors that are closest to line.
 **  Scans as long as neighbors have the same value of tNode::public1.
 **  Makes an array of nodes up to and including the first nodes with
 **  different values of tNode::public1.
 **  Modified from (OSU CHILD) tDebrisFlow::FindBanks.
 **
 **  Called by: tDebrisFlow::Initialize.
 **  Takes: a pointer to a node and an (empty) array for pointers to the 
 **    scanned nodes.
 **  Calls: 
 **    - LinePerpendicularToFlowEdge( tLNode*, tArray<double>& ) (global)
 **    - FindNbrsOnLine( tLNode*, tArray<double>&, vector<tLNode*>& ) (global)
 **    - ScanLineNextNode( tLNode*, tLNode*, tArray<double>& ) (global)
 **  Changes: values in passed tLNode* array; gives the neighbors.
 **
 **  - STL, 7/2010
\***************************************************************************/
void ScanLineNodesByPublicFlag( tLNode* atPtr, 
				tPtrList<tLNode> &scanLineNodesList )
{
  const int numCluster = atPtr->public1;
  vector< tLNode*> BVec(2);
  tArray<double> abc(3);
  // find the line perpendicular to the flowedge of the starting node:
  // (set values of tArray<double> abc):
  LinePerpendicularToFlowEdge( atPtr, abc ); // (global function)
  // find the nodes to left and right closest to that line:
  // (set values of vector<tLNode*> BVec):
  FindNbrsOnLine( atPtr, abc, BVec ); // (global function)
  // find nodes to left up to and including the first node outside cluster;
  // leave nodes in left-to-right order:
  tPtrList<tLNode> leftList;
  leftList.insertAtFront( BVec[0] );
  if( BVec[0]->public1 == numCluster )
    {
      tLNode* prvPtr = atPtr;
      tLNode* curPtr = BVec[0];
      bool inCluster;
      do
	{
	  inCluster = false;
	  tLNode* nxtPtr = ScanLineNextNode( prvPtr, curPtr, abc ); // global fn
	  leftList.insertAtFront( nxtPtr );
	  if( nxtPtr->public1 == numCluster ) 
	    {	 
	      prvPtr = curPtr;
	      curPtr = nxtPtr;
	      inCluster = true;
	    }
	} while( inCluster );
    }
  // find nodes to right up to and including first node outside cluster:
  tPtrList<tLNode> rightList;
  rightList.insertAtBack( BVec[1] );
  if( BVec[1]->public1 == numCluster )
    {
      tLNode* prvPtr = atPtr;
      tLNode* curPtr = BVec[1];
      bool inCluster;
      do
	{
	  inCluster = false;
	  tLNode* nxtPtr = ScanLineNextNode( prvPtr, curPtr, abc ); // global fn
	  rightList.insertAtBack( nxtPtr );
	  if( nxtPtr->public1 == numCluster ) 
	    {	 
	      prvPtr = curPtr;
	      curPtr = nxtPtr;
	      inCluster = true;
	    }
	} while( inCluster);
    }
//   const int numScanLine = leftList.getSize() + rightList.getSize() + 1;
  // take nodes from lists and center and put into array from left to right:
//   scanLineNode.resize( numScanLine );
  {
    //     int i=0;
    while( tLNode* node = leftList.removeFromFront() )
      scanLineNodesList.insertAtBack( node );
    scanLineNodesList.insertAtBack( atPtr );
    while( tLNode* node = rightList.removeFromFront() )
      scanLineNodesList.insertAtBack( node );
    //     while( ( scanLineNode[i] = leftList.removeFromFront() ) ) ++i;
    //     scanLineNode[i] = atPtr;
    //     ++i;
    //     while( ( scanLineNode[i] = rightList.removeFromFront() ) ) ++i;
    //     assert( i == numScanLine );
  }
}

/***************************************************************************\
 **  void ScanLineNodesByDistance( tLNode* cn, double maxDist, 
 **                                vector<tLNode*> &scanLineNode )
 **
 **  Global function to find a node's neighbors that are closest to line.
 **  Scans as long as neighbors are within the specified distance.
 **  Makes an array of nodes up to and including the first nodes beyond
 **  the prescibed distance, but not including boundary nodes.
 **  Modified from (OSU CHILD) tDebrisFlow::FindBanks.
 **
 **  Called by: tDebrisFlow::RunScourDeposit.
 **  Takes: a pointer to a node, a maximum distance, and an (empty) array 
 **    for pointers to the scanned nodes.
 **  Calls: 
 **    - LinePerpendicularToFlowEdge( tLNode*, tArray<double>& ) (global)
 **    - FindNbrsOnLine( tLNode*, tArray<double>&, vector<tLNode*>& ) (global)
 **    - LeftUnitVector( tEdge* ) (global)
 **    - RightUnitVector( tEdge* ) (global)
 **    - DotProduct2D( tLNode*, vector<tLNode*>&, tArray<double>& ) (global)
 **    - ScanLineNextNode( tLNode*, tLNode*, tArray<double>& ) (global)
 **  Changes: values in passed tLNode* array
 **  Returns: index to original node (cn) in array of nodes (scanLineNode).
 **
 **  - STL, 7/2010
\***************************************************************************/
int ScanLineNodesByDistance( tLNode *cn, double maxDist, 
			     tPtrList<tLNode> &scanLineNodesList )
{
  vector<tLNode*> BVec(2);
  tArray<double> abc(3);
  // find the line perpendicular to the flowedge of the downstream-most node
  // (set values of tArray<double> abc):
  LinePerpendicularToFlowEdge( cn, abc ); // (global function)
  // find the nodes to left and right closest to that line
  // (set values of vector<tLNode*> BVec):
  FindNbrsOnLine( cn, abc, BVec ); // (global function)
  // find nodes to left up to and including the first node beyond the
  // prescribed distance;
  // leave nodes in left-to-right order:
  tPtrList<tLNode> leftList;
  leftList.insertAtFront( BVec[0] );
  tArray<double> uVec = LeftUnitVector( cn->getFlowEdg() );
  double perpDist = fabs( DotProduct2D( cn, BVec[0], uVec ) );
  if( perpDist <= maxDist && BVec[0]->isNonBoundary() )
    {
      tLNode* prvPtr = cn;
      tLNode* curPtr = BVec[0];
      bool closeEnough;
      do
	{
	  closeEnough = false;
	  tLNode* nxtPtr = ScanLineNextNode( prvPtr, curPtr, abc ); // global fn
	  if( nxtPtr > 0 && nxtPtr->isNonBoundary() )
	    {
	      leftList.insertAtFront( nxtPtr );
	      perpDist += fabs( DotProduct2D( curPtr, nxtPtr, uVec ) );
	      if( perpDist <= maxDist ) 
		{	 
		  prvPtr = curPtr;
		  curPtr = nxtPtr;
		  closeEnough = true;
		}
	    }
	} while( closeEnough );
    }
  // find nodes to right up to and including first node beyond the 
  // prescribed distance:
  tPtrList<tLNode> rightList;
  rightList.insertAtBack( BVec[1] );
  uVec = RightUnitVector( cn->getFlowEdg() );
  perpDist = fabs( DotProduct2D( cn, BVec[0], uVec ) );
  if( perpDist <= maxDist && BVec[1]->isNonBoundary() )
    {
      tLNode* prvPtr = cn;
      tLNode* curPtr = BVec[1];
      bool closeEnough = false;
      do
	{
	  closeEnough = false;
	  tLNode* nxtPtr = ScanLineNextNode( prvPtr, curPtr, abc ); // global fn
	  if( nxtPtr > 0 && nxtPtr->isNonBoundary() )
	    {
	      rightList.insertAtBack( nxtPtr );
	      perpDist += fabs( DotProduct2D( curPtr, nxtPtr, uVec ) );
	      if( perpDist <= maxDist ) 
		{	 
		  prvPtr = curPtr;
		  curPtr = nxtPtr;
		  closeEnough = true;
		}
	    }
	} while( closeEnough);
    }
  const int numScanLine = leftList.getSize() + rightList.getSize() + 1;
//   scanLineNode.resize( numScanLine );
  int iStream;
  { // take nodes from lists and center and put into array from left to right:
    int i=0;
    while( tLNode* node = leftList.removeFromFront() )
      {
	scanLineNodesList.insertAtBack( node );
	++i;
      }
    // record index of "in stream" node:
    iStream = i;
    scanLineNodesList.insertAtBack( cn );
    ++i;
    while( tLNode* node = rightList.removeFromFront() )
      {
	scanLineNodesList.insertAtBack( node );
	++i;
      }
    //     while( ( scanLineNode[i] = leftList.removeFromFront() ) ) ++i;
    //     scanLineNode[i] = cn;
    //     // record index of "in stream" node:
    //     iStream = i;
    //     ++i;
    //     while( ( scanLineNode[i] = rightList.removeFromFront() ) ) ++i;
    assert( i == numScanLine );
  }
  return iStream;
}

/***************************************************************************\
 **  void FindNbrsOnLine( tLNode*, tArray<double> &, vector<tLNode*> & )
 **
 **  Global function to find a node's neighbors that are closest to line.
 **  Modified from (OSU CHILD) tDebrisFlow::FindBanks.
 **  Called by: 
 **    - ScanLineNodesByPublicFlag (global) and 
 **    - ScanLineNodesByDistance (global).
 **  Takes: a pointer to a node, reference to an array with the line eq,
 **    and an (empty) array for pointers to the neighbor nodes.
 **  Calls: 
 **    - DistanceToLine( x, y, a, b, c ) (global)
 **  Changes: values in passed tLNode* array; gives the neighbors.
 **
 **  - STL, 7/2010
\***************************************************************************/
void FindNbrsOnLine( tLNode* atPtr, tArray<double> &abc, 
		     vector<tLNode*> &BVec )
{
   tSpkIter spI( atPtr );
   const int numspokes = spI.getNumSpokes();
   const int arrsize = numspokes * 2;
   tArray< double > spD( arrsize ),
       spR( arrsize );

   tEdge* fe = atPtr->getFlowEdg();
   tEdge* ce = fe;
   int i=0;
   // find distance and remainders of pts wrt line perpendicular 
   // to downstream direction
   do
     {
       tArray2<double> xy;
       ce->getDestinationPtr()->get2DCoords( xy );
       spD[i] = spD[numspokes + i] = abc[0] * xy.at(0) + abc[1] * xy.at(1) + abc[2];
       spR[i] = spR[numspokes + i] = 
	 DistanceToLine( xy.at(0), xy.at(1), abc[0], abc[1], abc[2] );
       ce = ce->getCCWEdg();
       ++i;
     } while( ce != fe );
   // find point pairs:
   i=0;
   int j=0;
   ce = fe;
   do
     {
       // find signs of 'this' and 'next' remainders
       double s1 = 0.0;
       if( spD[i] != 0.0 ) s1 = spD[i] / fabs( spD[i] );
       double s2 = 0.0;
       if( spD[i + 1] != 0.0 ) s2 = spD[i + 1] / fabs( spD[i + 1] );
       if( s1 != s2 ) // points are on opposite sides of the perp.
	 {
	   assert( j<2 );
	   // find shorter distance and use that node for bank:
	   if( spR[i] < spR[i + 1] )
	     BVec[j] = static_cast<tLNode*>( ce->getDestinationPtrNC() );
	   else
	     BVec[j] = static_cast<tLNode*>( ce->getCCWEdg()->getDestinationPtrNC() );
	   // rare case: one of the remainders is zero:
	   // need to make sure we don't choose the same node for both banks
	   if( s1 == 0.0 || s2 == 0.0  )
             ++i; // advance an extra step
	   ++j;
	 }
       ++i;
       ce = ce->getCCWEdg();
     } while( ce != fe );
}

/***************************************************************************\
 **  void LinePerpendicularToFlowEdge( tLNode* atPtr, tArray<double> &abc )
 **
 **  Global function to find line perpendicular to node's flowedge.
 **  Called by:
 **    - ScanLineNodesByPublicFlag (global) and 
 **    - ScanLineNodesByDistance (global).
 **  Takes: a pointer to a node, reference to an array (preferably empty).
 **  Calls: 
 **    - UnitVector( tEdge* ).
 **  Changes: values in passed array; gives the equation of the line.
 **
 **  - STL, 7/2010
\***************************************************************************/
void LinePerpendicularToFlowEdge( tLNode* atPtr, tArray<double> &abc )
{
  tArray<double> dxy = UnitVector( atPtr->getFlowEdg() );
  const double x = atPtr->getX();
  const double y = atPtr->getY();
  LinePerpendicularToVector( x, y, dxy, abc );
}

/***************************************************************************\
 **  tArray<double> RightUnitVector( tEdge* )
 **
 **  Global function to return an array with the magnitudes of the x and y
 **  components of the unit vector perpendicular to the pointed-to edge and 
 **  to its right.
 **  Called by:
 **    - tDebrisFlow::Initialize 
 **    - ScanLineNodesByDistance (global)
 **    - FillCrossSection (global)
 **  Takes: a pointer to an edge.
 **  Calls: 
 **  Changes: 
 **  Returns: array with components of unit vector
 **
 **  - STL, 2/99 (incorporated into CU CHILD, 9/2010)
\***************************************************************************/
tArray< double > RightUnitVector( tEdge* ePtr )
{
  //assert( ePtr != 0 );
   tNode* on = ePtr->getOriginPtrNC();
   tNode* dn = ePtr->getDestinationPtrNC();
   double dx = dn->getY() - on->getY();
   double dy = on->getX() - dn->getX();
   return UnitVector( dx, dy );
}

/***************************************************************************\
 **  tArray<double> RightUnitVector( tEdge* )
 **
 **  Global function to return an array with the magnitudes of the x and y
 **  components of the unit vector perpendicular to the pointed-to edge and 
 **  to its left.
 **  Called by:
 **    - tDebrisFlow::Initialize 
 **    - ScanLineNodesByDistance (global).
 **  Takes: a pointer to an edge.
 **  Calls: 
 **  Changes: 
 **  Returns: array with components of unit vector
 **
 **  - STL, 2/99 (incorporated into CU CHILD, 9/2010)
\***************************************************************************/
tArray< double > LeftUnitVector( tEdge* ePtr )
{
  //assert( ePtr != 0 );
   tNode* on = ePtr->getOriginPtrNC();
   tNode* dn = ePtr->getDestinationPtrNC();
   double dx = on->getY() - dn->getY();
   double dy = dn->getX() - on->getX();
   return UnitVector( dx, dy );
}

/******************************************************************************\
 **  ScanLineNextNode: This is modified from tStreamMeander::FindBankErody.
 **    Check all nbrs; find distances to line. Find pair of
 **    consecutive nbrs which fall on either side of line (going ccw from
 **    edge to previous point on line), and then find which point of the 
 **    pair is closest to the line. 
 **
 **  Takes:
 **    - prvPtr: pointer to previous node on scanline
 **    - atPtr: pointer to current node on scanline; this and prvPtr define
 **      scanning direction
 **    - abc: reference to three member array containing, a, b, and c, 
 **      such that for a point (x,y) on the line, a*x + b*y + c = 0
 **  Calls: global function DistanceToLine(x, y, a, b, c)
 **  Changes: nothing
 **  Returns: pointer to next node (tLNode*) on scanline
 **
 **  Created: 7/2010 SL
\*****************************************************************************/
tLNode* ScanLineNextNode( tLNode* prvPtr, tLNode* atPtr, tArray<double> &abc )
{
   tSpkIter spI( atPtr );
   const int numspokes = spI.getNumSpokes();
   const int arrsize = numspokes * 2;
   tArray< double > spD( arrsize ),
       spR( arrsize );
   tEdge* pe = atPtr->EdgToNod( prvPtr );
   tEdge* ce = pe;
   int i=0;
   // find distance and remainders of pts wrt line perpendicular 
   // to downstream direction
   do
     {
       tArray2<double> xy;
       ce->getDestinationPtr()->get2DCoords( xy );
       spD[i] = spD[numspokes + i] = abc[0] * xy.at(0) + abc[1] * xy.at(1) + abc[2];
       spR[i] = spR[numspokes + i] = 
	 DistanceToLine( xy.at(0), xy.at(1), abc[0], abc[1], abc[2] );
       ce = ce->getCCWEdg();
       ++i;
     } while( ce != pe );
   // find point pairs:
   i=0;
//    tLNode *nxtPtr=0;
   // start at edge to previous node:
   ce = pe;
   do
     {
       // find signs of 'this' and 'next' remainders
       double s1 = 0.0;
       if( spD[i] != 0.0 ) s1 = spD[i] / fabs( spD[i] );
       double s2 = 0.0;
       if( spD[i + 1] != 0.0 ) s2 = spD[i + 1] / fabs( spD[i + 1] );
       if( s1 != s2 ) // points are on opposite sides of the perp.
	 {
	   if( i > 0 ) // not still at edge to previous node
	     {
	       // find shorter distance and return that node:
	       if( spR[i] < spR[i + 1] )
		 return static_cast<tLNode*>( ce->getDestinationPtrNC() );
// 		 nxtPtr = static_cast<tLNode*>( ce->getDestinationPtrNC() );
	       else
		 return static_cast<tLNode*>( ce->getCCWEdg()->getDestinationPtrNC() );
// 		 nxtPtr = static_cast<tLNode*>( ce->getCCWEdg()->getDestinationPtrNC() );
	       // rare case: one of the remainders is zero:
	       // need to make sure we don't choose the same node for both banks
	       if( s1 == 0.0 || s2 == 0.0  )
		 ++i; // advance an extra step
	     }
	 }
       ++i;
       ce = ce->getCCWEdg();
     } while( ce != pe );
   return 0;
}

/***************************************************************************\
 **  void FillCrossSection( int iStream, double areaFrontal, 
 **                         double widthFrontal,
 **		            vector<tLNode*> &scanLineNode, 
 **		            vector<double> &leftXYZ, 
 **                         vector<double> &rightXYZ )
 **
 **  Global function to find vertices of defining area filled by a given
 **  cross-sectional area. Starting with stream node, finds bank heights,
 **  and finds whether given cross-sectional area will cover them. 
 **  Finally, finds the points halfway along edges between last filled nodes
 **  and next nodes on cross-section.
 **
 **  Called by: tDebrisFlow::RunScourDeposit.
 **  Takes: 
 **    - iStream: index to stream node;
 **    - areaFrontal: area to fill cross-section;
 **    - widthFrontal: maximum spreading width;
 **    - scanLineNode: nodes along scan line;
 **    - leftXYZ, rightXYZ: 3D coords of left and right endpoints, 
 **        respectively, of filled cross-section.
 **  Calls: 
 **    - RightUnitVector( tEdge* ), LeftUnitVector (tEdge* ).
 **  Changes: values in passed arrays leftXYZ and rightXYZ.
 **
 **  - STL, 7/2010
\***************************************************************************/
void FillCrossSection( int iStream, double areaFrontal, double widthFrontal,
		       vector<tLNode*> &scanLineNode, 
		       tArray<double> &leftXYZ, tArray<double> &rightXYZ )
{ // find depths by spreading flow to successive neighbor nodes
  // until no longer over-topping banks
  // find perpendicular distances from left to right with dot product:
  const int numScanLine = scanLineNode.size();
  tArray<double> perpDistance( numScanLine - 1 );
  tArray<double> rightUnitVec = 
    RightUnitVector( scanLineNode[iStream]->getFlowEdg() );
  for( int i=0; i<numScanLine-1; ++i )
    // dot product of edge and right unit vector:
    perpDistance[i] = 
      fabs( DotProduct2D( scanLineNode[i], scanLineNode[i+1], rightUnitVec ) );
  // fill cross-section:
  int iL = iStream-1;
  int iR = iStream+1;
  int iNew = iStream;
  double zStream = scanLineNode[iStream]->getZ();
  double leftPerpDist = perpDistance[iNew-1] / 2.0;
  double rightPerpDist = perpDistance[iNew] / 2.0;
  double newWidth = leftPerpDist + rightPerpDist;
  double newHeight = areaFrontal / newWidth;
  double maxBankHeight = 0.0;
  const double kLargeBankHeight = newHeight * 100.;
  double lBankHeight = scanLineNode[iL]->getZ() - zStream;
  if( iL == 0 ) lBankHeight = kLargeBankHeight;
  double rBankHeight = scanLineNode[iR]->getZ() - zStream;
  if( iR == numScanLine-1 ) rBankHeight = kLargeBankHeight;
  while( newWidth < widthFrontal 
	 && ( newHeight > lBankHeight
	      || newHeight > rBankHeight ) )
    {
      double bankHeight = 0.0;
      double addWidth = 0.0;
      if( lBankHeight < rBankHeight )
	{
	  bankHeight = lBankHeight;
	  iNew = iL;
	  addWidth = ( perpDistance[iNew-1] + perpDistance[iNew] ) / 2.0;
	  leftPerpDist += addWidth;
	  --iL;
	  if( iL > 0 ) // don't consider last point
	    lBankHeight = scanLineNode[iL]->getZ() - zStream;
	  else
	    lBankHeight = kLargeBankHeight;
	}
      else
	{
	  bankHeight = rBankHeight;
	  iNew = iR;
	  addWidth = ( perpDistance[iNew-1] + perpDistance[iNew] ) / 2.0;
	  rightPerpDist += addWidth;
	  ++iR;
	  if( iR < numScanLine-1 )
	    rBankHeight = scanLineNode[iR]->getZ() - zStream;
	  else
	    rBankHeight = kLargeBankHeight;
	}
      double oldWidth = newWidth;
      newWidth += addWidth;
      double oldHeight = newHeight;
      if( bankHeight >= maxBankHeight ) 
	{
	  maxBankHeight = bankHeight;
	  newHeight =
	    ( oldHeight * oldWidth + maxBankHeight * addWidth ) / newWidth;
	}
      else
	{
	  // subtract pit depth from amount to spread
	  double pitFillDiff = 
	    ( oldHeight - maxBankHeight ) * oldWidth / addWidth 
	    -  maxBankHeight + bankHeight;
	  if( pitFillDiff > 0.0 )
	    newHeight = maxBankHeight + pitFillDiff / newWidth;
	  else
	    { // not enough to overtop pit
	      newHeight = maxBankHeight;
	      newWidth = widthFrontal + 1.0; // MAKE IT NOT SPREAD ANYMORE
	    }
	}
    }
  // coordinates for scour zone are stream averages of last inundated nodes and 
  // nodes next on scanline for x and y; for z use top of flow:
  tArray<double> leftUnitVec = 
    LeftUnitVector( scanLineNode[iStream]->getFlowEdg() );
  leftXYZ[0] = scanLineNode[iStream]->getX() + leftPerpDist * leftUnitVec[0];
  leftXYZ[1] = scanLineNode[iStream]->getY() + leftPerpDist * leftUnitVec[1];
  leftXYZ[2] = zStream + newHeight;
  rightXYZ[0] = scanLineNode[iStream]->getX() + rightPerpDist * rightUnitVec[0];
  rightXYZ[1] = scanLineNode[iStream]->getY() + rightPerpDist * rightUnitVec[1];
  rightXYZ[2] = zStream + newHeight;
}

