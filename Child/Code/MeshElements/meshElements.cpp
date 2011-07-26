/**************************************************************************/
/**
**  @file meshElements.cpp
**  @brief Functions for mesh element classes tNode, tEdge, and
**         tTriangle. (formerly called gridElements.cpp)
**
**  For a description of these 3 classes, see meshElements.h.
**
**  Modifications:
**   - gt added code to recreate the "edg" implementation as an alternative
**     to spokelist. Dec 1 '97.
**   - previously separate tNode, tEdge, and tTriangle files grouped into
**     "gridElements", 1/20/98 gt
**   - added tNode::AttachNewSpoke and tEdge::WelcomeCCWNeighbor gt 2/99
**   - 2/2/00: GT transferred get/set, constructors, and other small
**     functions to .h file to inline them
**   - 2/2000 GT added tNode functions getVoronoiVertexList and
**     getVoronoiVertexXYZList to support dynamic remeshing.
**
**  $Id: meshElements.cpp,v 1.76 2004-06-16 13:37:27 childcvs Exp $
*/
/**************************************************************************/

#include <assert.h>
#include "meshElements.h"
#include "../globalFns.h" // For PlaneFit; this could go in geometry; TODO

/**  GLOBAL FUNCTIONS  ****************************************************/

/*************************************************************************\
**
**  FindIntersectionCoords
**
**  Finds and returns intersection of line segments
**  defined by endpoints given as arguments; 1st seg. endpts are xy1, xy2;
**  2nd seg. endpts are xy3, xy4. (SL)  TODO: add inline doc
**  (Note: does this need to be global? can it be tNode mbr?)
**
**  Called by:  tNode::ComputeVoronoiArea
**
\*************************************************************************/
tArray2< double > FindIntersectionCoords( tArray2< double > const &xy1,
					  tArray2< double > const &xy2,
					  tArray2< double > const &xy3,
					  tArray2< double > const &xy4 )
{
   tArray2< double > intxy;

   const double dxa = xy2.at(0) - xy1.at(0);
   const double dxb = xy4.at(0) - xy3.at(0);
   const double dya = xy2.at(1) - xy1.at(1);
   const double dyb = xy4.at(1) - xy3.at(1);
   const double &a = dya;
   const double b = -dxa;
   const double c = dxa * xy1.at(1) - dya * xy1.at(0);
   const double &f = dyb;
   const double g = -dxb;
   //const double h = dxb * xy3.at(1) - dyb * xy4.at(0);
   // seems to be a bug above; fixed here?
   const double h = dxb * xy3.at(1) - dyb * xy3.at(0);
   if( fabs(dxa) > 0 && fabs(dxb) > 0 )
   {
      if( fabs(f - g * a / b) > 0 )
      {
         intxy.at(0) = (g * c / b - h) / (f - g * a / b);
         intxy.at(1) = (-c - a * intxy.at(0)) / b;
      }
   }
   else
   {
      if( fabs(dya) > 0 && fabs(dyb) > 0 )
      {
         if( fabs(g - f * b / a) > 0 )
         {
            intxy.at(1) = (f * c / a - h) / (g - f * b / a);
            intxy.at(0) = (-c - b * intxy.at(1)) / a;
         }
      }
      else //one horiz. line and one vert. line:
      {
         if( fabs(dya) == 0 )
         {
            intxy.at(0) = xy3.at(0);
            intxy.at(1) = xy1.at(1);
         }
         else
         {
            intxy.at(0) = xy1.at(0);
            intxy.at(1) = xy3.at(1);
         }
      }
   }
   return intxy;
}

// set static boolean for tNode:
bool tNode::freezeElevations = false;

/*****************************************************************************\
**
**  tNode::Dist
**
**  Computes the distance of the node from the line formed
**  by points p0 p1 using x y. (replaces dis)
**
**      Data members updated: (none)
**      Called by:
**      Calls: (none)
**
\*****************************************************************************/
double tNode::Dist( tNode const * p0, tNode const * p1 ) const
{
  const double a=(p1->y)-(p0->y);
  const double b=-((p1->x)-(p0->x));
  const double c=-((a*(p0->x))+(b*(p0->y)));
  return fabs((a*x + b*y + c) / sqrt(a*a + b*b));
}


/*****************************************************************************\
**
**  tNode::EdgToNod
**
**  Finds and returns the spoke (edge) that connects the current node to _nod_,
**  or zero if no such spoke is found.
**
\*****************************************************************************/
tEdge *tNode::EdgToNod( tNode const * nod )
{
   tSpkIter spokIter( this );
   tEdge * ce;

   for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
   {
      if( ce->getDestinationPtr() == nod ) return ce;
   }
   return 0;
}


/*****************************************************************************\
**
**  tNode::ComputeVoronoiArea
**
**  Computes the node's Voronoi area by summing the area of embedded
**  triangles, and also calls CalcSpokeVEdgLengths to compute the length of
**  the sides of the Voronoi cell (the length of cell sides is needed if
**  for example transport between adjacent cells depends on the width of
**  their shared face).
**
**  The basic Voronoi polygon is described by the set of "right-hand
**  Voronoi vertices" associated with each spoke (edge). These vertices
**  are computed by setVoronoiVertices() as the intersection of the
**  perpendicular bisectors of consecutive spokes. However, in some cases
**  the basic polygon will be distorted, with consecutive vertices NOT
**  being counter-clockwise from one another. (This seems to be the result
**  of numerical errors (possibly in estimating the intersection of two
**  nearly parallel bisectors); according to Sugihara and Iri
**  [J. Comp. Geom. & Appl., 1994, v. 4, p. 179], each Delaunay triangle
**  should be associated with one Voronoi vertex). This is handled by
**  detecting "loops" in the Voronoi polygon and cutting them off by
**  taking the area of the closest (counterclockwise) intersection of
**  perpendicular bisectors.
**
**  TODO: for speed, return to triangle-based method of computing nodes.
**        (could save space too by only storing VVert in tri's not edges)
**
**  5/7/98, SL: Found a bug in this procedure that requires a not-so-simple
**  fix. When a migrating node approaches the boundary, part of its calculated
**  Voronoi area may lie outside the mesh domain. The area can even blow up
**  if the migrating node comes very close to the boundary.
**
**  If spokes connect to boundary nodes, find which side of the edge connecting
**  those boundary nodes the node in question falls on, and make sure all
**  Voronoi vertices fall on the same side. If not, find intersections of
**  Voronoi edges with the boundary edge and put those intersections in vertex
**  list when calculating Vor. areas. To implement this, put in a list of
**  vertex coordinates from which the area will finally be calculated.
**
\*****************************************************************************/
double tNode::ComputeVoronoiArea()
{
  double area = 0;
  double dx, dy, dx0, dx1, dy0, dy1, dx2, dy2;
  double vx, vy, x0, y0, x1, y1, x2, y2, m1, m2;
  tEdge * ce, *edgptr;
  int i;

  // the following is a temporary hack, which should be replaced by a more
  // proper handling but i'm going away for a week tomorrow and don't
  // have time to deal w/ it right now (gt, aug 02)
  if( getBoundaryFlag()!=kNonBoundary ) // if not an interior node, something's wrong
    {
      std::cout << "Warning: attempt to compute Voronoi area for a boundary node: "
	   << id << " " << x << " " << y
	   << " " << BoundName(getBoundaryFlag()) << std::endl;
      return 0.0;
    }

  // Create a duplicate list of edges; we will modify this list to obtain
  // the correct vertices. In some cases, we may need to delete an edge
  // to get the correct list of vertices; we don't want to delete the
  // spoke ptr, so we make a duplicate list.
  //if( id==83 ) std::cout << "NODE 83: " << x << "," << y << std::endl;

  tPtrList< tEdge > vedgList;
  tPtrListIter< tEdge > vtxIter( vedgList );
  ce = getEdg();
  do
    {
      assert( ce!=0 );
      vedgList.insertAtBack( ce );
      if (0) {//DEBUG
	tArray2< double > const &xy1 = ce->getRVtx();
	tArray2< double > const &xy2 =
	  vedgList.getLast()->getPtr()->getRVtx();
	std::cout
	  << xy1.at(0) << " " << xy1.at(1) << "; "
	  << xy2.at(0) << " " << xy2.at(1) << std::endl;
	//if( id==83) std::cout << " " << ce->getDestinationPtr()->getX() << ","
	//                 << ce->getDestinationPtr()->getY() << std::endl;
      }
      ce = ce->getCCWEdg();
    } while( ce != getEdg() );
  vedgList.makeCircular();
  //std::cout << std::endl;
  // Check boundary status: Voronoi area only defined for non-boundary nodes
  if( getBoundaryFlag() == kNonBoundary )
    {
      {
	tArray2< double > xy1, xy2;
	bool cw = true;
	//std::cout << "find clockwise loops" << std::endl;
	do
	  {
	    // go through the list; we want the vertex list to run CCW;
	    // in some cases of long skinny triangles, the 'unimproved'
	    // v. polygon sides may form loops; loops are detected by
	    // finding two (2) consecutive 'CCW' vertices; i.e., where
	    // the 'curvature' is CW rather than CCW. In such cases,
	    // we delete one of the edges from the vertex list and find
	    // the new vertex at the intersection of the perp. bisectors
	    // of the edges to either 'side' of the deleted edge. Iterate.
	    // Really. It works.
	    cw = false;
	    for( ce=vtxIter.FirstP(); !( vtxIter.AtEnd() ); ce=vtxIter.NextP() )
	      {
		tArray2< double > const &xy = ce->getRVtx();
		tArray2< double > const &xyn = vtxIter.NextP()->getRVtx();
		tArray2< double > const &xynn = vtxIter.NextP()->getRVtx();
		dx0 = xynn.at(0) - xyn.at(0);
		dy0 = xynn.at(1) - xyn.at(1);
		dx1 = xy.at(0) - xyn.at(0);
		dy1 = xy.at(1) - xyn.at(1);
		if( dy0 * dx1 > dx0 * dy1 ) //clockwise
		  {
		    tArray2< double > const &xynnn = vtxIter.NextP()->getRVtx();
		    dx0 = xynnn.at(0) - xynn.at(0);
		    dy0 = xynnn.at(1) - xynn.at(1);
		    dx1 = xyn.at(0) - xynn.at(0);
		    dy1 = xyn.at(1) - xynn.at(1);
		    if( dy0 * dx1 > dx0 * dy1 ) //clockwise
		      {
			//two consecutive clockwise vertices=>want intersection
			//of bisectors of ce->nextedg and
			//ce->nextedg->nextedg->nextedg:
			cw = true;
			x0 = x; //node.x
			y0 = y; //node.y
			ce->getDestinationPtr()->get2DCoords( xy1 );
			vtxIter.PrevP()->getDestinationPtr()->get2DCoords( xy2 );
			x1 = ( x0 + xy1.at(0) ) / 2;
			y1 = ( y0 + xy1.at(1) ) / 2;
			x2 = ( x0 + xy2.at(0) ) / 2;
			y2 = ( y0 + xy2.at(1) ) / 2;
			dx1 = x1 - x0;
			dy1 = y1 - y0;
			dx2 = x2 - x0;
			dy2 = y2 - y0;
			if( fabs(dy1)>0 && fabs(dy2) > 0 )
			  {
			    m1 = -dx1/dy1;
			    m2 = -dx2/dy2;
			    vx = (y2-m2*x2-y1+m1*x1) / (m1-m2);
			    vy = m1*(vx-x1)+y1;
			  }
			else
			  {
			    if( fabs(dx1) > 0 && fabs(dx2) > 0 )
			      {
				m1 = dy1/dx1;
				m2 = dy2/dx2;
				vy=(m1*y1+x1-m2*y2-x2)/(m1-m2);
				vx= -vy*m1+m1*y1+x1;
			      }
			    //otherwise one vert., one horiz. line:
			    else if( fabs(dx1) > 0 )
			      {
				vx = x2;
				vy = y1;
			      }
			    else
			      {
				vx = x1;
				vy = y2;
			      }
			  }
			edgptr = vtxIter.PrevP();
			dx = xy.at(0) - vx;
			dy = xy.at(1) - vy;
			//std::cout << "reset vedglen and rvtx for edge "
			//     << edgptr->getID() << " to len "
			//     << sqrt( dx*dx + dy*dy )
			//     << ", x, y, " << xyn.at(0) << ", " << xyn.at(1) << std::endl;
			//reset 'next' edge's vertex to newly found intersection,
			//length adjusted accordingly
			edgptr->setVEdgLen( sqrt( dx*dx + dy*dy ) );
			edgptr->setRVtx( tArray2<double>( vx, vy ) );
			edgptr = vtxIter.ReportNextP();
			//std::cout << "reset vedglen and rvtx for edge "
			//     << edgptr->getID()
			//     << " to len 0.0, x, y, " << xynnn[0] << ", "
			//     << xynnn[1] << std::endl;
			//reset 'next-next' edge's vertex to the coordinates
			//of the 'next-next-next' edge's vertex; length to zero
			edgptr->setVEdgLen(0.0);
			edgptr->setRVtx( xynnn );
			//delete the offending vertex's edge from list
			/* edgptr =*/ vedgList.removeNext( vtxIter.NodePtr() );
		      }
		  }
		vtxIter.Get( ce->getID() );
	      }
	  } while( cw ); //while we're still finding loops in the polygon
      }

      //Before the next step, make a list of V. vertex coord. arrays.
      //In doing so, check for parts of the V. area lying outside the
      //mesh domain and cut them off by substituting coordinates of
      //two intersections of V. edges with boundary edge for the V.
      //vertex lying outside the boundary. This should take care of any
      //outlying area as long as all boundaries are convex.
      // Go through spokes and put RVtx of ccw edge in coord list, but
      // first check that the vtx lies within the bndies
      {
	tArray2< double > xy2, xy3;
	tList< tArray2< double > > vcL; // list of vertex coordinates
	tEdge *ne, *nne;
	for( ce = vtxIter.FirstP(); !(vtxIter.AtEnd()); ce = vtxIter.NextP() )
	  {
	    ne = ce->getCCWEdg();
	    tArray2<double> const &xy1 = ne->getRVtx();
	    //checking polygon edge is on boundary and ccw edge's RVtx is on
	    //wrong side of bndy edge...
	    if( ce->getBoundaryFlag() != kNonBoundary &&
		ne->getBoundaryFlag() != kNonBoundary)
	      {
		//if( id==83 ) std::cout << " CASE A\n";
		tNode *bn0, *bn1;
		bn0 = ce->getDestinationPtrNC();
		bn1 = ne->getDestinationPtrNC();
		bn0->get2DCoords( xy2 );
		bn1->get2DCoords( xy3 );
		if( !PointsCCW( xy1, xy2, xy3 ) )
		  {
		    //"cut off" portion of V. area outside bndy by finding intersections
		    //of V. edges and bndy edge:
		    //if( id==83 ) std::cout << " CASE B\n";
		    vcL.insertAtBack(
				     FindIntersectionCoords( ce->getRVtx(),
							     xy1, xy2, xy3 )
				     );
		    nne = ne->getCCWEdg();
		    vcL.insertAtBack(
				     FindIntersectionCoords( xy1,
							     nne->getRVtx(), xy2, xy3 )
				     );
		  }
		else vcL.insertAtBack( xy1 );
	      }
	    else vcL.insertAtBack( xy1 );
	  }

	// Now that we've found the correct vertices, make triangles to
	// fill the polygon; the sum of the tri areas is the v. area.
	// For a convex polygon, we can compute the total area as the
	// sum of the area of triangles [P(1) P(i) P(i+1)] for i=2,3...N-1.
	//std::cout << "find polygon area" << std::endl;
	// coords of first vertex:
	tListIter< tArray2< double > > vcI( vcL ); // iterator for coord list
	tArray2< double > const * const xy = vcI.FirstP(); //ce = vtxIter.FirstP();
	//xy = ce->getRVtx();
	// Find out # of vertices in polygon:
	const int nverts = vcL.getSize(); //vedgList.getSize();
	for( i=2; i<=nverts-1; i++ )
	  {
	    double a, b, c;

	    tArray2<double> const * const xyn = vcI.NextP(); //xyn = vtxIter.NextP()->getRVtx();// Vertex i
	    tArray2<double> const * const xynn = vcI.NextP();//vtxIter.ReportNextP()->getRVtx(); // Vertex i+1
	    {
	      const double dx = (*xyn).at(0) - (*xy).at(0);
	      const double dy = (*xyn).at(1) - (*xy).at(1);
	      a = sqrt( dx*dx + dy*dy );
	    }
	    {
	      const double dx = (*xynn).at(0) - (*xyn).at(0);
	      const double dy = (*xynn).at(1) - (*xyn).at(1);
	      b = sqrt( dx*dx + dy*dy );
	    }
	    {
	      const double dx = (*xynn).at(0) - (*xy).at(0);
	      const double dy = (*xynn).at(1) - (*xy).at(1);
	      c = sqrt( dx*dx + dy*dy );
	    }
	    // Kahan, W. 1986. Calculating Area and Angle of a Needle-like
	    // Triangle, unpublished manuscript
	    // Goldberg, David, What Every Computer Scientist Should Know about
	    // Floating-Point arithmetic, ACM Computing Surveys, Vol. 23, #1,
	    // March 1991, pp. 5-48
	    {
	      // order a, b, c such as a >= b >= c
#define ORDER(A,B) if (A<B) { const double t_ = A; A = B; B = t_; }
	      ORDER(a,b)
		ORDER(b,c)
		ORDER(a,b)
#undef ORDER
		assert(a>=b && b>=c);
	      area += sqrt(
			   (a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
			   )/4;
	    }
	    vcI.Prev();
	  }
      }
    }
  varea = area;
  if( varea<=0.0 ) { // debug
    std::cout << "Error: zero or negative varea = " << varea << std::endl;
    std::cout << "Node: " << id << " " << x << " " << y << " "
	 << BoundName(boundary) << std::endl;
    getEdg()->TellCoords();
  }
  assert( varea>0.0 );
  varea_rcp = 1.0/varea;

  // debug
  if (0){ //DEBUG
    if( id==83 ) {
      std::cout << " reading list: ";
      for( ce = vtxIter.FirstP(); !(vtxIter.AtEnd()); ce = vtxIter.NextP() )
	{
	  tArray2<double> const &xy = ce->getRVtx();
	  std::cout << xy.at(0) << " " << xy.at(1) << "; ";
	}
      std::cout << std::endl;
      std::cout << "reading spokes: ";
      ce = getEdg();
      do
	{
	  assert( ce!=0 );
	  tArray2<double> const &xy = ce->getRVtx();
	  std::cout << xy.at(0) << " " << xy.at(1) << "; ";
	  ce = ce->getCCWEdg();
	} while( ce != getEdg() );
      std::cout << std::endl;
    }
  }

  return area;
}

/*******************************************************************\
**
**  tNode::getVoronoiVertexList
**
**  Creates and returns a list of (x,y) coordinates for the
**  Voronoi vertices associated with the node. The list is
**  created by moving around the spokes and adding each spoke's
**  right-hand Voronoi vertex to the back of the list.
**    A pointer to the vertex list is passed as a parameter; any
**  prior contents are flushed before the list of points is created.
**
**    Created: 2/2000 GT
**
\*******************************************************************/
void tNode::getVoronoiVertexList( tList<Point2D> * vertexList )
{
   assert( !boundary );

   vertexList->Flush();

   // Loop around spokes, adding the right-hand Voronoi vertex of
   // each to the list, until we've gone all the way around
   tEdge *ce = edg;
   do
   {
      tArray2<double> const & vtxarr = ce->getRVtx();
      const Point2D vtx( vtxarr.at(0), vtxarr.at(1) );
      //std::cout << "ADDING TO LIST: x " << vtx.x << " y " << vtx.y << std::endl;
      vertexList->insertAtBack( vtx );
      ce = ce->getCCWEdg();
   }
   while( ce!=edg );

   assert( vertexList->getSize()!=0 );
}
/*******************************************************************\
**
**  tNode::getVoronoiVertexXYZList
**
**  Creates and returns a list of (x,y,z) coordinates for the
**  Voronoi vertices associated with the node. The list is
**  created by moving around the spokes and adding each spoke's
**  right-hand Voronoi vertex to the back of the list. The z
**  coordinate is obtained by linear interpolation from the 3
**  points of the triangle of which the vertex is the circumcenter.
**    A pointer to the vertex list is passed as a parameter; any
**  prior contents are flushed before the list of points is created.
**
**    Created: 2/2000 GT
**
\*******************************************************************/
void tNode::getVoronoiVertexXYZList( tList<Point3D> * vertexList )
{
   tNode *n1, *n2;

   assert( !boundary );

   vertexList->Flush();

   // Loop around spokes, adding the right-hand Voronoi vertex of
   // each to the list, until we've gone all the way around
   tEdge *ce = edg;
   n2 = ce->getDestinationPtrNC();
   do
   {
      ce = ce->getCCWEdg();
      n1 = n2;
      n2 = ce->getDestinationPtrNC();
      tArray2<double> const &vtxarr = ce->getRVtx();
      const tArray<double> zvals( this->z, n1->getZ(), n2->getZ());
      const double zz =
	PlaneFit(vtxarr.at(0), vtxarr.at(1), this->get2DCoords(),
		 n1->get2DCoords(), n2->get2DCoords(), zvals );
      const Point3D vtx(vtxarr.at(0), vtxarr.at(1), zz);
      if (0) //DEBUG
	std::cout << "ADDING TO LIST: x " << vtx.x << " y " << vtx.y << " z " << vtx.z << std::endl;
      vertexList->insertAtBack( vtx );
   }
   while( ce!=edg );

   assert( vertexList->getSize()!=0 );
}

/*******************************************************************\
**
**  tNode::ConvertToClosedBoundary
**
**  Makes the node into a closed boundary by setting its boundary
**  status flag. The function also updates the boundary ("flow
**  allowed") status of the node's spokes and their complements.
**
\*******************************************************************/
void tNode::ConvertToClosedBoundary()
{
   tEdge *ce;   // an edge and its complement

   // Reset boundary flag
   boundary = kClosedBoundary;

   // Signal all connected edges and their complements to become no-flow
   ce = edg;
   do
   {
      assert( ce!=0 );
      if( ce->getBoundaryFlag()!=kNonBoundary )
      {
         ce->setFlowAllowed( tEdge::kFlowNotAllowed );
         // get complement and change it too TODO
      }

   } while( (ce=ce->getCCWEdg()) != edg );

}

/**************************************************************************\
**
**  TellAll
**
**  ...is a debugging routine that prints out a bunch of info about a node.
**
\**************************************************************************/
#ifndef NDEBUG
void tNode::TellAll() const
{
   std::cout << " NODE " << id << ":\n";
   std::cout << "  x=" << x << " y=" << y << " z=" << z;
   std::cout << "  boundary: " << BoundName(boundary)
        << "\n  varea: " << varea << std::endl;
   if( edg )
       std::cout << "  points to edg #" << edg->getID() << std::endl;
   else std::cout << "  edg is undefined!\n";

}
#endif


/**************************************************************************\
\***  Functions for class tEdge  ******************************************/

/**************************************************************************\
**
**  tEdge::TellCoords
**
**  Debugging routine that reports edge ID and coords of endpoints.
**
\**************************************************************************/
void tEdge::TellCoords()
{
   std::cout << "EDGE " << id << ":\n";
   std::cout << "  " << org->getID() << " (" << org->getX() << ","
        << org->getY() << ") -> " << dest->getID() << " ("
        << dest->getX() << "," << dest->getY() << ")" << std::endl;
}



/**************************************************************************\
**
**  tEdge::FindComplement
**
**  Finds and returns the edge's complement edge. Does this by checking
**  the spokes connected to its destination node and returning the one
**  that connects back to its origin.
**
**  Data mbrs modified:  none
**  Returns:  ptr to the complement edge
**  Assumes:  valid connectivity (one of destination's spokes connects
**            back to origin)
**  Created: 4/29/98 GT
**
\**************************************************************************/
tEdge * tEdge::FindComplement()
{
   assert( org!=0 && dest!=0 && dest->getEdg()!=0 );

   tEdge * ce = dest->getEdg();
   while( ce->getDestinationPtrNC() != org )
   {
      ce = ce->getCCWEdg();
      assert( ce!=0 && ce!=dest->getEdg() );
   }
   return ce;
   // TODO: test for infinite loop using assert

}

/**************************************************************************\
**
**  tEdge::CheckConsistency
**
**  Check for valid origin, destination, and ccwedg
**  Previously in tMesh.
**
**  AD - March 2004
**
\**************************************************************************/
bool tEdge::CheckConsistency(){
  tNode * org, * dest;

  if( (org=getOriginPtrNC() ) == NULL)
    {
      std::cerr << "EDGE #" << getID()
	   << " does not have a valid origin point\n";
      return false;
    }
  if( (dest=getDestinationPtrNC() ) == NULL)
    {
      std::cerr << "EDGE #" << getID()
	   << " does not have a valid destination point\n";
      return false;
    }
  if( (ccwedg=getCCWEdg() ) == NULL)
    {
      std::cerr << "EDGE #" << getID()
	   << " does not point to a valid counter-clockwise edge\n";
      return false;
    }
  if( ccwedg->getOriginPtrNC()!=org )
    {
      std::cerr << "EDGE #" << getID()
	   << " points to a CCW edge with a different origin\n";
      return false;
    }
  if( ccwedg->getDestinationPtrNC()==dest )
    {
      std::cerr << "EDGE #" << getID()
	   << " points to a CCW edge with the same destination\n";
      return false;
    }
  tEdge *cwedg;
  if( (cwedg=getCWEdg() ) == NULL)
    {
      std::cerr << "EDGE #" << getID()
	   << " does not point to a valid clockwise edge\n";
      return false;
    }
  if( cwedg->getOriginPtrNC()!=org )
    {
      std::cerr << "EDGE #" << getID()
	   << " points to a CW edge with a different origin\n";
      return false;
    }
  if( cwedg->getDestinationPtrNC()==dest )
    {
      std::cerr << "EDGE #" << getID()
	   << " points to a CCW edge with the same destination\n";
      return false;
    }
  if ( getCCWEdg()->getCWEdg() != this ||
       getCWEdg()->getCCWEdg() != this)
    {
      std::cerr << "EDGE #" << getID()
	   << ": inconsistency related to pointers to CW and CCW edges\n";
      return false;
    }

  if( org==dest )
    {
      std::cerr << "EDGE #" << getID()
	   << " has the same origin and destination nodes\n";
      return false;
    }
  return true;
}

//***************************************************************************
// InitializeEdge: Function to initialize edge by setting its origin
//   to n1, dest to n2, and putting it in spoke list of n1. The third node,
//   n3, is the third node of the triangle, and, so, n1, n2, and n3 are
//   _clockwise_ rather than ccw as are the similar arguments in
//   tMesh::AddEdge(...)
// 3/99 SL
// 4/2003 AD
// 8/2003 SL: added flag for "usefuturePosn" so that if geometry is used to 
//   place edge in spoke list, new coordinates can be used if nodes are 
//   moving. Has default value of "false", so existing calls need not
//   be changed.
//***************************************************************************
void tEdge::InitializeEdge( tNode* n1, tNode* n2, tNode const * n3, bool useFuturePosn )
{
   assert( n1!=0 && n2!=0 && n3!=0 );
   setOriginPtr( n1 );
   setDestinationPtr( n2 );
   setFlowAllowed(n1, n2);
   len = slope = vedglen = 0.0;
   tSpkIter sI( n1 );
   //add pointers to the new edge to node's spokeLists:
   if( sI.isEmpty() )
       sI.insertAtFront( this );
   else if( sI.ReportNextP() == sI.CurSpoke() )
       sI.insertAtFront( this );
   else
   {
      tEdge* ce;
      for( ce = sI.FirstP();
           ce->getDestinationPtr() != n3 && !( sI.AtEnd() );
           ce = sI.NextP() ) ;
      //make sure we found the right spoke; if not:
      if( sI.AtEnd() )
      {
         for( ce = sI.FirstP();
              !( sI.AtEnd() );
              ce = sI.NextP() )
         {
            if( !useFuturePosn )
            {
                if( PointsCCW( UnitVector( ce ),
                               UnitVector( this ),
                               UnitVector( sI.ReportNextP() )
                               )
                    )
                    break;
            }
            else
            {
               if( PointsCCW( NewUnitVector( ce ),
                              NewUnitVector( this ),
                              NewUnitVector( sI.ReportNextP() )
                              )
                   )
                   break;
            }
         }
      }
      //put edge2 in SPOKELIST:
      sI.insertAtNext( this );
   }
}

/**************************************************************************\
**
**  Functions for class tTriangle.
**
**  For a description of tTriangle objects, see meshElements.h.
**
**  Modifications:
**    - added tTriangle::FindCircumcenter() 1/11/98 gt
**
\**************************************************************************/

//***************************************************************************
// InitializeTriangle: Function to initialize triangle by setting node ptrs to
//   n0, n1, n2; edge ptrs to appropriate edges; and nbr triangles. Also
//   resets edges' tri ptrs. Nodes in argument are given in ccw order.
// 3/99 SL
// 4/2003 AD
//***************************************************************************
void tTriangle::InitializeTriangle( tNode* n0, tNode* n1, tNode* n2 )
{
   assert( n0 != 0 && n1 != 0 && n2 != 0 );
   p[0] = n0;
   p[1] = n1;
   p[2] = n2;
   // setEPtr also sets tri ptr of edge:
   setEPtr( 0, n0->EdgToNod( n2 ) );
   setEPtr( 1, n1->EdgToNod( n0 ) );
   setEPtr( 2, n2->EdgToNod( n1 ) );
   // Now we assign the neighbor triangle pointers. The loop successively
   // gets the spokelist for (p0,p1,p2) and sets cn to the next ccw point
   // (p1,p2,p0). It then finds the edge (spoke) that joins the two points
   // (p0->p1, p1->p2, p2->p0). These are the edges that are shared with
   // neighboring triangles (t2,t0,t1) and are pointed to by the neighboring
   // triangles. This means that in order to find neighboring triangle t2,
   // we need to find the triangle that points to edge (p0->p1), and so on.
   // In general, t((j+2)%3) is the triangle that points to edge
   // p(j)->p((j+1)%3).
   for( int j=0; j<3; j++ )
   {
      // Find edge ce that connects p(j)->p(j+1)
      tEdge* ce = p[j]->EdgToNod( p[(j+1)%3] );
      // Find the triangle, if any, that shares (points to) this edge
      // and assign it as the neighbor triangle t((j+2)%3).
      tTriangle* nbrtriPtr = ce->TriWithEdgePtr();
      t[(j+2)%3] = nbrtriPtr;  //set tri TRI ptr (j+2)%3
      // If a neighboring triangle was found, tell it that the current
      // new triangle is its neighbor too. We need to tell it which
      // neighbor we are (0, 1, or 2), and the mapping is like this:
      // if the nbr tri calls the shared edge (0,1,2) then we are its
      // nbr (1,2,0). (ie, tri_number = (edg_number+1)%3 )
      if( nbrtriPtr != 0 )
      {
	int i;
	for( i=0; i<3; ++i )
	  {
            assert( nbrtriPtr->e[i] != 0 );
            assert( ce != 0 );
            if( nbrtriPtr->e[i] == ce ) goto found;
	  }
	assert( 0 );
	::abort();
      found:
	nbrtriPtr->t[(i+1)%3] = this;  //set NBR TRI ptr to tri
      }
   }
}

//destructor
tTriangle::~tTriangle()
{
   for( int i=0; i<3; i++ )
   {
      p[i] = 0;
      setEPtr( i, 0 );
      t[i] = 0;
   }
}

/*****************************************************************************\
**
**  tTriangle::FindCircumcenter
**
**  Finds the circumcenter of the triangle by finding the intersection of
**  the perpendicular bisectors of sides (p0,p1) and (p0,p2). Returns the
**  coordinates of the circumcenter as a 2-element array. Note that the
**  circumcenter is also the Voronoi cell vertex associated with the
**  triangle's three nodes (that's the point of computing it).
**
\*****************************************************************************/
tArray2< double >
tTriangle::FindCircumcenter() const
{
   assert( pPtr(0) && pPtr(1) && pPtr(2) );

   // Coordinates of triangle's nodes p0, p1, and p2
   tArray2< double > xyo; pPtr(0)->get2DCoords( xyo );
   tArray2< double > xyd1; pPtr(1)->get2DCoords( xyd1 );
   tArray2< double > xyd2; pPtr(2)->get2DCoords( xyd2 );

   // Find the midpoints of the two sides (p0,p1) and (p0,p2) and store them
   // in (x1,y1) & (x2,y2). Then get the distance between p0 and the
   // midpoints of each side
   const double x1 = (xyo.at(0) + xyd1.at(0)) / 2;
   const double y1 = (xyo.at(1) + xyd1.at(1)) / 2;
   const double x2 = (xyo.at(0) + xyd2.at(0)) / 2;
   const double y2 = (xyo.at(1) + xyd2.at(1)) / 2;
   const double dx1 = x1-xyo.at(0);
   const double dy1 = y1-xyo.at(1);
   const double dx2 = x2-xyo.at(0);
   const double dy2 = y2-xyo.at(1);

   double XX, YY;

   // Compute the intercept of the bisectors of the two sides:
   // Case: neither spoke is horizontal (ok to divide by dy1 and dy2)
   if( fabs(dy1)>0 && fabs(dy2)>0 )
   {
      assert( dy1!=0 && dy2!=0 );
      const double m1 = -dx1/dy1;
      const double m2 = -dx2/dy2;
      assert( m1!=m2 ); // should never happen; means edges are parallel
      XX = (y2 - m2 * x2 - y1 + m1 * x1) / (m1 - m2);
      YY = m1 * (XX - x1) + y1;
   }
   // Case: one spoke is horizontal, but neither are vertical
   else if( dx1!=0 && dx2!=0 )
   {
      assert( dx1!=0 && dx2!=0 );
      const double m1 = dy1/dx1;
      const double m2 = dy2/dx2;
      assert( m1!=m2 );
      YY = (m1 * y1 + x1 - m2 * y2 - x2) / (m1 - m2);
      XX = -YY * m1 + m1 * y1 + x1;
   }
   // Special case: one is vertical, the other horizontal
   else
   {
      if( dx1!=0 )
      {
         assert( dx2==0 && dy1==0 );
	 XX = xyo.at(0) + dx1;
	 YY = xyo.at(1) + dy2;
      }
      else
      {
         assert( dx1==0 && dy2==0 );
	 XX = xyo.at(0) + dx2;
	 YY = xyo.at(1) + dy1;
      }
   }
   return tArray2< double >(XX, YY);
}


/*****************************************************************************\
**
**  tTriangle::SetIndexByOrder
**
**  Build a circular array so that pPtr(index[0])->getID() is the lowest
**
\*****************************************************************************/
void tTriangle::SetIndexIDOrdered()
{
  const int ID[] = { pPtr(0)->getID(),
		     pPtr(1)->getID(),
		     pPtr(2)->getID() };
  // find the vertex with the lowest ID, can be 0 or 1 or 2
  int ID_;
  index_[0] = 0; ID_ = ID[0];
  if (ID[1] < ID_ ) { index_[0] = 1; ID_ = ID[1]; }
  if (ID[2] < ID_ ) { index_[0] = 2; }
  // complete the sequence
  index_[1] = (index_[0]+1)%3;
  index_[2] = (index_[1]+1)%3;
}

/*****************************************************************************\
**
**  tTriangle::NbrToward
**
**  Geometric calculations to find triangle neighbor
**
\*****************************************************************************/
tTriangle* tTriangle::NbrToward( double x, double y )
{
  tTriangle* ct = 0;
  tArray2< double > xy0(x, y);
  tArray2< double > xy1;
  tArray2< double > xy2;
  tArray2< double > xy3;
  p[0]->get2DCoords( xy1 );
  p[1]->get2DCoords( xy2 );
  p[2]->get2DCoords( xy3 );
  // use Predicates::orient2d:
  double c = predicate.orient2d( xy1.getArrayPtr(), xy2.getArrayPtr(),
				 xy0.getArrayPtr() );
  if( c < 0.0 )
    return ct = t[2];
  c = predicate.orient2d( xy2.getArrayPtr(), xy3.getArrayPtr(), xy0.getArrayPtr() );
  if( c < 0.0 )
    return ct = t[0];
  c = predicate.orient2d( xy3.getArrayPtr(), xy1.getArrayPtr(), xy0.getArrayPtr() );
  if( c < 0.0 )
    return ct = t[1];
  return ct;
}

/*****************************************************************************\
**
**  tTriangle::containsPoint
**
**  Does the current triangle contain the point (x,y)
**
\*****************************************************************************/
bool tTriangle::containsPoint(double x, double y) const
{
  const double P[] = {x, y};
  const double P0[] = {p[0]->getX(), p[0]->getY()};
  const double P1[] = {p[1]->getX(), p[1]->getY()};
  if (predicate.orient2d(P0, P1, P) <0)
    return false;
  const double P2[] = {p[2]->getX(), p[2]->getY()};
  if (predicate.orient2d(P1, P2, P) <0)
    return false;
  if (predicate.orient2d(P2, P0, P) <0)
    return false;
  return true;
}

/* TellAll: debugging output routine */
#ifndef NDEBUG
void tTriangle::TellAll() const
{
   int i;

   assert( this!=0 );
   std::cout << "TRIANGLE #" << id << ":\n";
   for( i=0; i<3; i++ )
   {
      std::cout << "  P" << i << " ";
      if( p[i]!=0 ) std::cout << p[i]->getID() << " (" << p[i]->getX() << ","
                         << p[i]->getY() << ")";
      else std::cout << "(ndef)";
      std::cout << "  E" << i << " ";
      if( e[i]!=0 ) std::cout << e[i]->getID();
      else std::cout << "(ndef)";
      std::cout << "  T" << i << " ";
      if( t[i]!=0 ) std::cout << t[i]->getID();
      else std::cout << "(ndef)";
      std::cout << std::endl;
   }
}
#endif


