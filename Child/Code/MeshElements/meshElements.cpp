/**************************************************************************\
**
**  meshElements.cpp: Functions for mesh element classes tNode, tEdge, and
**                    tTriangle. (formerly called gridElements.cpp)
**
**  For a description of these 3 classes, see meshElements.h.
**
**  Modifications:
**   - gt added code to recreate the "edg" implementation as an alternative
**     to spokelist. Dec 1 '97.
**   - previously separate tNode, tEdge, and tTriangle files grouped into
**     "gridElements", 1/20/98 gt
**   - added tNode::AttachNewSpoke and tEdge::WelcomeCCWNeighbor gt 2/99
**
**  $Id: meshElements.cpp,v 1.33 1999-09-09 21:31:07 gtucker Exp $
\**************************************************************************/

#include <assert.h>
#include <math.h>
#include "meshElements.h"

int PointsCCW( tArray< double > &, tArray< double > &, tArray< double > & );


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
tArray< double > FindIntersectionCoords( tArray< double > xy1,
                                         tArray< double > xy2,
                                         tArray< double > xy3,
                                         tArray< double > xy4 )
{
   double dxa, dxb, dya, dyb, a, b, c, f, g, h;
   //Xx, y;
   tArray< double > intxy(2);

   dxa = xy2[0] - xy1[0];
   dxb = xy4[0] - xy3[0];
   dya = xy2[1] - xy1[1];
   dyb = xy4[1] - xy3[1];
   a = dya;
   b = -dxa;
   c = dxa * xy1[1] - dya * xy1[0];
   f = dyb;
   g = -dxb;
   //h = dxb * xy3[1] - dyb * xy4[0];
   h = dxb * xy3[1] - dyb * xy3[0];  // seems to be a bug above; fixed here?
   if( fabs(dxa) > 0 && fabs(dxb) > 0 )
   {
      if( fabs(f - g * a / b) > 0 )
      {
         intxy[0] = (g * c / b - h) / (f - g * a / b);
         intxy[1] = (-c - a * intxy[0]) / b;
      }
   }
   else
   {
      if( fabs(dya) > 0 && fabs(dyb) > 0 )
      {
         if( fabs(g - f * b / a) > 0 )
         {
            intxy[1] = (f * c / a - h) / (g - f * b / a);
            intxy[0] = (-c - b * intxy[1]) / a;
         }
      }
      else //one horiz. line and one vert. line:
      {
         if( fabs(dya) == 0 )
         {
            intxy[0] = xy3[0];
            intxy[1] = xy1[1];
         }
         else
         {
            intxy[0] = xy1[0];
            intxy[1] = xy3[1];
         }
      }
   }
   return intxy;
}



/***************************************************************************\
\**  Functions for class tNode  ********************************************/

/***********************************************************************\
**
**  Constructors & destructors:
**
**  Default:  initializes values to zero.
**  Copy:  copies all values and makes duplicate spoke list
**  Destructor:  no longer used
**
\***********************************************************************/

//default constructor
tNode::tNode()                                                      //tNode
{
   id = boundary = 0;
   x = y = z = varea = varea_rcp = 0.0;
   edg = 0;
   
     //cout << "tNode()" << endl;
}

//copy constructor
tNode::tNode( const tNode &original )                               //tNode
{
   if( &original != 0 )
   {
      id = original.id;
      x = original.x;
      y = original.y;
      z = original.z;
      boundary = original.boundary;
      varea = original.varea;
      varea_rcp = original.varea_rcp;
      edg = original.edg;
      if( &(original.spokeList) != 0 )
      {
         for( int i=0; i<original.spokeList.getSize(); i++ )
         {
            insertBackSpokeList( original.spokeList.getIthPtrNC( i ) );
         }
      }
        //else spokeList = 0;
   }
     //cout << "tNode( original )" << endl;
}

tNode::~tNode()                                                      //tNode
{
     //if( spokeList != 0 ) delete spokeList;
     //cout << "    ~tNode()" << endl;
}


/***********************************************************************\
**
**  Overloaded operators:
**
**    assignment: copies all values (spokelist's assignment operator
**                creates duplicate copy of list)
**    right shift: takes input for x, y and z values from input stream
**    left shift: sends the following data to the output stream:
**                node ID, x, y, z values, and IDs of neighboring nodes,
**                which are obtained through the spokelist.
**
\***********************************************************************/

//assignment
const tNode &tNode::operator=( const tNode &right )                  //tNode
{
   if( &right != this )
   {
      tPtrListIter< tEdge > spokIter;
      //XtEdge *ce;
      id = right.id;
      x = right.x;
      y = right.y;
      z = right.z;
      boundary = right.boundary;
      varea = right.varea;
      varea_rcp = right.varea_rcp;
      edg = right.edg;
      spokeList = right.spokeList;
   }
   return *this;
}

//right shift
istream &operator>>( istream &input, tNode &node )                   //'tNode'
{
   cout << "x y z:" << endl;
   input >> node.x >> node.y >> node.z;
   return input;
}

//left shift
ostream &operator<<( ostream &output, tNode &node )            //'tNode'
{
   //tPtrListIter< tEdge > spokeIter( node.getSpokeListNC() );
   
   output << node.id << ": " << node.x << " " << node.y << " "
          << node.z << ";";
   output <<
       node.spokeList.getFirstNC()->getPtrNC()->getDestinationPtrNC()->getID()
   //output << spokeIter.DatPtr()->getDestinationPtrNC()->getID()
          << " ";
   output << endl;
   return output;
}


/***********************************************************************\
**
**  tNode "get" functions:
**
**  get3DCoords - returns x, y, z as a 3-element array
**  get2DCoords - returns x & y coords as a 2-element array
**  getID - returns ID #
**  getX - returns node's x coord
**  getY - returns node's y coord
**  getZ - returns node's z value
**  getVArea - returns Voronoi area
**  getVArea_Rcp - returns 1 / Voronoi area
**  getBoundaryFlag - returns boundary code
**  getEdg - returns pointer to one spoke
**  getSpokeList - returns const reference to spoke list
**  getSpokeListNC - returns non-const reference to spoke list
**  getFirstSpokeNode - returns ptr to 1st spokelist item 
**                      (return type tPtrListNode *)
**  getFirstSpokeNodeNC - non-const version of the above
**  getNextSpokeNode - returns ptr to next spokelist item that follows
**                       prevedg
**  getNextSpokeNodeNC - non-const version of the above
**
\***********************************************************************/

tArray< double >                                                       //tNode
tNode::get3DCoords() const
{
   tArray< double > xyz(3);
   xyz[0] = x;
   xyz[1] = y;
   xyz[2] = z;
   return xyz;
}

tArray< double >                                                   //tNode
tNode::get2DCoords() const
{
   tArray< double > xy(2);
   xy[0] = x;
   xy[1] = y;
   return xy;
}

int tNode::getID() const {return id;}                    //tNode
double tNode::getX() const {return x;}
double tNode::getY() const {return y;}
double tNode::getZ() const {return z;}
double tNode::getVArea() const {return varea;}             //tNode
double tNode::getVArea_Rcp() const {return varea_rcp;}     //tNode
int tNode::getBoundaryFlag() const {return boundary;}      //tNode
tEdge * tNode::getEdg() {return edg;}

const tPtrList< tEdge > &                                     //tNode
tNode::getSpokeList() const {assert( &spokeList != 0 ); return spokeList;}

tPtrList< tEdge > &                                           //tNode
tNode::getSpokeListNC() {assert( &spokeList != 0 ); return spokeList;}

const tPtrListNode< tEdge > *                                 //tNode
tNode::getFirstSpokeNode() const {return spokeList.getFirst();}

tPtrListNode< tEdge > *                                       //tNode
tNode::getFirstSpokeNodeNC() {return spokeList.getFirstNC();}

const tPtrListNode< tEdge > *tNode::                          //tNode
getNextSpokeNode( const tPtrListNode< tEdge > *prevedg ) const
{return prevedg->getNext();}

tPtrListNode< tEdge > *tNode::                                //tNode
getNextSpokeNodeNC( tPtrListNode< tEdge > *prevedg ) const
{return prevedg->getNextNC();}


/***********************************************************************\
**
**  tNode "set" functions:
**
**  setID - sets ID number to val
**  setX - sets x coord to val
**  setY - sets y coord to val
**  setZ - sets z value to val
**  setVArea - sets Voronoi area to val
**  setVArea_Rcp - sets 1/Voronoi area to val
**  setBoundaryFlag - returns boundary code
**  set3DCoords - sets x, y, z to val1, val2, val3
**  set2DCoords - sets x & y to val1 and val2
**  setEdg - sets edge ptr to theEdg
**
**  Note: unless otherwise noted, no runtime value checking is done
**        (aside from assert statements)
**
\***********************************************************************/

void tNode::setID( int val ) {id = val;}                             //tNode
void tNode::setX( double val ) {x = val;}                             //tNode
void tNode::setY( double val ) {y = val;}                             //tNode
void tNode::setZ( double val ) {z = val;}                             //tNode

void tNode::setVArea( double val )                                    //tNode
{
  assert( val>=0.0 );
  varea = val;
  /*varea = ( val >= 0.0 ) ? val : 0.0;*/
}

void tNode::setVArea_Rcp( double val )                                //tNode
{
  assert( val>=0.0 );
  varea_rcp = val;
  /*varea_rcp =  ( val >= 0.0 ) ? val : 0.0;*/
}

void tNode::setBoundaryFlag( int val )                               //tNode
{
  assert( val>=0 && val<=2 );
  boundary = val;
  /*boundary = (val >=0 && val <= 2) ? val : 0;*/
}

void tNode::set2DCoords( double val1, double val2 )                    //tNode
{
   setX( val1 );
   setY( val2 );
}

void tNode::set3DCoords( double val1, double val2, double val3 )        //tNode
{
   setX( val1 );
   setY( val2 );
   setZ( val3 );
}

void tNode::setEdg( tEdge * theEdg )
{
   assert( theEdg > 0 );
   edg = theEdg;
   //cout << "Assigning edge " << theEdg->getID() << " to node " << getID() << endl;
}


/***********************************************************************\
**
**  tNode::ChangeZ:  Adds delz to current z value
**
\***********************************************************************/
void tNode::ChangeZ( double delz ) { z += delz; }                    //tNode


/***********************************************************************\
**
**  tNode::insertFrontSpokeList
**  tNode::insertBackSpokeList
**
**  Places eptr at the front or back of the spoke list (respectively)
**  and makes the list circular (is the latter step necessary? TODO).
**
\***********************************************************************/
void tNode::insertFrontSpokeList( tEdge *eptr )                      //tNode
{
   spokeList.insertAtFront( eptr );
   assert( spokeList.getFirst() != 0 );
   makeWheel();
}

void tNode::insertBackSpokeList( tEdge *eptr )                       //tNode
{
   spokeList.insertAtBack( eptr );
   assert( spokeList.getFirst() != 0 );
   makeWheel();
}


/***********************************************************************\
**
**  tNode::makeWheel:  makes the spoke list circular
**
\***********************************************************************/
void tNode::makeWheel() {spokeList.makeCircular();}                  //tNode

const tEdge *tNode::NextSpoke( tPtrListNode< tEdge > * current )//tNode
{
   if( current == 0 || current == spokeList.getLast() )
       current = spokeList.getFirstNC();
   else current = current->getNextNC();
   return current->getPtr();
}


/*****************************************************************************\
**
**  tNode::AttachFirstSpoke
**
**  Attaches the first spoke to the node by pointing edg to that spoke,
**  and then telling the spoke to point to itself. thespoke is the edge
**  being added.
**
**      Data members updated: edg, thespoke's ccw edge
**      Called by: tMesh::AddEdge
**      Calls: (none)
**      Created: 2/4/99 GT
**
\*****************************************************************************/
void tNode::AttachFirstSpoke( tEdge *thespoke )
{
   assert( thespoke!=0 );
   assert( thespoke->getOriginPtr()==this );
   //assert( edg==0 );
   edg = thespoke;
   thespoke->setCCWEdg( thespoke );
}


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
double tNode::Dist( tNode * p0, tNode * p1 )
{
  double a,b,c,res;

  a=(p1->y)-(p0->y);  
  b=-((p1->x)-(p0->x));
  c=-((a*(p0->x))+(b*(p0->y)));
  res=(a*x + b*y + c) / sqrt(a*a + b*b);
  if (res<0) res=-res;
  return(res);
}


/*****************************************************************************\
**
**  tNode::EdgToNod
**
**  Finds and returns the spoke (edge) that connects the current node to _nod_,
**  or zero if no such spoke is found.
**
\*****************************************************************************/
tEdge *tNode::EdgToNod( tNode * nod )
{
   tPtrListIter< tEdge > spokIter( this->spokeList );
   tEdge * ce;

   for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
   {
      if( ce->getDestinationPtr()->getID() == nod->getID() ) return ce;
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
   int cw;
   double area = 0;
   double a, b, c, dx, dy, dx0, dx1, dy0, dy1, dx2, dy2;
   double vx, vy, x0, y0, x1, y1, x2, y2, m1, m2;
   tEdge * ce, *edgptr;
   tPtrList< tEdge > vedgList /*= spokeList*/;
   tPtrListIter< tEdge > //spokIter( spokeList ),
       vtxIter( vedgList );
   tList< tArray< double > > vcL; // list of vertex coordinates
   tListIter< tArray< double > > vcI( vcL ); // iterator for coord list
   tArray< double > xy, xyn, xynn, xynnn, xy1, xy2, xy3, xy4;
   int i;

   // Create a duplicate list of edges; we will modify this list to obtain
   // the correct vertices. In some cases, we may need to delete an edge
   // to get the correct list of vertices; we don't want to delete the
   // spoke ptr, so we make a duplicate list.
   //if( id==83 ) cout << "NODE 83: " << x << "," << y << endl;
   ce = edg;
   do
   {
      assert( ce>0 );
      vedgList.insertAtBack( ce );
      //xy = ce->getRVtx();
      //cout << xy[0] << " " << xy[1] << "; " << flush;
      //xy = vedgList.getLast()->getPtrNC()->getRVtx();
      //cout << xy[0] << " " << xy[1] << endl << flush;
      //if( id==83) cout << " " << ce->getDestinationPtr()->getX() << ","
      //                 << ce->getDestinationPtr()->getY() << endl;
      ce = ce->getCCWEdg();
   } while( ce != edg );
   vedgList.makeCircular();
   //cout << endl << flush;
   // Check boundary status: Voronoi area only defined for non-boundary nodes
   if( boundary == kNonBoundary )
   {
      cw = TRUE;
      //cout << "find clockwise loops" << endl << flush;
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
         cw = FALSE;
         for( ce=vtxIter.FirstP(); !( vtxIter.AtEnd() ); ce=vtxIter.NextP() )
         {
            xy = ce->getRVtx();
            xyn = vtxIter.NextP()->getRVtx();
            xynn = vtxIter.NextP()->getRVtx();
            dx0 = xynn[0] - xyn[0];
            dy0 = xynn[1] - xyn[1];
            dx1 = xy[0] - xyn[0];
            dy1 = xy[1] - xyn[1];
            if( dy0 * dx1 > dx0 * dy1 ) //clockwise
            {
               xynnn = vtxIter.NextP()->getRVtx();
               dx0 = xynnn[0] - xynn[0];
               dy0 = xynnn[1] - xynn[1];
               dx1 = xyn[0] - xynn[0];
               dy1 = xyn[1] - xynn[1];
               if( dy0 * dx1 > dx0 * dy1 ) //clockwise
               {
                  //two consecutive clockwise vertices=>want intersection
                  //of bisectors of ce->nextedg and
                  //ce->nextedg->nextedg->nextedg:
                  cw = TRUE;
                  x0 = x; //node.x
                  y0 = y; //node.y
                  xy1 = ce->getDestinationPtr()->get2DCoords();
                  //vtxIter.Prev();
                  xy2 = vtxIter.PrevP()->getDestinationPtr()->get2DCoords();
                  x1 = ( x0 + xy1[0] ) / 2;
                  y1 = ( y0 + xy1[1] ) / 2;
                  x2 = ( x0 + xy2[0] ) / 2;
                  y2 = ( y0 + xy2[1] ) / 2;
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
                  xyn[0] = vx;
                  xyn[1] = vy;
                  dx = xy[0] - vx;
                  dy = xy[1] - vy;
                  //cout << "reset vedglen and rvtx for edge "
                  //     << edgptr->getID() << " to len "
                  //     << sqrt( dx*dx + dy*dy )
                  //     << ", x, y, " << xyn[0] << ", " << xyn[1] << endl << flush;
                    //reset 'next' edge's vertex to newly found intersection,
                    //length adjusted accordingly
                  edgptr->setVEdgLen( sqrt( dx*dx + dy*dy ) );
                  edgptr->setRVtx( xyn );
                  edgptr = vtxIter.ReportNextP();
                  //cout << "reset vedglen and rvtx for edge "
                  //     << edgptr->getID()
                  //     << " to len 0.0, x, y, " << xynnn[0] << ", "
                  //     << xynnn[1] << endl << flush;
                    //reset 'next-next' edge's vertex to the coordinates
                    //of the 'next-next-next' edge's vertex; length to zero
                  edgptr->setVEdgLen(0.0);
                  edgptr->setRVtx( xynnn );
                  edgptr = 0;
                    //delete the offending vertex's edge from list
                  vedgList.removeNext( edgptr, vtxIter.NodePtr() );
               }
            }
            vtxIter.Get( ce->getID() );
         }
      } while( cw ); //while we're still finding loops in the polygon

      //Before the next step, make a list of V. vertex coord. arrays.
      //In doing so, check for parts of the V. area lying outside the
      //mesh domain and cut them off by substituting coordinates of
      //two intersections of V. edges with boundary edge for the V.
      //vertex lying outside the boundary. This should take care of any
      //outlying area as long as all boundaries are convex.
      // Go through spokes and put RVtx of ccw edge in coord list, but
      // first check that the vtx lies within the bndies
      tEdge *ne, *nne;
      tNode *bn0, *bn1;
      for( ce = vtxIter.FirstP(); !(vtxIter.AtEnd()); ce = vtxIter.NextP() )
      {
         ne = ce->getCCWEdg();
         xy1 = ne->getRVtx();
         //checking polygon edge is on boundary and ccw edge's RVtx is on
         //wrong side of bndy edge...
         if( ce->getBoundaryFlag() && ne->getBoundaryFlag() )
         {
            //if( id==83 ) cout << " CASE A\n";
            bn0 = ce->getDestinationPtrNC();
            bn1 = ne->getDestinationPtrNC();
            xy2 = bn0->get2DCoords();
            xy3 = bn1->get2DCoords();
            if( !PointsCCW( xy1, xy2, xy3 ) )
            {
               //"cut off" portion of V. area outside bndy by finding intersections
               //of V. edges and bndy edge:
               //if( id==83 ) cout << " CASE B\n";
               xy = FindIntersectionCoords( ce->getRVtx(), xy1, xy2, xy3 );
               vcL.insertAtBack( xy );
               nne = ne->getCCWEdg();
               xy = FindIntersectionCoords( xy1, nne->getRVtx(), xy2, xy3 );
               vcL.insertAtBack( xy );
            }
            else vcL.insertAtBack( xy1 );
         }
         else vcL.insertAtBack( xy1 );
      }
      
      // Now that we've found the correct vertices, make triangles to
      // fill the polygon; the sum of the tri areas is the v. area.
      // For a convex polygon, we can compute the total area as the
      // sum of the area of triangles [P(1) P(i) P(i+1)] for i=2,3...N-1.
      //cout << "find polygon area" << endl << flush;
      // coords of first vertex:
      xy = *(vcI.FirstP()); //ce = vtxIter.FirstP();
      //if( id==83 ) cout << "starting pt " << xy[0] << "," << xy[1] <<endl;
      //xy = ce->getRVtx(); 
      // Find out # of vertices in polygon:
      int nverts = vcL.getSize(); //vedgList.getSize(); 
      for( i=2; i<=nverts-1; i++ )
      {
         xyn = *(vcI.NextP()); //xyn = vtxIter.NextP()->getRVtx();// Vertex i
         xynn = *(vcI.NextP());//vtxIter.ReportNextP()->getRVtx(); // Vertex i+1
         //if( id==83 ) cout << "other two: (" << xyn[0] << "," << xyn[1] << "), (" << 
         //    xynn[0] << "," << xynn[1] << ")\n";
         dx = xyn[0] - xy[0];
         dy = xyn[1] - xy[1];
         //if( id==83 ) cout << "dx: " << dx << " dy: " << dy << endl;
         a = sqrt( dx*dx + dy*dy );
         dx = xynn[0] - xyn[0];
         dy = xynn[1] - xyn[1];
         //if( id==83 ) cout << "dx: " << dx << " dy: " << dy << endl;
         b = sqrt( dx*dx + dy*dy );
         dx = xynn[0] - xy[0];
         dy = xynn[1] - xy[1];
         //if( id==83 ) cout << "dx: " << dx << " dy: " << dy << endl;
         c = sqrt( dx*dx + dy*dy );
         //TODO: check for sqrt neg # w/ assert
         area += 0.25*sqrt( 4*a*a*b*b -
                            (c*c - (b*b + a*a))*(c*c - (b*b + a*a)));
         /*if( id==83 ) {
            cout<<"ND "<<id<<" V_a 3: a "<<a<<", b "<<b<<", c "<<c<<endl<<flush;
            cout<<" Acum (1,"<<i<<","<<i+1<<") = " <<area<<endl;
            }*/
         vcI.Prev();
      }
   }
   varea = area;
   varea_rcp = 1.0/varea;

   // debug
   /*
   if( id==83 ) {
      cout << " reading list: ";
      for( ce = vtxIter.FirstP(); !(vtxIter.AtEnd()); ce = vtxIter.NextP() )
      {
         xy = ce->getRVtx();
         cout << xy[0] << " " << xy[1] << "; " << flush;
      }
      cout << endl << flush;
      cout << "reading spokes: ";
      ce = edg;
      do
      {
         assert( ce>0 );
         xy = ce->getRVtx();
         cout << xy[0] << " " << xy[1] << "; " << flush;
         ce = ce->getCCWEdg();
      } while( ce != edg );
      cout << endl << flush;
      }*/
   
   return area;
}


/*******************************************************************\
**
**  tNode::makeCCWEdges
**
**  This function provides for compatibility between the CCW Edge
**  data structure and the Spoke List data structure. It sets up
**  CCW edge connectivity from the spoke list data (which is 
**  assumed to be up to date) by: (1) setting the node's edg 
**  pointer to the first spoke on the list, and (2) setting the
**  ccwedg pointer for each spoke.
**
\*******************************************************************/
void tNode::makeCCWEdges()
{
   tEdge *ce, *ccwe;
   tPtrListIter< tEdge > spokIter( spokeList );
   
   ce = spokIter.FirstP();
   assert( ce != 0 );
   setEdg( ce );

//     if(id==793){
//        tNode * nbr = ce->getDestinationPtrNC();
//        cout<<"makeCCWEdges() node "<<id<<endl;
//        cout<<"edge "<<ce->getID()<<" dstn "<<nbr->getID()<<endl;
//     }
   for( ; !(spokIter.AtEnd()); ce = spokIter.NextP() )
   {
      ccwe = spokIter.ReportNextP();
      assert( ccwe != 0 );
      ce->setCCWEdg( ccwe );
//        if(id==793){
//           tNode * nbr = ccwe->getDestinationPtrNC();
//           cout<<"edge "<<ccwe->getID()<<" dstn "<<nbr->getID()<<endl;
//        }
   }
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
      if( ce->getBoundaryFlag()==kFlowAllowed )
      {
         ce->setFlowAllowed( 0 );
         // get complement and change it too TODO
      }
      
   } while( (ce=ce->getCCWEdg()) != edg );
   
}

/*******************************************************************\
**
**  tNode::WarnSpokeLeaving( tEdge * edglvingptr )
**
**  This function is called when an edge is being removed from the edge list.
**  If edg (the edge pointer member of tNode) is pointing to the edge
**  which will be removed, this edg must be updated.
**
**  edglvingptr is as it says, a pointer to the edge which will be
**  removed.
**
**  Called from tMesh::ExtricateEdge
**
**  9/98 NG and GT
\*******************************************************************/
void tNode::WarnSpokeLeaving( tEdge * edglvingptr )
{
   if( edglvingptr == edg )
       edg = edg->getCCWEdg();

}

/**********************************************************************\
 **
 **  tNode::InitializeNode()
 **
 **  A virtual function.
 **  This functions doesn't do anything here, only in inherited classes.
 **  Used for initializing things in newly created nodes that are set up
 **  for the rest of the nodes when the mesh is created.
 **
 **  1/1999  NG
 \**********************************************************************/
void tNode::InitializeNode()
{
}


/**************************************************************************\
**
**  TellAll
**
**  ...is a debugging routine that prints out a bunch of info about a node.
**
\**************************************************************************/
#ifndef NDEBUG
void tNode::TellAll()
{
   cout << " NODE " << id << ":\n";
   cout << "  x=" << x << " y=" << y << " z=" << z;
   cout << "  boundary: " << boundary
        << "\n  varea: " << varea << endl;
   if( edg )
       cout << "  points to edg #" << edg->getID() << endl;
   else cout << "  edg is undefined!\n";
   
}
#endif


/**************************************************************************\
\***  Functions for class tEdge  ******************************************/

/***********************************************************************\
**
**  Constructors & destructors:
**
**  Default:  initializes values to zero and makes rvtx a 2-elem array
**  Copy:  copies all values
**  Destructor:  no longer used
**
\***********************************************************************/

//default constructor
tEdge::tEdge()                                                       //tEdge
        : rvtx(2)
{
   id = 0;
   len = 0;
   slope = 0;
   vedglen = 0;
   org = dest = 0;
   ccwedg = 0;
   flowAllowed = 0;
     //cout << "tEdge()" << endl;
}

//copy constructor
tEdge::tEdge( const tEdge &original )                                //tEdge
{
   if( &original != 0 )
   {
      id = original.id;
      len = original.len;
      slope = original.slope;
      rvtx = original.rvtx;
      vedglen = original.vedglen;
      org = original.org;
      dest = original.dest;
      ccwedg = original.ccwedg;
      flowAllowed = original.flowAllowed;
   }
     //cout << "tEdge( orig )" << endl;
}

tEdge::~tEdge() {/*cout << "    ~tEdge()" << endl;*/}      //tEdge


/***********************************************************************\
**
**  Overloaded operators:
**
**    assignment: copies all values (spokelist's assignment operator
**                creates duplicate copy of list)
**    left shift: sends the following data to the output stream:
**                edge ID, length, slope, and origin and destination IDs
**
\***********************************************************************/
const tEdge &tEdge::operator=( const tEdge &original )               //tEdge
{
   if( &original != this )
   {
      id = original.id;
      len = original.len;
      slope = original.slope;
      rvtx = original.rvtx;
      vedglen = original.vedglen;
      org = original.org;
      dest = original.dest;
      ccwedg = original.ccwedg;
      flowAllowed = original.flowAllowed;
   }
   return *this;
}

//left shift
ostream &operator<<( ostream &output, const tEdge &edge )            //'tEdge'
{
   output << edge.id << " " << edge.len << " " << edge.slope << " " 
          << edge.org->getID()
          << " " << edge.dest->getID() << endl;
   return output;
}

//istream &operator>>( istream &input, tEdge &edge )                 //'tEdge'
//{
     //int oid, did;
     //cout << "edge id, origin id, dest id:";
     //input >> edge.id >> edge.org >> edge.dest; //temporarily assign id vals to ptrs
  // return input;
//}


/***********************************************************************\
**
**  tEdge "get" functions:
**
**  getID - returns ID #
**  getBoundaryFlag - returns boundary code
**  getLength - returns projectd length
**  getSlope - returns slope
**  getOriginPtr - returns const ptr to origin node
**  getDestinationPtr - returns const ptr to destination node
**  getOriginPtrNC - returns non-const ptr to origin node
**  getDestinationPtrNC - returns non-const ptr to destination node
**  getOrgZ- returns z value of origin node
**  getDestZ - returns z value of destination node
**  getCCWEdg - returns ptr to counterclockwise neighboring edge
**  FlowAllowed - returns the boundary flag, which indicates whether
**                or not the edge is an active flow conduit (which is
**                true as long as neither endpoint is a closed bdy node)
**  getRVtx - returns coordinates of right-hand Voronoi vertex as a
**            2-element array
**  getVEdgLen - returns the length of the corresponding Voronoi edge
**
\***********************************************************************/

int tEdge::getID() const {return id;}                                //tEdge

//return 0 if flow allowed to match kNonBoundary:
int tEdge::getBoundaryFlag() const                                   //tEdge
{return !( flowAllowed == kFlowAllowed );} 

double tEdge::getLength() const {return len;}                         //tEdge

double tEdge::getSlope() const {return slope;}                         //tEdge

const tNode *tEdge::getOriginPtr() const {return org;}                //tEdge

const tNode *tEdge::getDestinationPtr() const {return dest;}          //tEdge

tNode *tEdge::getOriginPtrNC() {return org;}                     //tEdge

tNode *tEdge::getDestinationPtrNC() {return dest;}               //tEdge

double tEdge::getOrgZ() 
{
   //Xconst tNode * org = getOriginPtr(); 5/99
   assert( org!=0 );
   return( org->getZ() );
}

double tEdge::getDestZ()
{
   //Xconst tNode * dest = getDestinationPtr(); 5/99
   assert( dest!=0 );
   return( dest->getZ() );
}

tEdge * tEdge::getCCWEdg() 
{
   return ccwedg;
}

int tEdge::FlowAllowed() 
{
   return flowAllowed;
}

tArray< double >
tEdge::getRVtx() const
{
   //tArray< double >  xy( rvtx );
   //return xy;
   //cout << "getRVtx: ";
   return rvtx;
}

double tEdge::getVEdgLen() const {return vedglen;}


/***********************************************************************\
**
**  tEdge "set" functions:
**
**  setID - sets ID # to val
**  setLength - sets length to val
**  setSlope - sets slope to slp
**  setOriginPtr - sets origin pointer to ptr (if ptr is nonzero)
**  setDestinationPtr - sets destination pointer to ptr (if nonzero)
**  setFlowAllowed - sets flowAllowed status to val
**  setCCWEdg - sets ptr to counter-clockwise neighbor to edg
**  setRVtx - sets the coordinates of the right-hand Voronoi vertex
**            (ie, the Voronoi vertex at the circumcenter of the RH
**            triangle) to the 1st two elements in arr, which is
**            assumed to be a 2-element array
**  setVEdgLen - sets vedglen to val (vedglen is the length of the
**               corresponding Voronoi cell edge)
**
**  Note: unless otherwise noted, no checking of range or validity is
**        performed in these routines (aside from assert statements)
**
\***********************************************************************/

void tEdge::setID( int val ) {
   assert( id>=0 );
   id = val;
   /*id = ( val >=0 ) ? val : 0;*/}           //tEdge

void tEdge::setLength( double val )                                   //tEdge
{
   assert( val>=0.0 );
   len = val;
   /*len = ( val >= 0.0 ) ? val : 0.0;*/}

void tEdge::setSlope( double slp )
{ slope = slp; }

void tEdge::setOriginPtr( tNode * ptr ) {if( ptr != 0 ) org = ptr;}  //tEdge

void tEdge::setDestinationPtr( tNode * ptr )                         //tEdge
{if( ptr != 0 ) dest = ptr;}

void tEdge::setFlowAllowed( int val )                                //tEdge
{
   assert( val==0 || val==1 );
   flowAllowed = val;
   /*flowAllowed = ( val == 0 || val == 1 ) ? val : 0;*/}

void tEdge::setCCWEdg( tEdge * edg )
{
   assert( edg > 0 );
     //assert( ccwedg > 0 );
   ccwedg = edg;
}

void tEdge::setRVtx( tArray< double > arr )
{
   assert( &arr != 0 );
   assert( arr.getSize() == 2 );
     //cout << "setRVtx for edge " << id
     //   << " to x, y, " << arr[0] << ", " << arr[1] << endl;
   rvtx = arr;
}

void tEdge::setVEdgLen( double val )
{
   assert( val>=0.0 );
   vedglen = val;
   /*vedglen = ( val > 0 ) ? val : 0;*/
}


/**************************************************************************\
**
**  tEdge::CalcLength
**
**  Computes the edge length and returns it. (Length is the projected
**  on the x,y plane). Assumes org and dest are valid.
**
\**************************************************************************/
double tEdge::CalcLength()
{
   //Xconst tNode * org = getOriginPtr(); 5/99
   //Xconst tNode * dest = getDestinationPtr(); 5/99
   assert( org!=0 );  // Failure = edge has no origin and/or destination node
   assert( dest!=0 );
   
   double dx = org->getX() - dest->getX();
   double dy = org->getY() - dest->getY();
   len = sqrt( dx*dx + dy*dy );
   return len;
}

/**************************************************************************\
**
**  tEdge::CalcSlope
**
**  Computes the slope of the edge as ( Zorg - Zdest ) / length.
**
**  Returns: the slope
**  Modifies: slope (data mbr)
**  Assumes: length >0; org and dest valid.
**
\**************************************************************************/
double tEdge::CalcSlope()
{
   //Xconst tNode * org = getOriginPtr(); 5/99
   //Xconst tNode * dest = getDestinationPtr(); 5/99
   
   assert( org!=0 );  // Failure = edge has no origin and/or destination node
   assert( dest!=0 );
   assert( len>0.0 );

   slope = ( org->getZ() - dest->getZ() ) / len;
   return slope;
}


/**************************************************************************\
**
**  tEdge::CalcVEdgLen
**
**  Calculates the length of the Voronoi cell edge associated with the
**  current triangle edge. The Voronoi cell edge length is equal to the
**  distance between the Voronoi vertex of the right-hand triangle and
**  the Voronoi vertex of the left-hand triangle. The vertex for the
**  right-hand triangle is stored with rvtx[] (and is assumed to be up to
**  date), and the vertex for the left-hand triangle is stored in the
**  edge's counter-clockwise (left-hand) neighbor (also assumed valid and
**  up to date).
**
**  Data mbrs modified:  vedglen
**  Returns:  the Voronoi edge length
**  Assumes:  ccwedg valid, rvtx[] up to date
**
\**************************************************************************/
double tEdge::CalcVEdgLen()
{
	assert( ccwedg!=0 );
	
	double dx, dy;
	
	dx = rvtx[0] - ccwedg->rvtx[0];
	dy = rvtx[1] - ccwedg->rvtx[1];
	vedglen = sqrt( dx*dx + dy*dy );
	return( vedglen );
}


/**************************************************************************\
**
**  tEdge::TellCoords
**
**  Debugging routine that reports edge ID and coords of endpoints.
**
\**************************************************************************/
void tEdge::TellCoords()
{
   cout << "EDGE " << id << ":\n";
   cout << "  " << org->getID() << " (" << org->getX() << ","
        << org->getY() << ") -> " << dest->getID() << " ("
        << dest->getX() << "," << dest->getY() << ")" << endl;
}


/**************************************************************************\
**
**  tEdge::WelcomeCCWNeighbor
**
**  Welcomes a new spoke to the neighborhood! neighbor is a new edge to
**  be inserted counter-clockwise from this one. We point neighbor at
**  the edge we're currently pointing to, and then point ourself to
**  neighbor, thus maintaining the edge connectivity.
**
**  Data mbrs modified:  ccwedg, neighbor->ccwedg
**  Created: 2/4/99 GT
**
\**************************************************************************/
void tEdge::WelcomeCCWNeighbor( tEdge * neighbor )
{
   assert( neighbor!=0 );
   assert( neighbor->org == org );
   neighbor->ccwedg = ccwedg;
   ccwedg = neighbor;
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
**  Functions for class tTriangle.
**
**  For a description of tTriangle objects, see meshElements.h.
**
**  Modifications:
**    - added tTriangle::FindCircumcenter() 1/11/98 gt
**
\**************************************************************************/

/***********************************************************************\
**
**  Constructors & destructors:
**
**  Default:  initializes node, edge, and triangle ptrs to zero.
**  Copy:  copies all values
**  Destructor:  no longer used
**
\***********************************************************************/

//default
tTriangle::tTriangle()                                           //tTriangle
{
   assert( p != 0 && e != 0 && t != 0 );
   for( int i=0; i<3; i++ )
   {
      p[i] = 0;
      e[i] = 0;
      t[i] = 0;
   }
     //cout << "tTriangle()" << endl;
}

//copy constructor
tTriangle::tTriangle( const tTriangle &init )                    //tTriangle
{
   assert( p != 0 && e != 0 && t != 0 );
   if( &init != 0 )
   {
      id = init.id;
      for( int i=0; i<3; i++ )
      {
         p[i] = init.p[i];
         e[i] = init.e[i];
         t[i] = init.t[i];
      }
   }
     //cout << "tTriangle( orig )" << endl;
}

//destructor
tTriangle::~tTriangle()                                          //tTriangle
{
     //cout << "    ~tTriangle()" << endl;
}


/***********************************************************************\
**
**  Overloaded operators:
**
**    assignment: copies all values 
**    left shift: sends the following data to the output stream:
**                triangle ID and the IDs of its 3 nodes, clockwise
**                edges, ad neighboring triangles (or -1 if no
**                neighboring triangle exists across a given face)
**    right shift: reads triangle ID and 3 other unspecified IDs from
**                 the input stream (the latter are not currently used
**                 for anything)
**
\***********************************************************************/

//overloaded assignment operator
const tTriangle &tTriangle::operator=( const tTriangle &init )   //tTriangle
{
   if( &init != this )
   {
      id = init.id;
      for( int i=0; i<3; i++ )
      {
         p[i] = init.p[i];
         e[i] = init.e[i];
         t[i] = init.t[i];
      }
   }
   return *this;
}

//left shift
ostream &operator<<( ostream &output, const tTriangle &tri )     //'tTriangle'
{
   int i;
   output << tri.id << ":";
   for( i=0; i<3; i++ )
       output << " " << tri.p[i]->getID();
   output << ";";
   for( i=0; i<3; i++ )
       output << " " << tri.e[i]->getID();
   output << ";";
   for( i=0; i<3; i++ )
   {
      if( tri.t[i] != 0 ) output << " " << tri.t[i]->getID();
      else  output << " -1";
   }
   output << endl;
   return output;
}

istream &operator>>( istream &input, tTriangle &tri )            //'tTriangle'
{
   int id1, id2, id3;
   cout << "triangle id, origin id, dest id:";
   input >> tri.id >> id1 >> id2 >> id3; //temporarily assign id vals to ptrs
     //tri.setPPtr( tMesh::h.getList().
   return input;
}


/***********************************************************************\
**
**  tTriangle "get" functions:
**
**  getID - returns ID #
**  pPtr - returns ptr to one of the 3 vertex nodes, as specified by
**         _index_ (index is 0, 1, or 2)
**  ePtr - returns ptr to one of the 3 clockwise edges, as specified by
**         _index_ (index is 0, 1, or 2)
**  tPtr - returns ptr to one of the 3 adjacent triangles, as specified
**         by _index_ (index is 0, 1, or 2)
**
\***********************************************************************/
int tTriangle::getID() const {return id;}                        //tTriangle

tNode *tTriangle::pPtr( int index )                              //tTriangle
{
   assert( index >= 0 && index <= 3 );
   return p[index];
}

tEdge *tTriangle::ePtr( int index )                              //tTriangle
{
   assert( index >= 0 && index <= 3 );
   return e[index];
}

tTriangle *tTriangle::tPtr( int index )                          //tTriangle
{
   assert( index >= 0 && index <= 3 );
   return t[index];
}


/***********************************************************************\
**
**  tTriangle "get" functions:
**
**  setID - sets ID #
**  setPPtr - sets pointer to one of the 3 vertex nodes, as specified
**            by _index_ (index is 0, 1, or 2)
**  setEPtr - sets pointer to one of the 3 clockwise edges, as specified
**            by _index_ (index is 0, 1, or 2)
**  setTPtr - sets pointer to one of the 3 adjacent triangles, as
**            specified by _index_ (index is 0, 1, or 2)
**
\***********************************************************************/
void tTriangle::setID( int val ) {id = ( val >= 0 ) ? val : 0;}  //tTriangle

void tTriangle::setPPtr( int index, tNode * ndptr )              //tTriangle
{
   assert( index >= 0 && index <= 3 );
   p[index] = ndptr;
}

void tTriangle::setEPtr( int index, tEdge * egptr )              //tTriangle
{
   assert( index >= 0 && index <= 3 );
   e[index] = egptr;
}

void tTriangle::setTPtr( int index, tTriangle * trptr )          //tTriangle
{
   assert( index >= 0 && index <= 3 );
   t[index] = trptr;
}


/**************************************************************************\
**
**  tTriangle::nVOp
**
**  Returns the side number (0, 1, or 2) of the neighboring triangle ct.
**  Assumes that ct _is_ one of the neighboring triangles.
**
** NOTE: for runtime error checking, may want to take out the assert
** and instead generate a runtime error when i>2.
\**************************************************************************/
int tTriangle::nVOp( tTriangle *ct )
{
   int i;

   for( i=0; i<4; i++ )
   {
      assert( i<3 );
      if( t[i] == ct ) return i;
   }
   return i;
}


/**************************************************************************\
**
**  tTriangle::nVtx
**
**  Returns the vertex number (0, 1, or 2) associated with node cn.
**  (In other words, it says whether cn is vertex 0, 1, or 2 in the 
**  triangle).
**  Assumes that cn _is_ one of the triangle's vertices.
**
** NOTE: for runtime error checking, may want to take out the assert
** and instead generate a runtime error when i>2.
\**************************************************************************/
int tTriangle::nVtx( tNode *cn )
{
   int i;
   for( i=0; i<4; i++ )
   {
      assert( i<3 );
      if( p[i] == cn ) return i;
   }
   return i;
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
tArray< double >
tTriangle::FindCircumcenter()
{
   double x1, y1, x2, y2, dx1, dy1, dx2, dy2, m1, m2;
   tArray< double > xyo, xyd1, xyd2, xy(2);

   assert( pPtr(0) && pPtr(1) && pPtr(2) );
   
   // Coordinates of triangle's nodes p0, p1, and p2
   xyo = pPtr(0)->get2DCoords();
   xyd1 = pPtr(1)->get2DCoords();
   xyd2 = pPtr(2)->get2DCoords();

   // Find the midpoints of the two sides (p0,p1) and (p0,p2) and store them 
   // in (x1,y1) & (x2,y2). Then get the distance between p0 and the 
   // midpoints of each side
   x1 = (xyo[0] + xyd1[0]) / 2;
   y1 = (xyo[1] + xyd1[1]) / 2;
   x2 = (xyo[0] + xyd2[0]) / 2;
   y2 = (xyo[1] + xyd2[1]) / 2;
   dx1 = x1-xyo[0];
   dy1 = y1-xyo[1];
   dx2 = x2-xyo[0];
   dy2 = y2-xyo[1];

   // Compute the intercept of the bisectors of the two sides:
   // Case: neither spoke is horizontal (ok to divide by dy1 and dy2)
   if( fabs(dy1)>0 && fabs(dy2)>0 )
   {
      assert( dy1!=0 && dy2!=0 );
      m1= -dx1/dy1;
      m2 = -dx2/dy2;
      assert( m1!=m2 ); // should never happen; means edges are parallel
      xy[0] = (y2 - m2 * x2 - y1 + m1 * x1) / (m1 - m2);
      xy[1] = m1 * (xy[0] - x1) + y1;
   }
   // Case: one spoke is horizontal, but neither are vertical
   else if( dx1!=0 && dx2!=0 )
   {
      assert( dx1!=0 && dx2!=0 );
      m1 = dy1/dx1;
      m2 = dy2/dx2;
      assert( m1!=m2 );
      xy[1] = (m1 * y1 + x1 - m2 * y2 - x2) / (m1 - m2);
      xy[0] = -xy[1] * m1 + m1 * y1 + x1;
   }
   // Special case: one is vertical, the other horizontal
   else
   {
      if( dx1!=0 )
      {
         xy[0] = xyo[0] + dx1;
         xy[1] = xyo[1] + dy2;
         assert( dx2==0 && dy1==0 );
      }
      else
      {
         xy[0] = xyo[0] + dx2;
         xy[1] = xyo[1] + dy1;
         assert( dx1==0 && dy2==0 );
      }
   }
   assert( &xy != 0 );

   return xy;
}


/* TellAll: debugging output routine */
#ifndef NDEBUG
void tTriangle::TellAll()
{
   int i;
   
   assert( this!=0 );
   cout << "TRIANGLE #" << id << ":\n";
   for( i=0; i<3; i++ )
   {
      cout << "  P" << i << " ";
      if( p[i]!=0 ) cout << p[i]->getID() << " (" << p[i]->getX() << ","
                         << p[i]->getY() << ")";
      else cout << "(ndef)";
      cout << "  E" << i << " ";
      if( e[i]!=0 ) cout << e[i]->getID();
      else cout << "(ndef)";
      cout << "  T" << i << " ";
      if( t[i]!=0 ) cout << t[i]->getID();
      else cout << "(ndef)";
      cout << endl;
   }
}
#endif

      
