/**************************************************************************\
**
**  gridElements.cpp: Functions for mesh element classes tNode, tEdge, and
**                    tTriangle.
**
**  Modifications:
**   - gt added code to recreate the "edg" implementation as an alternative
**     to spokelist. Dec 1 '97.
**   - previously separate tNode, tEdge, and tTriangle files grouped into
**     "gridElements", 1/20/98 gt
**
**  $Id: meshElements.cpp,v 1.2 1998-01-21 20:15:31 gtucker Exp $
\**************************************************************************/

#include <assert.h>
#include <math.h>
#include "gridElements.h"

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
//assignment
const tNode &tNode::operator=( const tNode &right )                  //tNode
{
   if( &right != this )
   {
      tPtrListIter< tEdge > spokIter;
      tEdge *ce;
      id = right.id;
      x = right.x;
      y = right.y;
      z = right.z;
      boundary = right.boundary;
      varea = right.varea;
      varea_rcp = right.varea_rcp;
      edg = right.edg;
      spokeList = right.spokeList;
        /*spokeList.Flush();
      if( &(right.spokeList) != 0 )
      {
         spokIter.Reset( right.getSpokeListNC() );
         for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
         {
            insertBackSpokeList( ce );
         }
      }*/
   }
   return *this;
}
//get
tArray< float >                                                       //tNode
tNode::get3DCoords() const
{
   tArray< float > xyz(3);
   xyz[0] = x;
   xyz[1] = y;
   xyz[2] = z;
   return xyz;
}

tArray< float >                                                   //tNode
tNode::get2DCoords() const
{
   tArray< float > xy(2);
   xy[0] = x;
   xy[1] = y;
   return xy;
}

int tNode::getID() const {return id;}                                //tNode
float tNode::getX() const {return x;}
float tNode::getY() const {return y;}
float tNode::getZ() const {return z;}
float tNode::getVArea() const {return varea;}                        //tNode
float tNode::getVArea_Rcp() const {return varea_rcp;}                //tNode
int tNode::getBoundaryFlag() const {return boundary;}                //tNode

tEdge * tNode::GetEdg() 
{
   return edg;
}

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

//set
void tNode::setID( int val ) {id = val;}                             //tNode
void tNode::setX( float val ) {x = val;}                             //tNode
void tNode::setY( float val ) {y = val;}                             //tNode
void tNode::setZ( float val ) {z = val;}                             //tNode
void tNode::setVArea( float val )                                    //tNode
{varea = ( val >= 0.0 ) ? val : 0.0;}

void tNode::setVArea_Rcp( float val )                                //tNode
{varea_rcp =  ( val >= 0.0 ) ? val : 0.0;}

void tNode::setBoundaryFlag( int val )                               //tNode
{boundary = (val >=0 && val <= 2) ? val : 0;}

void tNode::set2DCoords( float val1, float val2 )                    //tNode
{
   setX( val1 );
   setY( val2 );
}
void tNode::set3DCoords( float val1, float val2, float val3 )        //tNode
{
   setX( val1 );
   setY( val2 );
   setZ( val3 );
}

void tNode::SetEdg( tEdge * theEdg )
{
   assert( theEdg > 0 );
   edg = theEdg;
   cout << "Assigning edge " << theEdg->getID() << " to node " << getID() << endl;
}


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

void tNode::makeWheel() {spokeList.makeCircular();}                  //tNode

istream &operator>>( istream &input, tNode &node )                   //'tNode'
{
   cout << "x y z:" << endl;
   input >> node.x >> node.y >> node.z;
   return input;
}
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

const tEdge *tNode::NextSpoke( tPtrListNode< tEdge > * current )//tNode
{
   if( current == 0 || current == spokeList.getLast() )
       current = spokeList.getFirstNC();
   else current = current->getNextNC();
   return current->getPtr();
}

/*****************************************************************************\
**
**      Dist: replaces dis; the distance of the node from the line formed
**                by points p0 p1 
**                using x y 
**
**
**      Data members updated: 
**      Called by: 
**      Calls:  
**
\*****************************************************************************/
float tNode::Dist( tNode * p0, tNode * p1 )
{
  float a,b,c,res;

  a=(p1->y)-(p0->y);  
  b=-((p1->x)-(p0->x));
  c=-((a*(p0->x))+(b*(p0->y)));
  res=(a*x + b*y + c) / sqrt(a*a + b*b);
  if (res<0) res=-res;
  return(res);
}

void tNode::CalcSpokeVEdgLengths()
{
   tEdge *curedg, *nextedg;
   tArray< float > cvtx, nvtx;
   tPtrListIter< tEdge > spokIter( spokeList );
   
   for( curedg = spokIter.FirstP(); !( spokIter.AtEnd() );
        curedg = spokIter.NextP() )
   {
      nextedg = spokIter.ReportNextP();
      cvtx = curedg->getRVtx();
      nvtx = nextedg->getRVtx();
      float dx = nvtx[0] - cvtx[0];
      float dy = nvtx[1] - cvtx[1];
      nextedg->setVEdgLen( sqrt( dx * dx + dy * dy ) );
   }
}

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
**  ComputeVoronoiArea
**
**  Computes the node's Voronoi area by summing the area of embedded
**  triangles, and also calls CalcSpokeVEdgLengths to compute the length of
**  the sides of the Voronoi cell (the length of cell sides is needed if
**  for example transport between adjacent cells depends on the width of
**  their shared face).
**
**  The basic Voronoi polygon is described by the set of "right-hand
**  Voronoi vertices" associated with each spoke (edge). These vertices
**  are computed by SetVoronoiVertices() as the intersection of the
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
\*****************************************************************************/
float tNode::ComputeVoronoiArea()
{
     //calculates area by defining the voronoi polygon vertices, dividing that
     //polygon into triangles, and calc'ing the area of each triangle, sum of
     //tri areas is v. area.
//#if DEBUG
   //cout<<"VoronoiArea(...)...";
//#endif
   int cw;
   float area = 0;
   float a, b, c, dx, dy, dx0, dx1, dy0, dy1, dx2, dy2;
   float vx, vy, x0, y0, x1, y1, x2, y2, m1, m2;
   tEdge * ce, *edgptr;
   tPtrList< tEdge > vedgList/*( centre->getSpokeListNC() )*/;
   tPtrListIter< tEdge > //XspokIter( centre->getSpokeListNC() ),
       vtxIter( vedgList );
   tArray< float > xy, xyn, xynn, xynnn, xy1, xy2;
   int i;

/*X   for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.NextP() )
   {
      cout << "INSerting edg " << ce->getID() << endl;
      vedgList.insertAtBack( ce );
   }*/

   // Create a duplicate list of edges; we will modify this list to obtain
     //the correct vertices. In some cases, we may need to delete an edge
     //to get the correct list of vertices; we don't want to delete the
     //spoke ptr, so we make a duplicate list.
   ce = edg;
   do
   {
      assert( ce>0 );
      vedgList.insertAtBack( ce );
      ce = ce->GetCCWEdg();
   } while( ce != edg );
   vedgList.makeCircular();

   // Check boundary status: Voronoi area only defined for non-boundary nodes
   if( boundary == kNonBoundary )
   {
      CalcSpokeVEdgLengths();
      cw = TRUE;
      do
      {
           //go through the list; we want the vertex list to run CCW;
           //in some cases of long skinny triangles, the 'unimproved'
           //v. polygon sides may form loops; loops are detected by
           //finding two (2) consecutive 'CCW' vertices; i.e., where
           //the 'curvature' is CW rather than CCW. In such cases,
           //we delete one of the edges from the vertex list and find
           //the new vertex at the intersection of the perp. bisectors
           //of the edges to either 'side' of the deleted edge. Iterate.
           //Really. It works.
         cw = FALSE;
         for( ce = vtxIter.FirstP(); !( vtxIter.AtEnd() ); ce = vtxIter.NextP() )
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
                  assert( vtxIter.Prev() );
                  xy2 = vtxIter.DatPtr()->getDestinationPtr()->get2DCoords();
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
                     m1 = dy1/dx1;
                     m2 = dy2/dx2;
                     vy=(m1*y1+x1-m2*y2-x2)/(m1-m2);
                     vx= -vy*m1+m1*y1+x1;
                  }
                  assert( vtxIter.Prev() );
                  xyn[0] = vx;
                  xyn[1] = vy;
                  dx = xy[0] - vx;
                  dy = xy[1] - vy;
                  cout << "reset vedglen and rvtx for edge "
                       << vtxIter.DatPtr()->getID() << " to len "
                       << sqrt( dx*dx + dy*dy )
                       << ", x, y, " << xyn[0] << ", " << xyn[1] << endl;
                    //reset 'next' edge's vertex to newly found intersection,
                    //length adjusted accordingly
                  vtxIter.DatPtr()->setVEdgLen( sqrt( dx*dx + dy*dy ) );
                  vtxIter.DatPtr()->setRVtx( xyn );
                  cout << "reset vedglen and rvtx for edge "
                       << vtxIter.ReportNextP()->getID()
                       << " to len 0.0, x, y, " << xynnn[0] << ", "
                       << xynnn[1] << endl;
                    //reset 'next-next' edge's vertex to the coordinates
                    //of the 'next-next-next' edge's vertex; length to zero
                  vtxIter.ReportNextP()->setVEdgLen(0.0);
                  vtxIter.ReportNextP()->setRVtx( xynnn );
                    //delete the offending vertex's edge from list
                  assert( vedgList.removeNext( edgptr, vtxIter.NodePtr() ) );
               }
            }
            assert( vtxIter.Get( ce->getID() ) );
         }
      } while( cw ); //while we're still finding loops in the polygon

      // Now that we've found the correct vertices, make triangles to
      // fill the polygon; the sum of the tri areas is the v. area.
      // For a convex polygon, we can compute the total area as the
      // sum of the area of triangles [P(1) P(i) P(i+1)] for i=2,3...N-1.
      ce = vtxIter.FirstP();
      xy = ce->getRVtx(); // coords of first vertex
      int nverts = vedgList.getSize(); // Find out # of vertices in polygon
      for( i=2; i<=nverts-1; i++ )
      {
         xyn = vtxIter.NextP()->getRVtx();        // Vertex i
         xynn = vtxIter.ReportNextP()->getRVtx(); // Vertex i+1
         dx = xyn[0] - xy[0];
         dy = xyn[1] - xy[1];
         a = sqrt( dx*dx + dy*dy );
         dx = xynn[0] - xyn[0];
         dy = xynn[1] - xyn[1];
         b = sqrt( dx*dx + dy*dy );
         dx = xynn[0] - xy[0];
         dy = xynn[1] - xy[1];
         c = sqrt( dx*dx + dy*dy );
           //cout<<"V_a 3: a "<<a<<", b "<<b<<", c "<<c<<endl<<flush;
         area += 0.25*sqrt( 4*a*a*b*b -
                            (c*c - (b*b + a*a))*(c*c - (b*b + a*a)));
      }
   }
   
   cout<<"Varea: node id "<<id<<", varea "
       <<area<<endl;
     //cout<<"V_area done"<<endl<<flush;
   //cout << "finished: voronoi area = " << area << endl;
   varea = area;
   varea_rcp = 1.0/varea;
   return area;
}




/**************************************************************************\
**
**  Functions for class tEdge
**
\**************************************************************************/

/*#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "../tArray/tArray.h"
#include "../tPtrList/tPtrList.h"
#include "../tNode/tNode.h"
#include "tEdge.h"
#include <math.h>  // For sqrt() call in CalcLength
*/

//default constructor
tEdge::tEdge()                                                       //tEdge
{
   id = 0;
   len = 0;
   slope = 0;
   org = dest = 0;
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

tEdge::~tEdge() {/*cout << "    ~tEdge()" << endl;*/}                    //tEdge

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

int tEdge::getID() const {return id;}                                //tEdge

//return 0 if flow allowed to match kNonBoundary:
int tEdge::getBoundaryFlag() const                                   //tEdge
{return !( flowAllowed == kFlowAllowed );} 

float tEdge::getLength() const {return len;}                         //tEdge

float tEdge::getSlope() const {return slope;}                         //tEdge

const tNode *tEdge::getOriginPtr() const {return org;}                     //tEdge

const tNode *tEdge::getDestinationPtr() const {return dest;}               //tEdge

tNode *tEdge::getOriginPtrNC() {return org;}                     //tEdge

tNode *tEdge::getDestinationPtrNC() {return dest;}               //tEdge

float tEdge::getOrgZ() 
{
   const tNode * org = getOriginPtr();
   assert( org!=0 );
   return( org->getZ() );
}

float tEdge::getDestZ()
{
   const tNode * dest = getDestinationPtr();
   assert( dest!=0 );
   return( dest->getZ() );
}

tEdge * tEdge::GetCCWEdg() 
{
   return ccwedg;
}

int tEdge::FlowAllowed() 
{
   return flowAllowed;
}


void tEdge::setID( int val ) {id = ( val >=0 ) ? val : 0;}           //tEdge

void tEdge::setLength( float val )                                   //tEdge
{len = ( val >= 0.0 ) ? val : 0.0;}

void tEdge::setSlope( float slp )
{ slope = slp; }

void tEdge::setOriginPtr( tNode * ptr ) {if( ptr != 0 ) org = ptr;}  //tEdge

void tEdge::setDestinationPtr( tNode * ptr )                         //tEdge
{if( ptr != 0 ) dest = ptr;}

void tEdge::setFlowAllowed( int val )                                //tEdge
{flowAllowed = ( val == 0 || val == 1 ) ? val : 0;}

void tEdge::SetCCWEdg( tEdge * edg )
{
   assert( edg > 0 );
     //assert( ccwedg > 0 );
   ccwedg = edg;
}


/**************************************************************************\
**
**  CalcLength
**
**  Computes the edge length and returns it. (Length is the projected
**  on the x,y plane). Assumes org and dest are valid.
**
\**************************************************************************/
float tEdge::CalcLength()
{
   const tNode * org = getOriginPtr();
   const tNode * dest = getDestinationPtr();
   assert( org!=0 );  // Failure = edge has no origin and/or destination node
   assert( dest!=0 );
   
   float dx = org->getX() - dest->getX();
   float dy = org->getY() - dest->getY();
   len = sqrt( dx*dx + dy*dy );
   return len;
}


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

tArray< float >
tEdge::getRVtx() const
{
   tArray< float >  xy( rvtx );
   return xy;
}

float tEdge::getVEdgLen() const {return vedglen;}

void tEdge::setRVtx( tArray< float > arr )
{
   assert( &arr != 0 );
     //cout << "setRVtx for edge " << id
     //   << " to x, y, " << arr[0] << ", " << arr[1] << endl;
   rvtx = arr;
}

void tEdge::setVEdgLen( float val ) {vedglen = ( val > 0 ) ? val : 0;}



/**************************************************************************\
**
**  Functions for class tTriangle.
**
**  For a description of tTriangle objects, see gridElements.h.
**
**  Modifications:
**    - added tTriangle::FindCircumcenter() 1/11/98 gt
**
**
\**************************************************************************/

/*#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "../tArray/tArray.h"
#include "../tPtrListNode/tPtrListNode.h"
#include "../tPtrList/tPtrList.h"
#include "../tNode/tNode.h"
#include "../tEdge/tEdge.h"
#include "tTriangle.h"
*/

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

tTriangle::~tTriangle()                                          //tTriangle
{
     //cout << "    ~tTriangle()" << endl;
}

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
     //tri.setPPtr( tGrid::h.getList().
   return input;
}

/******************************************************************/
// NOTE: for runtime error checking, may want to take out the assert
// and instead generate a runtime error when i>2.
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


/*****************************************************************************\
**
**  FindCircumcenter
**
**  Finds the circumcenter of the triangle by finding the intersection of
**  the perpendicular bisectors of sides (p0,p1) and (p0,p2).
**
\*****************************************************************************/
tArray< float >
tTriangle::FindCircumcenter()
{
   float x, y, x1, y1, x2, y2, dx1, dy1, dx2, dy2, m1, m2;
   tArray< float > xyo, xyd1, xyd2, xy(2);

   assert( pPtr(0) && pPtr(1) && pPtr(2) );
   
   // Coordinates of triangle's nodes p0, p1, and p2
   xyo = pPtr(0)->get2DCoords();
   xyd1 = pPtr(1)->get2DCoords();
   xyd2 = pPtr(2)->get2DCoords();

      // Find the midpoints of the two sides (p0,p1) and (p0,p2) and store in 
      // (x1,y1) & (x2,y2), then get the distance between p0 and the 
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
