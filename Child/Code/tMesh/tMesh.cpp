/***************************************************************************\
**
**  tGrid.cpp: Functions for class tGrid
**
**  $Id: tMesh.cpp,v 1.6 1998-01-29 19:59:03 gtucker Exp $
\***************************************************************************/

#include <assert.h>
#include <math.h>
#include "../tLNode/tLNode.h"
#include "../Mathutil/mathutil.h"
#include "tGrid.h"


     //global functions:


//global function; determines whether test point violates "Delaunay-ness"
//of other three; i.e., does the triangle 'pass the test' against the other
int TriPasses( tArray< double > &ptest,
               tArray< double > &p0,
               tArray< double > &p1,
               tArray< double > &p2 )
{
   assert( (&ptest != 0) && (&p0 != 0) && (&p1 != 0) && (&p1 != 0) );
     //cout << "TriPasses? ";
   double dx0, dx1, dy0, dy1;
   double crossp, dotp, angle0_2_1, angle0_test_1;
   dx0 = p0[0] - p2[0];
   dx1 = p1[0] - p2[0];
   dy0 = p0[1] - p2[1];
   dy1 = p1[1] - p2[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   dotp = dx0 * dx1 + dy0 * dy1;
   angle0_2_1 = atan2( crossp, dotp );
   dx0 = p0[0] - ptest[0];
   dx1 = p1[0] - ptest[0];
   dy0 = p0[1] - ptest[1];
   dy1 = p1[1] - ptest[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   dotp = dx0 * dx1 + dy0 * dy1;
   angle0_test_1 = atan2( crossp, dotp );
   //determine:
   if( angle0_2_1 < angle0_test_1 )
   {
        //cout << "Yes" << endl;
      return 0;
   }
   else
   {
        //cout << "No" << endl;
      return 1;
   }
}

//global function; determines whether points are counter-clockwise:
int PointsCCW( tArray< double > &p0,
               tArray< double > &p1,
               tArray< double > &p2 )
{
   assert( &p0 != 0 && &p1 != 0 && &p1 != 0 );
     //cout << "PointsCCW? ";
   int i;
   double dx0, dx1, dy0, dy1;
   double crossp;
   dx0 = p1[0] - p0[0];
   dx1 = p2[0] - p0[0];
   dy0 = p1[1] - p0[1];
   dy1 = p2[1] - p0[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   if( crossp > 0 )
   {
        //cout << "Yes" << endl;
      return 1;
   }
   else
   {
        //cout << "No" << endl;
      return 0;
   }
}

//global function; determines whether triangle's "new" coords are CCW
int NewTriCCW( tTriangle *ct )
{
   assert( ct != 0 );
   tLNode *cn;
   cn = (tLNode *) ct->pPtr(0);
   tArray< double > p0( cn->get2DCoords() );
   if( cn->Meanders() ) p0 = cn->getNew2DCoords();
   cn = (tLNode *) ct->pPtr(1);
   tArray< double > p1( cn->get2DCoords() );
   if( cn->Meanders() ) p1 = cn->getNew2DCoords();
   cn = (tLNode *) ct->pPtr(2);
   tArray< double > p2( cn->get2DCoords() );
   if( cn->Meanders() ) p2 = cn->getNew2DCoords();
   if( PointsCCW( p0, p1, p2 ) ) return 1;
   else return 0;
}

//global function; determines whether coords are in "new" triangle
int InNewTri( tArray< double > &xy, tTriangle *ct )
{
   tLNode *vtx;
   tArray< double > xy1, xy2;
   for( int j=0; j<3; j++ )
   {
      vtx = (tLNode *) ct->pPtr(j);
      if( vtx->Meanders() ) xy1 = vtx->getNew2DCoords();
      else xy1 = vtx->get2DCoords();
      vtx = (tLNode *) ct->pPtr( (j+1)%3 );
      if( vtx->Meanders() ) xy2 = vtx->getNew2DCoords();
      else xy2 = vtx->get2DCoords();
      if ( ( (xy1[1] - xy[1]) * (xy2[0] - xy[0]) ) >
           ( (xy1[0] - xy[0]) * (xy2[1] - xy[1])) )
          break;
   }
   if( j == 3) return 1;
   else return 0;
}


//global function; determines whether nbr node currently pointed to
//by iterator and the next two in the nbr list form a Delaunay triangle:
template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &nbrList,
                   tPtrListIter< tSubNode > &nbrIter )
{
   static ncalls = 0;
   ncalls++;
   
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   tSubNode *cn, *nbrnd;
   nbrnd = nbrIter.DatPtr();
     //cout << "Next3Delaunay? no calls =  " << ncalls << endl;
     //global function prototypes
     //int PointsCCW( tArray< double > &, tArray< double > &,
     //           tArray< double > & );
     //int TriPasses( tArray< double > &, tArray< double > &,
     //           tArray< double > &, tArray< double > & );
     //cout << "N3D: pt a\n";
   tPtrListIter< tSubNode > nbrIterCopy( nbrList );
     //cout << "N3D: pt b\n";
   int i = nbrIter.Where();
     //cout << "N3D: pt c\n";
   nbrIterCopy.Get(i);
     //cout << "N3D: point d\n";
   tPtrListNode< tSubNode > *tempptrnode = nbrIter.NodePtr();
   tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   tArray< double > p1( nbrIterCopy.NextP()->get2DCoords() );
   tArray< double > p2( nbrIterCopy.NextP()->get2DCoords() );
     //cout << "N3D: point B\n";
   if( !PointsCCW( p0, p1, p2 ) ) return 0;
   tArray< double > ptest;
   cn = nbrIterCopy.NextP();
   while( cn != nbrnd )
   {
      ptest = cn->get2DCoords();
      if( !TriPasses( ptest, p0, p1, p2 ) )
      {
           //cout << "Next3Delaunay? No" << endl;
         return 0;
      }
        //else cout << "Next3Del? this tri passed..\n";
      
      cn = nbrIterCopy.NextP();
   }
   //cout << "Next3Delaunay? Yes" << endl;
   return 1;
}

//global function; determines whether nbr node currently pointed to
//by iterator and the next two in the nbr list form a Delaunay triangle:
template< class tSubNode >
int PointAndNext2Delaunay( tSubNode &testNode, tPtrList< tSubNode > &nbrList,
                           tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) && (&testNode != 0) );
   //cout << "PointAndNext2Delaunay?" << endl;
     //global function prototypes
     //int PointsCCW( tArray< double > &, tArray< double > &,
     //           tArray< double > & );
     //int TriPasses( tArray< double > &, tArray< double > &,
     //           tArray< double > &, tArray< double > & );
   tPtrListIter< tSubNode > nbrIterCopy( nbrList );
   int i = nbrIter.Where();
     //cout << "Where: " << i << endl;
   nbrIterCopy.Get( i );
     //assert( nbrIterCopy.Get( i ) );
   assert( nbrIterCopy.DatPtr() == nbrIter.DatPtr() );
   tPtrListNode< tSubNode > *tempptrnode = nbrIter.NodePtr();
   tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   assert( nbrIterCopy.Next() );
   tArray< double > p1( nbrIterCopy.DatPtr()->get2DCoords() );
   tArray< double > p2( testNode.get2DCoords() );
   if( !PointsCCW( p0, p1, p2 ) ) return 0;
   tArray< double > ptest;
   assert( nbrIterCopy.Next() );
   while( nbrIterCopy.DatPtr() != nbrIter.DatPtr() )
   {
      ptest = nbrIterCopy.DatPtr()->get2DCoords();
      if( !TriPasses( ptest, p0, p1, p2 ) )
      {
         //cout << "No" << endl;
         return 0;
      }
      assert( nbrIterCopy.Next() );
   }
   //cout << "Yes" << endl;
   return 1;
}

/*****************************************************************************\
**
**      Intersect: uses newx, newy
**      Data members updated: Grid
**      Called by: 
**      Calls:  
**
\*****************************************************************************/
int Intersect( tEdge * ae, tEdge * be )
{
   //cout << "Intersect(...)..." << endl;
   double dxa, dxb, dya, dyb,
       a, b, c, f, g, h, x, y,
       rangemina, rangeminb, rangemaxa, rangemaxb, rangemin, rangemax;
   double smallnum = 0.0000000001;
   tLNode * lnode;
   if( !ae || !be )
   {
      cout<<"Intersect: Warning: invalid edge(s)"<<endl<<flush;
      return( NULL );
   }
   if( !ae->getOriginPtr() || !ae->getDestinationPtr() ||
       !be->getOriginPtr() || !be->getOriginPtr() )
   {
      cout<<"Intersect: Warning: invalid org or dest"<<endl<<flush;
      return( NULL );
   }
   lnode = (tLNode *) ae->getOriginPtrNC();
   tArray< double > aoc( lnode->get2DCoords() );
   if( lnode->Meanders() ) aoc = lnode->getNew2DCoords();
   lnode = (tLNode *) ae->getDestinationPtrNC();
   tArray< double > adc( lnode->get2DCoords() );
   if( lnode->Meanders() ) adc = lnode->getNew2DCoords();
   lnode = (tLNode *) be->getOriginPtrNC();
   tArray< double > boc( lnode->get2DCoords() );
   if( lnode->Meanders() ) boc = lnode->getNew2DCoords();
   lnode = (tLNode *) be->getDestinationPtrNC();
   tArray< double > bdc( lnode->get2DCoords() );
   if( lnode->Meanders() ) bdc = lnode->getNew2DCoords();
   
   dxa = adc[0] - aoc[0] + smallnum;
   dxb = bdc[0] - boc[0] + smallnum;
   dya = adc[1] - aoc[1] + smallnum;
   dyb = bdc[1] - boc[1] + smallnum;
   a = dya;
   b = -dxa;
   c = dxa * aoc[1] - dya * aoc[0];
   f = dyb;
   g = -dxb;
   h = dxb * boc[1] - dyb * boc[0];
   if( fabs(dxa) > 0 && fabs(dxb) > 0 )
   {
      if( fabs(f - g * a / b) > smallnum )
      {
         x = (g * c / b - h) / (f - g * a / b);
         y = (-c - a * x) / b;
         if( adc[0] < aoc[0] )
         {
            rangemina = adc[0];
            rangemaxa = aoc[0];
         }
         else
         {
            rangemina = aoc[0];
            rangemaxa = adc[0];
         }
         if( bdc[0] < boc[0] )
         {
            rangeminb = bdc[0];
            rangemaxb = boc[0];
         }
         else
         {
            rangeminb = boc[0];
            rangemaxb = bdc[0];
         }
         if( rangemina < rangeminb ) rangemin = rangeminb + smallnum;
         else rangemin = rangemina + smallnum;
         if( rangemaxa > rangemaxb ) rangemax = rangemaxb - smallnum;
         else rangemax = rangemaxa - smallnum;
         if( x > rangemin && x < rangemax )
             return( 1 );
         else
             return( 0 );
      }
      else return (0);
   }
   else
   {
      if( fabs(g - f * b / a) > smallnum )
      {
         y = (f * c / a - h) / (g - f * b / a);
         x = (-c - b * y) / a;
         if( adc[1] < aoc[1] )
         {
            rangemina = adc[1];
            rangemaxa = aoc[1];
         }
         else
         {
            rangemina = aoc[1];
            rangemaxa = adc[1];
         }
         if( bdc[1] < boc[1] )
         {
            rangeminb = bdc[1];
            rangemaxb = boc[1];
         }
         else
         {
            rangeminb = boc[1];
            rangemaxb = bdc[1];
         }
         if( rangemina < rangeminb ) rangemin = rangeminb + smallnum;
         else rangemin = rangemina + smallnum;
         if( rangemaxa > rangemaxb ) rangemax = rangemaxb - smallnum;
         else rangemax = rangemaxa - smallnum;
         if( y > rangemin && y < rangemax ) return( 1 );
         else return( 0 );
      }
      else return 0;
   }
}

/**************************************************************************\
**
**         Utilities for class tGrid
**
\**************************************************************************/
//default constructor
template< class tSubNode >
tGrid< tSubNode >::
tGrid() {nnodes = nedges = ntri = seed = 0;cout<<"tGrid()"<<endl;}     //tGrid

template< class tSubNode >
tGrid< tSubNode >::
tGrid( tListInputData< tSubNode > &input )                           //tGrid
{
   seed = 0;
   // Assign number of nodes, edges, and triangles
   nnodes = input.x.getSize();
   nedges = input.orgid.getSize();
   ntri = input.p0.getSize();
   assert( nnodes > 0 );
   assert( nedges > 0 );
   assert( ntri > 0 );

   // Create the node list by creating a temporary node and then iteratively
   // (1) assigning it values from the input data and (2) inserting it onto
   // the back of the node list.
   // TODO: insertion position should depend on boundary flag
   tSubNode tempnode;
   for( int i = 0; i< nnodes; i++ )
   {
      tempnode.set3DCoords( input.x[i], input.y[i], input.z[i] );
      tempnode.setID( i );
      tempnode.setBoundaryFlag( input.boundflag[i] );
      nodeList.insertAtBack( tempnode );
   }

   // Create and initialize the edge list by creating two temporary edges
   // (which are complementary, ie share the same endpoints) and then
   // iteratively assigning values to the pair and inserting them onto the
   // back of the edgeList
   tGridListIter< tSubNode > nodIter( nodeList );
   tEdge tempedge1, tempedge2;
   for( i = 0; i< nedges-1; i+=2 )
   {
      tempedge1.setID( i );
      tempedge2.setID( i );
      if( nodIter.Get( input.orgid[i] ) )
      {
         tempedge1.setOriginPtr( &(nodIter.DatRef()) );
         tempedge2.setOriginPtr( &(nodIter.DatRef()) );
      }
      if( nodIter.Get( input.destid[i] ) )
      {
         tempedge1.setDestinationPtr( &(nodIter.DatRef()) );
         tempedge2.setDestinationPtr( &(nodIter.DatRef()) );
      }
      edgeList.insertAtBack( tempedge1 );
      edgeList.insertAtBack( tempedge2 );
   }

   // Set up the lists of edges (spokes) connected to each node
   // (GT added code to also assign the 1st edge to "edg" as an alternative
   // to spokelist implementation)
   int e1;
   int ne;
   tSubNode * curnode;
   tGridListIter< tEdge > edgIter( edgeList );
   i = 0;
   do                                        //for( i=0; i<nnodes; i++ )
   {
      e1 = input.edgid[i]; //THIS REF TO "i" IS UNDEFINED! fixed
      curnode = nodIter->DatPtr();
      if( edgIter.Get( e1 ) ) // TODO: report error if fails
      {
          curnode->insertBackSpokeList( &(edgIter.DatRef()) );
          //nodIter.DatRef().insertBackSpokeList( &(edgIter.DatRef()) );
          curnode->SetEdg( edgIter.DatPtr() );
      }
      for( ne = input.nextid[e1]; ne != e1; ne = input.nextid[ne] )
          if( edgIter.Get( ne ) ) // TODO: report error if fails
              //nodIter.DatRef().insertBackSpokeList( &(edgIter.DatRef()) );
              curnode->insertBackSpokeList( &(edgIter.DatRef()) );
      i++;
      
   } while( nodIter.Next() );

   // Assign ccwedg connectivity (that is, tell each edge about its neighbor
   // immediately counterclockwise) (added by gt Dec 97 to re-implement
   // previous data structure)
   tEdge * curedg, * ccwedg;
   for( i=0; i<nedges; i++ )
   {
      curedg = edgIter.Get( i );
      ccwedgid = input.nextid[i];
      ccwedg = edgIter.Get( ccwedgid );
      curedg->SetCCWEdg( ccwedg );
      cout << "Edg " << ccwedgid << " (" << ccwedg->getOriginPtr()->getID()
           << " " << ccwedg->getDestinationPtr()->getID() << ") is ccw from "
           << curedg->getID() << " ("
           << curedg->getOriginPtr()->getID()
           << " " << curedg->getDestinationPtr()->getID() << ") " << endl;
      
   }
   
   tTriangle temptri;
   for ( i=0; i<ntri; i++ )
   {
      temptri.setID( i );
      if( nodIter.Get( input.p0[i] ) )
          temptri.setPPtr( 0, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p1[i] ) )
          temptri.setPPtr( 1, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p2[i] ) )
          temptri.setPPtr( 2, &(nodIter.DatRef()) );
      if( edgIter.Get( input.e0[i] ) )
          temptri.setEPtr( 0, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e1[i] ) )
          temptri.setEPtr( 1, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e2[i] ) )
          temptri.setEPtr( 2, &(edgIter.DatRef()) );
      triList.insertAtBack( temptri );
   }
   
   tListIter< tTriangle > trIter1( triList );
   tListIter< tTriangle > trIter2( triList );
   cout << "tGrid: set tri ptrs" << endl;
   i = 0;
   do                                                //for ( i=0; i<ntri; i++ )
   {
      cout << "tGrid: Tri loop, i: " << i << endl;
      if( trIter2.Get( input.t0[i] ) )
          trIter1.DatRef().setTPtr( 0, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 0, 0 );
      if( trIter2.Get( input.t1[i] ) )
          trIter1.DatRef().setTPtr( 1, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 1, 0 );
      if( trIter2.Get( input.t2[i] ) )
          trIter1.DatRef().setTPtr( 2, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 2, 0 );
      i++;
   }
   while( trIter1.Next() );
   cout<<"set list data values"<<endl;

   CheckMeshConsistency();
   
}

//constructor specifically for tLNodes
tGrid< tLNode >::
tGrid( tListInputData< tLNode > &input, int lnodflag = 0 )                                         //tGrid
{
   seed = 0;
   // Set the number of nodes, edges, and triangles in the grid mesh
   assert( lnodflag );
   nnodes = input.x.getSize();
   nedges = input.orgid.getSize();
   ntri = input.p0.getSize();
   cout << "nnodes, nedges, ntri: " << nnodes << " " << nedges << " " << ntri << endl;
   assert( nnodes > 0 );
   assert( nedges > 0 );
   assert( ntri > 0 );

   // Create the node list by creating a temporary node and then iteratively
   // (1) assigning it values from the input data and (2) inserting it onto
   // the back of the node list.
   tLNode tempnode;
   int bound;
   for( int i = 0; i< nnodes; i++ )
   {
      tempnode.set3DCoords( input.x[i], input.y[i], input.z[i] );
      tempnode.setID( i );
      bound = input.boundflag[i];
      assert( bound >= 0 && bound <= 2 );
      tempnode.setBoundaryFlag( bound );
      if( bound == kNonBoundary )
          nodeList.insertAtActiveBack( tempnode );
      else if( bound == kOpenBoundary )
          nodeList.insertAtBoundFront( tempnode );
      else
          nodeList.insertAtBack( tempnode );       //kClosedBoundary
      //cout << input.x[i] << input.y[i] << input.z[i]
      //     << input.boundflag[i] << endl;
      cout << tempnode.getBoundaryFlag() << " ";
      cout << nodeList.getLast()->getDataPtr()->getBoundaryFlag() << endl;
   }
   
   // Create and initialize the edge list by creating two temporary edges
   // (which are complementary, ie share the same endpoints) and then
   // iteratively assigning values to the pair and inserting them onto the
   // back of the edgeList
   tGridListIter< tLNode > nodIter( nodeList );
   tEdge tempedge1, tempedge2;
   int obnd, dbnd;
   for( i = 0; i< nedges-1; i+=2 )
   {
      // Assign values: ID, origin and destination pointers
      tempedge1.setID( i );
      tempedge2.setID( i + 1 );
      cout << input.orgid[i] << " " << input.destid[i] << endl;
      cout << nodIter.Get( input.orgid[i] ) << " ";
      cout << nodIter.Get( input.destid[i] ) << endl;
      assert( nodIter.Get( input.orgid[i] ) );
          //{
         tempedge1.setOriginPtr( &(nodIter.DatRef()) );
         tempedge2.setDestinationPtr( &(nodIter.DatRef()) );
         obnd = nodIter.DatRef().getBoundaryFlag();
         cout << nodIter.DatRef().getID() << "->";
         //}
         assert( nodIter.Get( input.destid[i] ) );
          //{
         tempedge1.setDestinationPtr( &(nodIter.DatRef()) );
         tempedge2.setOriginPtr( &(nodIter.DatRef()) );
         dbnd = nodIter.DatRef().getBoundaryFlag();
         cout << nodIter.DatRef().getID() << endl;
         //}
      // Set the "flowallowed" status (FALSE if either endpoint is a
         // closed boundary) and insert edge pair onto the list --- active
         // part of list if flow is allowed, inactive if not
         cout << "BND: " << obnd << " " << dbnd << " " << kClosedBoundary
              << endl;
      if( obnd == kClosedBoundary || dbnd == kClosedBoundary )
      {
         cout << "setting edges " << tempedge1.getID() << " and "
              << tempedge2.getID() << " as no-flux" << endl;
         tempedge1.setFlowAllowed( 0 );
         tempedge2.setFlowAllowed( 0 );
         edgeList.insertAtBack( tempedge1 );
         edgeList.insertAtBack( tempedge2 );
      }
      else
      {
         cout << "setting edges " << tempedge1.getID() << " and "
              << tempedge2.getID() << " as OPEN" << endl;
         tempedge1.setFlowAllowed( 1 );
         tempedge2.setFlowAllowed( 1 );
         edgeList.insertAtActiveBack( tempedge1 );
         edgeList.insertAtActiveBack( tempedge2 );
         cout << "EDGFA " << tempedge2.FlowAllowed() << endl;
      }
   }

   // Set up the lists of edges (spokes) connected to each node
   // (GT added code to also assign the 1st edge to "edg" as an alternative
   // to spokelist implementation)
   int e1;
   int ne;
   tGridListIter< tEdge >
       edgIter( edgeList );
   tLNode * curnode;
   assert( nodIter.First() );
   i = 0;
   do                                        //for( i=0; i<nnodes; i++ )
   {
      //e1 = input.edgid[i]; //error: assumes nodes are in ID order!
      curnode = nodIter.DatPtr();
      e1 = input.edgid[curnode->getID()];  //fix of above error
      if( edgIter.Get( e1 ) )
      {
          curnode->insertBackSpokeList( &(edgIter.DatRef()) );
          //nodIter.DatRef().insertBackSpokeList( &(edgIter.DatRef()) );
          curnode->SetEdg( edgIter.DatPtr() );
          cout << "Node " << curnode->getID() << " has edg "
               << (curnode->GetEdg())->getID() << endl;
      }
      for( ne = input.nextid[e1]; ne != e1; ne = input.nextid[ne] )
      {
         if( ne>=nedges )
         {
             cerr << "Warning: edge " << e1 << " has non-existant ccw edge "
                  << ne << endl;
             cerr << "This is likely to be a problem in the edge input file"
                  << endl;
         }
         if( edgIter.Get( ne ) )
              curnode->insertBackSpokeList( &(edgIter.DatRef()) );
          //nodIter.DatRef().insertBackSpokeList( &(edgIter.DatRef()) );
      }
      i++;
   }
   while( nodIter.Next() );

   // Assign ccwedg connectivity (that is, tell each edge about its neighbor
   // immediately counterclockwise) (added by gt Dec 97 to re-implement
   // previous data structure)
   tEdge * curedg, * ccwedg;
   int ccwedgid;
   for( i=0; i<nedges; i++ )
   {
      curedg = edgIter.GetP( i );
      ccwedgid = input.nextid[i];
      ccwedg = edgIter.GetP( ccwedgid );
      curedg->SetCCWEdg( ccwedg );
      cout << "Edg " << ccwedgid << " (" << ccwedg->getOriginPtr()->getID()
           << " " << ccwedg->getDestinationPtr()->getID() << ") is ccw from "
           << curedg->getID() << " ("
           << curedg->getOriginPtr()->getID()
           << " " << curedg->getDestinationPtr()->getID() << ") " << endl;
      
   }
   
   tPtrListIter< tEdge > spokIter;
   assert( nodIter.First() );
   do                                        //for( i=0; i<nnodes; i++ )
   {
      spokIter.Reset( nodIter.DatRef().getSpokeListNC() );
      cout << " node " << nodIter.DatRef().getID() << " with spoke edges";
      i = 0;
      do
      {
         if( i > 0 ) spokIter.Next();
         cout << " " << spokIter.DatPtr()->getID();
         i++;
      }
      while( spokIter.NextIsNotFirst() );
      cout << endl;
   }
   while( nodIter.Next() );
   

   tTriangle temptri;
   for ( i=0; i<ntri; i++ )
   {
      temptri.setID( i );
      if( nodIter.Get( input.p0[i] ) )
          temptri.setPPtr( 0, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p1[i] ) )
          temptri.setPPtr( 1, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p2[i] ) )
          temptri.setPPtr( 2, &(nodIter.DatRef()) );
      if( edgIter.Get( input.e0[i] ) )
          temptri.setEPtr( 0, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e1[i] ) )
          temptri.setEPtr( 1, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e2[i] ) )
          temptri.setEPtr( 2, &(edgIter.DatRef()) );
      triList.insertAtBack( temptri );
   }
   
   tListIter< tTriangle >
       trIter1( triList );
   tListIter< tTriangle >
       trIter2( triList );
   cout << "tGrid: set tri ptrs" << endl;
   i = 0;
   do                                                //for ( i=0; i<ntri; i++ )
   {
      cout << "tGrid: Tri loop, i: " << i << endl;
      if( trIter2.Get( input.t0[i] ) )
          trIter1.DatRef().setTPtr( 0, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 0, 0 );
      if( trIter2.Get( input.t1[i] ) )
          trIter1.DatRef().setTPtr( 1, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 1, 0 );
      if( trIter2.Get( input.t2[i] ) )
          trIter1.DatRef().setTPtr( 2, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 2, 0 );
      i++;
   }
   while( trIter1.Next() );
   cout<<"set list data values"<<endl;
   cout << "tGrid( input )" << endl;

   CheckMeshConsistency();
}
/*
**   tGrid constructor to make grid from scratch; reads parameters
**        from input file to get grid size, spacing, method of point
**        placement
*/
template< class tSubNode >
tGrid< tSubNode >::
tGrid( tInputFile &infile )
{
   int i, j, id, n, nx, ny;
   int numPts;
   double dist;
   double delGrid;
   tArray< double > xyz(3);
   tSubNode tempnode, *cn, *node0, *node1, *node2, *node3;
   tEdge *ce;
   tGridListIter< tEdge > edgIter( edgeList );
   tGridListIter< tSubNode > nodIter( nodeList );
   tPtrList< tSubNode > bndList;
   cout << "THIS IS TEST CONTRUCTOR FOR tGrid\n";
   seed = infile.ReadItem( seed, "SEED" );
     //reads in size of grid (meters)
   double xGrid = infile.ReadItem( xGrid, "X_GRID_SIZE" );
   double yGrid = infile.ReadItem( yGrid, "Y_GRID_SIZE" );
     //read type of open boundary:
     //  0 = corner outlet (lower left)
     //  1 = one open side (lower)
     //  2 = two opposite sides (upper and lower)
     //  3 = all sides
     //  4 = specify outlet coordinates
   int boundType = infile.ReadItem( boundType, "TYP_BOUND" );
     //read mean elevation (also determines elev variation)
   double mElev = infile.ReadItem( mElev, "MEAN_ELEV" );
     //reads method of point placement:
     //  0 = uniform grid;
     //  1 = perturbed grid;
     //  2 = random scatter;
   int ptPlace = infile.ReadItem( ptPlace, "OPT_PT_PLACE" );
     //read grid spacing or number of points (excluding four boundary points)
   if( ptPlace == kUniformGrid || ptPlace == kPerturbedGrid )
   {
      delGrid = infile.ReadItem( delGrid, "GRID_SPACING" );
   }
   else
   {
      numPts = infile.ReadItem( numPts, "NUM_PTS" );
      delGrid = sqrt( xGrid * yGrid / numPts );
   }

   //MAKE BOUNDARY
   if( boundType == kCornerOutlet )
   {
      id = 0;
      tempnode.setBoundaryFlag( kOpenBoundary );
      tempnode.set3DCoords( 0, 0, 0 );
      tempnode.setID( id );
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      nodeList.insertAtBack( tempnode );
      bndList.insertAtBack( nodIter.LastP() );
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=1, id++; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( boundType == kOpenSide )
   {
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=1, id=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( i - 1 );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=n; i>=0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   if( boundType == kOppositeSidesOpen )
   {
      double upperZ = infile.ReadItem( upperZ, "UPPER_BOUND_Z" );
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=1, id=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = yGrid / delGrid;
      for( i=0; i<=n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=n-1; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, upperZ );
         tempnode.setID( id );
         nodeList.insertAtBoundFront( tempnode );
         bndList.insertAtBack( nodIter.GetP( id ) );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = yGrid / delGrid;
      for( i=n; i>=0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( boundType == kAllSidesOpen )
   {
      id = 0;
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( boundType == kSpecifyOutlet )
   {
      double xout = infile.ReadItem( xout, "OUTLET_X_COORD" );
      double yout = infile.ReadItem( yout, "OUTLET_Y_COORD" );
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=0, id=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( yout == 0 && xout > dist && xout < dist + delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( xout == xGrid && yout > dist && yout < dist + delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( yout == yGrid && xout < dist && xout > dist - delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
      n = yGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( xout == 0 && yout < dist && yout > dist - delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
   }
   bndList.makeCircular();
   cout << "made points; now adding edges\n";
   
   tPtrListIter< tSubNode > bndIter( bndList );
   for( node0 = bndIter.FirstP(); !( bndIter.AtEnd() ); node0 = bndIter.NextP() )
   {
      node1 = bndIter.ReportNextP();
      node2 = bndIter.ReportPrevP();
      AddEdge( node0, node1, node2 );
   }
   nnodes = nodeList.getSize();
   nedges = edgeList.getSize();
   ntri = 0;
     //DumpEdges:
     /*cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }*/
   cout << "calling repair mesh for initial boundary\n";
   assert( RepairMesh( bndList ) );
   cout << "filling in points\n";
   
     //FILL IN POINTS
   if( ptPlace == kUniformGrid || ptPlace == kPerturbedGrid )
   {
      nx = xGrid / delGrid;
      ny = yGrid / delGrid;
      for( i=1; i<nx; i++ )
      {
         for( j=1; j<ny; j++, id++ )
         {
            xyz[0] = i * delGrid;
            xyz[1] = j * delGrid;
            if( ptPlace == kPerturbedGrid )
            {
               xyz[0] += 0.01 * delGrid * ( ran3( &seed ) - 0.5 );
               xyz[1] += 0.01 * delGrid * ( ran3( &seed ) - 0.5 );
            }
            xyz[2] = mElev + mElev * 0.01 * ( ran3( &seed ) - 0.5 );
            AddNodeAt( xyz );
         }
      }
   }
   else if( ptPlace == kRandomGrid )
   {
      for( i=0; i<numPts; i++ )
      {
         xyz[0] = ran3(&seed) * xGrid;
         xyz[1] = ran3(&seed) * yGrid;
         xyz[2] = mElev + mElev * 0.01 * ( ran3( &seed ) - 0.5 );
         if( xyz[0] != 0 && xyz[0] != xGrid && xyz[1] != 0 && xyz[1] != yGrid )
             AddNodeAt( xyz );
      }
   }
   MakeCCWEdges();
   CheckMeshConsistency();
}

template< class tSubNode >
tGrid< tSubNode >::
~tGrid() {cout << "    ~tGrid()" << endl;}                    //tGrid
 


/*****************************************************************************\
**
**  CheckMeshConsistency
**
**  Performs a series of tests to make sure the mesh connectivity is correct.
**  Should be called immediately after reading in a user-defined mesh (but
**  of course can also be used for debugging).
**
**  The consistency checks include the following:
**
**  1) Each edge:
**     - Has valid origin and destination pointers
**     - Has a valid counter-clockwise edge, which shares the same origin but
**       not the same destination
**     - Is paired with its complement in the list
**
**  2) Each node:
**     - Points to a valid edge which has the node as its origin
**     - If the node is not a boundary, it has at least one neighbor that
**       is not a closed boundary
**     - Has a consistent spoke list (ie, you can go around the spokes and
**       get back to where you started)
**
**  3) Each triangle:
**     - Has 3 valid points and edges
**     - Each edge Ei has Pi as its origin and P((i+1)%3) as its
**       destination
**     - If an opposite triangle Ti exists, points P((i+1)%3) and
**       P((i+2)%3) are the same as points PO((n+2)%3) and PO((n+1)%3) in
**       the opposite triangle, where PO denotes a point in the opposite
**       triangle and n is the vertex ID (0, 1, or 2) of the point in the
**       opposite triangle that is opposite from the shared face.
**     - If an opposite triange Ti does not exist, points P((i+1)%3) and
**       and P((i+2)%3) should both be boundary points.
**
**      Data members updated: none
**      Called by: 
**      Calls:
**      Notes: does not check whether ID's are within the valid range;
**             that's assumed to be taken care of by the input routines
**        
**      Created: GT 1/98
**
\*****************************************************************************/
#define kMaxSpokes 100
template<class tSubNode>
void tGrid< tSubNode >::
CheckMeshConsistency( void )
{
   tGridListIter<tSubNode> nodIter( nodeList );
   tGridListIter<tEdge> edgIter( edgeList );
   tListIter<tTriangle> triIter( triList );
   tNode * cn, * org, * dest;
   tEdge * ce, * cne, * ccwedg;
   tTriangle * ct, * optr;
   int boundary_check_ok, i, nvop;

   // Edges: make sure complementary pairs are together in the list
   // (each pair Ei and Ei+1, for i=0,2,4,...nedges-1, should have the same
   // endpoints but the opposite orientation)
   for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() )
   {
      cne = edgIter.NextP();
      if( ce->getOriginPtrNC() != cne->getDestinationPtrNC()
          || ce->getDestinationPtrNC() != cne->getOriginPtrNC() )
      {
          cerr << "EDGE #" << ce->getID()
               << " must be followed by its complement in the list\n";
          goto error;
      }
      
   }

   // Edges: check for valid origin, destination, and ccwedg
   for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() )
   {
      if( !(org=ce->getOriginPtrNC() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not have a valid origin point\n";
         goto error;
      }
      if( !(dest=ce->getDestinationPtrNC() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not have a valid destination point\n";
         goto error;
      }
      if( !(ccwedg=ce->GetCCWEdg() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not point to a valid counter-clockwise edge\n";
         goto error;
      }
      if( ccwedg->getOriginPtrNC()!=org )
      {
         cerr << "EDGE #" << ce->getID()
              << " points to a CCW edge with a different origin\n";
         goto error;
      }
      if( ccwedg->getDestinationPtrNC()==dest )
      {
         cerr << "EDGE #" << ce->getID()
              << " points to a CCW edge with the same destination\n";
         goto error;
      }

   }
   cout << "EDGES PASSED\n";

   // Nodes: check for valid edg pointer, spoke connectivity, and connection
   // to at least one non-boundary or open boundary node
   for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP() )
   {
      // edg pointer
      if( !(ce = cn->GetEdg()) )
      {
         cerr << "NODE #" << cn->getID()
              << " does not point to a valid edge\n";
         goto error;
      }
      if( ce->getOriginPtrNC()!=cn )
      {
         cerr << "NODE #" << cn->getID()
              << " points to an edge that has a different origin\n";
         goto error;
      }

      // Boundary check and spoke consistency: if node is NOT a boundary,
      // it should be adjacent to at least one non-boundary or open boundary
      // point
      boundary_check_ok = ( cn->getBoundaryFlag()==kNonBoundary ) ? 0 : 1;
      i = 0;
      // Loop around the spokes until we're back at the beginning
      do
      {
         if( ce->getDestinationPtrNC()->getBoundaryFlag()!=kClosedBoundary )
             boundary_check_ok = 1;  // OK, there's at least one open nbr
         i++;
         if( i>kMaxSpokes ) // Uh-oh, an infinite loop
         {
            cerr << "NODE #" << cn->getID()
                 << ": infinite loop in spoke connectivity\n";
            goto error;
         }
      } while( (ce=ce->GetCCWEdg())!=cn->GetEdg() );
      if( !boundary_check_ok )
      {
         cerr << "NODE #" << cn->getID()
              << " is surrounded by closed boundary nodes\n";
         goto error;
      }
      
   }
   cout << "NODES PASSED\n";
   
   // Triangles: check for valid points and connectivity
   for( ct=triIter.FirstP(); !(triIter.AtEnd()); ct=triIter.NextP() )
   {
      for( i=0; i<=2; i++ )
      {
         // Valid point i?
         if( !(cn=ct->pPtr(i)) )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << " has an invalid point " << i << endl;
            goto error;
         }
         // Valid edge i?
         if( !(ce=ct->ePtr(i)) )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << " has an invalid edge " << i << endl;
            goto error;
         }
         // Edge and point consistency
         if( ce->getOriginPtrNC()!=cn )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << ": edge " << i << " does not have point " << i
                 << " as origin\n";
            goto error;
         }
         if( ce->getDestinationPtrNC()!=ct->pPtr((i+1)%3) )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << ": edge " << i << " does not have point " << (i+1)%3
                 << " as destination\n";
            goto error;
         }
         // Opposite triangle: if it exists, check common points
         if( (optr = ct->tPtr(i)) )
         {
            nvop = optr->nVOp(ct); // Num (0,1,2) of opposite vertex in optr
            if( ct->pPtr((i+1)%3) != optr->pPtr((nvop+2)%3)
                || ct->pPtr((i+2)%3) != optr->pPtr((nvop+1)%3) )
            {
               cerr << "TRIANGLE #" << ct->getID()
                    << ": opposite triangle " << i << " does not share nodes "
                    << (ct->pPtr((i+1)%3))->getID() << " and "
                    << (ct->pPtr((i+2)%3))->getID() << endl;
               goto error;
            }
         }
         // If no opposite triangle, make sure it really is a boundary
         else
         {
            if( (ct->pPtr((i+1)%3))->getBoundaryFlag()==kNonBoundary
                || (ct->pPtr((i+2)%3))->getBoundaryFlag()==kNonBoundary )
            {
               cerr << "TRIANGLE #" << ct->getID()
                    << ": there is no neighboring triangle opposite node "
                    << cn->getID() << " but one (or both) of the other nodes "
                    << "is a non-boundary point\n";
               goto error;
            }
         }       
      }
   }
   cout << "TRIANGLES PASSED\n";

   return;
   
  error:
   ReportFatalError( "Error in mesh consistency." );
   
}
#undef kMaxSpokes

template< class tSubNode >
void tGrid< tSubNode >::
Print()                                                  //tGrid
{
   triList.print();
   nodeList.print();
   edgeList.print();
}

template< class tSubNode >
void tGrid< tSubNode >::
MakeCCWEdges()
{
   tGridListIter< tSubNode > nodIter( nodeList );
   tPtrListIter< tEdge > spokIter;
   tSubNode *cn;
   tEdge *ce, *ccwe;
   for( cn = nodIter.FirstP(); !( nodIter.AtEnd() ); cn = nodIter.NextP() )
   {
      spokIter.Reset( cn->getSpokeListNC() );
      ce = spokIter.FirstP();
      assert( ce != 0 );
      cn->SetEdg( ce );
      for( ; !(spokIter.AtEnd()); ce = spokIter.NextP() )
      {
         ccwe = spokIter.ReportNextP();
         assert( ccwe != 0 );
         ce->SetCCWEdg( ccwe );
      }
   }
}

template< class tSubNode >
int tGrid< tSubNode >::
DeleteNode( tListNode< tSubNode > *nodPtr )
{
   //cout << "DeleteNode: " << nodPtr->getDataPtr()->getID() << endl;
   int i;
   tPtrList< tSubNode > nbrList;
     //tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode nodeVal;
   if( !( ExtricateNode( nodPtr->getDataPtrNC(), nbrList ) ) ) return 0;
   nbrList.makeCircular();
   if( nodPtr->getDataRefNC().getBoundaryFlag() )
   {
      nodeList.moveToBack( nodPtr );
      nodeList.removeFromBack( nodeVal );
   }
   else
   {
      nodeList.moveToFront( nodPtr );
      nodeList.removeFromFront( nodeVal );
   }
   
     //cout << "Removed node " << nodeVal.getID() << " at x, y "
     //   << nodeVal.getX() << ", " << nodeVal.getY() << "; " << endl;
   nnodes = nodeList.getSize();
   nedges = edgeList.getSize();
   ntri = triList.getSize();
   tPtrListIter< tSubNode > nbrIter( nbrList );
     /*cout << "leaving hole defined by " << endl << "   Node  x  y " << endl;
   for( i=0, nbrIter.First(); nbrIter.NextIsNotFirst(); i++ )
   {
      if( i>0 ) nbrIter.Next();
      cout << "   " << nbrIter.DatPtr()->getID() << "     "
           << nbrIter.DatPtr()->getX() << "  "
           << nbrIter.DatPtr()->getY() << endl;
   }*/
   if( !RepairMesh( nbrList ) ) return 0;
   //reset node id's
   tGridListIter< tSubNode > nodIter( nodeList );
   assert( nodIter.First() );
   i = 0;
   do
   {
      nodIter.DatRef().setID( i );
      i++;
   }
   while( nodIter.Next() );
   
   
     //cout << "Mesh repaired" << endl;
   return 1;
}

template< class tSubNode >
int tGrid< tSubNode >::
DeleteNode( tSubNode *node )
{
   //cout << "DeleteNode: " << node->getID() << endl;
   int i;
   tPtrList< tSubNode > nbrList;
     //tGridListIter< tSubNode > nodIter( nodeList );
   tListNode< tSubNode > *nodPtr;
   tGridListIter< tSubNode > nodIter( nodeList );
   nodIter.Get( node->getID() );
   nodPtr = nodIter.NodePtr();
   tSubNode nodeVal;
   if( !( ExtricateNode( node, nbrList ) ) ) return 0;
   nbrList.makeCircular();
   if( node->getBoundaryFlag() )
   {
      nodeList.moveToBack( nodPtr );
      nodeList.removeFromBack( nodeVal );
   }
   else
   {
      nodeList.moveToFront( nodPtr );
      nodeList.removeFromFront( nodeVal );
   }
   
     //cout << "Removed node " << nodeVal.getID() << " at x, y "
     //   << nodeVal.getX() << ", " << nodeVal.getY() << "; " << endl;
   nnodes = nodeList.getSize();
   nedges = edgeList.getSize();
   ntri = triList.getSize();
   tPtrListIter< tSubNode > nbrIter( nbrList );
     /*cout << "leaving hole defined by " << endl << "   Node  x  y " << endl;
   for( i=0, nbrIter.First(); nbrIter.NextIsNotFirst(); i++ )
   {
      if( i>0 ) nbrIter.Next();
      cout << "   " << nbrIter.DatPtr()->getID() << "     "
           << nbrIter.DatPtr()->getX() << "  "
           << nbrIter.DatPtr()->getY() << endl;
   }*/
   if( !RepairMesh( nbrList ) ) return 0;
   //reset node id's
   assert( nodIter.First() );
   i = 0;
   do
   {
      nodIter.DatRef().setID( i );
      i++;
   }
   while( nodIter.Next() );
     //cout << "Mesh repaired" << endl;
   return 1;
}


template< class tSubNode >
int tGrid< tSubNode >::
ExtricateNode( tSubNode *node, tPtrList< tSubNode > &nbrList )
{
   //cout << "ExtricateNode: " << node->getID() << endl;
   tPtrListIter< tEdge > spokIter( node->getSpokeListNC() );
   tEdge edgeVal1, edgeVal2;
   tSubNode *nodePtr;
     //cout << "Removing spokes: " << endl;
     //assert( ExtricateEdge( edgptrIter.DatPtr() ) );
   do
   {
      assert( spokIter.First() );
      nodePtr = ( tSubNode * ) spokIter.DatPtr()->getDestinationPtrNC();
      nbrList.insertAtBack( nodePtr );
      DeleteEdge( spokIter.DatPtr() );
   }  
   while( !(node->getSpokeList().isEmpty()) );
   nnodes--;
   if( node->getSpokeList().isEmpty() ) return 1;
   return 0;
}
   
template< class tSubNode >
int tGrid< tSubNode >::
DeleteEdge( tEdge * edgePtr )
{
   //cout << "DeleteEdge(...) " << edgePtr->getID() << endl;
   tEdge edgeVal1, edgeVal2;
   if( !ExtricateEdge( edgePtr ) ) return 0;
   if( edgePtr->getBoundaryFlag() )
   {
      assert( edgeList.removeFromBack( edgeVal1 ) );
      assert( edgeList.removeFromBack( edgeVal2 ) );
   }
   else
   {
      assert( edgeList.removeFromFront( edgeVal1 ) );
      assert( edgeList.removeFromFront( edgeVal2 ) );
   }
     //cout << "  edges " << edgeVal1.getID() << " and "
     //   <<  edgeVal2.getID() << " between nodes "
     //   << edgeVal1.getOriginPtr()->getID() << " and "
     //   << edgeVal2.getOriginPtr()->getID() << " removed" << endl;
   if( &edgeVal1 == 0 || &edgeVal2 == 0 ) return 0;
   return 1;
}

template< class tSubNode >
int tGrid< tSubNode >::
ExtricateEdge( tEdge * edgePtr )
{
   //cout << "ExtricateEdge: " << edgePtr->getID() << endl;
   assert( edgePtr != 0 );
     //temporary objects:
   tEdge *tempedgePtr, *ce, *cce, *spk;
   tGridListIter< tEdge > edgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   tPtrList< tEdge > *spkLPtr;
   tListNode< tEdge > *listnodePtr;
   tTriangle triVal1, triVal2;
   tArray< tTriangle * > triPtrArr(2);
     //cout << "find edge in list; " << flush;
   ce = edgIter.GetP( edgePtr->getID() );

     //cout << "update origin's spokelist if not done already; " << flush;
   spkLPtr = &( ce->getOriginPtrNC()->getSpokeListNC() );
   spokIter.Reset( *spkLPtr );
   for( spk = spokIter.FirstP(); spk != ce && !( spokIter.AtEnd() ); spk = spokIter.NextP() );
   if( spk == ce )
   {
      spk = spokIter.NextP();
      spkLPtr->removePrev( tempedgePtr, spokIter.NodePtr() );
   }
     //cout << "find triangle; " << flush;
   triPtrArr[0] = TriWithEdgePtr( edgePtr );
   listnodePtr = edgIter.NodePtr();
   assert( listnodePtr != 0 );
     //cout << "find compliment; " << flush;
   if( edgePtr->getID()%2 == 0 ) cce = edgIter.NextP();
   else if( edgePtr->getID()%2 == 1 ) cce = edgIter.PrevP();
   else return 0;
     //cout << "find other triangle; " << flush;
   triPtrArr[1] = TriWithEdgePtr( cce );
     //if triangles exist, delete them
   if( triPtrArr[0] != 0 )
       if( !DeleteTriangle( triPtrArr[0] ) ) return 0;
   if( triPtrArr[1] != 0 )
       if( !DeleteTriangle( triPtrArr[1] ) ) return 0;
     //update compliment's origin's spokelist
   spkLPtr = &(cce->getOriginPtrNC()->getSpokeListNC());
   spokIter.Reset( *spkLPtr );
   for( spk = spokIter.FirstP(); spk != cce && !( spokIter.AtEnd() );
        spk = spokIter.NextP() );
   if( spk == cce )
   {
      spk = spokIter.NextP();
      spkLPtr->removePrev( tempedgePtr, spokIter.NodePtr() );
   }
   if( ce->getBoundaryFlag() )
   {
        //move edges to back of list
      edgeList.moveToBack( listnodePtr );
      edgeList.moveToBack( edgIter.NodePtr() );
   }
   else
   {
        //move edges to front of list
      edgeList.moveToFront( edgIter.NodePtr() );
      edgeList.moveToFront( listnodePtr );
   }
   nedges-=2;
   return 1;
}

template< class tSubNode >
tTriangle * tGrid< tSubNode >::
LocateTriangle( double x, double y )
{
   //cout << "LocateTriangle" << endl;
   int n, lv=0;
   tListIter< tTriangle > triIter( triList );  //lt
   tTriangle *lt = &(triIter.DatRef());
   double a, b, c;
   tArray< double > xy1, xy2;
     /* it starts from the first triangle, 
        searches through the triangles until the point is on
        the same side of all the edges of a triangle.
        "lt" is the current triangle and "lv" is the edge number.
        */
   for (n=0 ;(lv!=3)&&(lt); n++)
   {
      xy1 = lt->pPtr(lv)->get2DCoords();
      xy2 = lt->pPtr( (lv+1)%3 )->get2DCoords();
      a = (xy1[1] - y) * (xy2[0] - x);
      b = (xy1[0] - x) * (xy2[1] - y);
      c = a - b;
      if ( c > 0.0 )
      {
         lt=lt->tPtr( (lv+2)%3 );
         lv=0;
      }
      else {lv++;}
        //if( !(n < ntri) )
      if( lt != 0 )
          cout << "find tri for point w/ x, y, " << x << ", " << y
               << "; no. tri's " << ntri << "; now at tri " << lt->getID() << endl;
      if( n >= ntri + 20 )
      {
         DumpTriangles();
         DumpNodes();
      }
      cout << flush;
      assert( n < ntri + 20 );
   }
   return(lt);
}

template< class tSubNode >
tTriangle * tGrid< tSubNode >::
LocateNewTriangle( double x, double y )
{
   //cout << "LocateTriangle" << endl;
   int n, lv=0;
   tListIter< tTriangle > triIter( triList );  //lt
   tTriangle *lt = triIter.FirstP();
   tSubNode *p1, *p2;

   tArray< double > xy1, xy2;
     /* it starts from the first triangle, 
        searches through the triangles until the point is on
        the same side of all the edges of a triangle.
        "lt" is the current triangle and "lv" is the edge number.
        */
   for (n=0 ;(lv!=3)&&(lt); n++)
   {
      p1 = (tSubNode *) lt->pPtr(lv);
      if( p1->Meanders() ) xy1 = p1->getNew2DCoords();
      else xy1 = p1->get2DCoords();
      p2 = (tSubNode *) lt->pPtr( (lv+1)%3 );
      if( p2->Meanders() ) xy1 = p2->getNew2DCoords();
      else xy2 = p2->get2DCoords();
      if ( ( (xy1[1] - y) * (xy2[0] - x) ) > ( (xy1[0] - x) * (xy2[1] - y)) )
      {
         lt=lt->tPtr( (lv+2)%3 );
         lv=0;
      }
      else {lv++;}
        /*if( !(n < ntri) )
          cout << "tri not found for point w/ x, y, " << x << ", " << y
               << "; no. tri's " << ntri << "; now at tri " << lt->getID() << endl;*/
        //assert( n < ntri + 20 );
   }
   return(lt);
}

template< class tSubNode >
tTriangle *tGrid< tSubNode >::
TriWithEdgePtr( tEdge *edgPtr )
{
   assert( edgPtr != 0 );
   tTriangle *ct;
     //cout << "TriWithEdgePtr " << edgPtr->getID();
   tListIter< tTriangle > triIter( triList ); 
   for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
       if( ct != 0 )
           if( ct->ePtr(0) == edgPtr ||
               ct->ePtr(1) == edgPtr ||
               ct->ePtr(2) == edgPtr ) return ct;
   return 0;
}

template< class tSubNode >
int tGrid< tSubNode >::
DeleteTriangle( tTriangle * triPtr )
{
   //cout << "DeleteTriangle(...) " << triPtr->getID() << endl;
   tTriangle triVal;
   if( !ExtricateTriangle( triPtr ) ) return 0;
   if( !( triList.removeFromFront( triVal ) ) ) return 0;
   if( &triVal == 0 ) return 0;
   return 1;
}
      
template< class tSubNode >
int tGrid< tSubNode >::
ExtricateTriangle( tTriangle *triPtr )
{
   //cout << "ExtricateTriangle" << endl;
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct;
   for( ct = triIter.FirstP(); ct != triPtr && !( triIter.AtEnd() ); ct = triIter.NextP() );
   if( ( triIter.AtEnd() ) ) return 0;
   int i, j;
   for( i=0; i<3; i++ ) for( j=0; j<3; j++ )
       if( triPtr->tPtr(i) != 0 )
           if( triPtr->tPtr(i)->tPtr(j) == triPtr ) triPtr->tPtr(i)->setTPtr( j, 0 );
   triList.moveToFront( triIter.NodePtr() );
   ntri--;
   return 1;
}

/*template< class tSubNode >
int tGrid< tSubNode >::
MakeMesh( tPtrList< tSubNode > &bndList )
{
   tPtrListIter< tSubNode > bndIter( bndList );
   tPtrList< tSubNode > nptrList;
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn, *cnn;
   for( cn = nodIter.FirstP(); !( nodIter.AtEnd() ); cn = nodIter.NextP() )
   {
      nptrList.insertAtBack( cn );
   }
   for( cn = bndIter.FirstP(); bndIter.NodePtr() != 0; cn = bndIter.NextP() )
   {
      cnn = bndIter.ReportNextP();
   }
   
}*/

template< class tSubNode >
int tGrid< tSubNode >::
RepairMesh( tPtrList< tSubNode > &nbrList )
{
   assert( &nbrList != 0 );
     //cout << "RepairMesh: " << endl;
   if( nbrList.getSize() < 3 ) return 0;
   int flowflag, i, j;
   tSubNode * gridnodePtr;
   nbrList.makeCircular();
   tPtrListIter< tSubNode > nbrIter( nbrList );
   
   while( nbrList.getSize() > 3 )
   {
        //cout << "in loop, nbr size = " << nbrList.getSize() << endl;
      
      flowflag = 1;
      if( Next3Delaunay( nbrList, nbrIter ) ) //checks for ccw and Del.
      {
           //cout << "found 3 Delaun!\n";
         
         AddEdgeAndMakeTriangle( nbrList, nbrIter );
           //remove "closed off" pt
         nbrList.removeNext( gridnodePtr, nbrIter.NodePtr() );
      }
        //else cout << "Not delaun\n";
      
      nbrIter.Next();                    //step forward once in nbrList
   }
   assert( nbrList.getSize() == 3 );
   assert( ntri == triList.getSize() );
   assert( nedges == edgeList.getSize() );
   assert( nnodes == nodeList.getSize() );       //make sure numbers are right
   MakeTriangle( nbrList, nbrIter );             //make final triangle
     //do some checking?
     //cout << "done" << endl;
   return 1;
}

//vertices of tri in ccw order; edges are added between node1 and node2
template< class tSubNode >
int tGrid< tSubNode >::
AddEdge( tSubNode *node1, tSubNode *node2, tSubNode *node3 ) 
{
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
     //cout << "AddEdge" << endl;
     //"between nodes " << node1->getID()
     //   << " and " << node2->getID() << " w/ ref to node " << node3->getID() << endl;
   int flowflag = 1;
   int i, j, newid;
   tEdge tempEdge1, tempEdge2;
   tEdge *ce, *le;
   tGridListIter< tEdge > edgIter( edgeList );
   tGridListIter< tSubNode > nodIter( nodeList );
   tPtrListIter< tEdge > spokIter;
        //DumpEdges:
     /*cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }*/

     //deal with new edges and new triangle:
   tempEdge1.setOriginPtr( node1 );   //set edge1 ORG
   tempEdge2.setDestinationPtr( node1 );//set edge2 DEST
   if( node1->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   tempEdge2.setOriginPtr( node2 );   //set edge2 ORG
   tempEdge1.setDestinationPtr( node2 );//set edge1 DEST
   if( node2->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   if( !( edgeList.isEmpty() ) )
       newid = edgIter.LastP()->getID() + 1;
   else newid = 0;
   tempEdge1.setID( newid );                     //set edge1 ID
   newid++;
   tempEdge2.setID( newid );                     //set edge2 ID
   tempEdge1.setFlowAllowed( flowflag );         //set edge1 FLOWALLOWED
   tempEdge2.setFlowAllowed( flowflag );         //set edge2 FLOWALLOWED
   if( flowflag == 1 )
   {
      edgeList.insertAtActiveBack( tempEdge1 );  //put edge1 active in list
      edgeList.insertAtActiveBack( tempEdge2 );  //put edge2 active in list
      le = edgIter.LastActiveP();                      //set edgIter to lastactive
   }
   else
   {
      edgeList.insertAtBack( tempEdge1 );        //put edge1 in list
      edgeList.insertAtBack( tempEdge2 );        //put edge2 in list
      le = edgIter.LastP();                            //set edgIter to last
   }
     //add pointers to the new edges to nodes' spokeLists:
   spokIter.Reset( node2->getSpokeListNC() );
     //cout << "node " << node2->getID() << ": ";
   if( node2->getSpokeListNC().isEmpty() )
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " in otherwise empty list" << endl;
      node2->insertFrontSpokeList( le);
      node2->getSpokeListNC().makeCircular();
   }
   else if( spokIter.ReportNextP() == spokIter.DatPtr() )
   {
      node2->insertFrontSpokeList( le);
        //node2->getSpokeListNC().makeCircular();
   }
   else
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " w/ reference to node " << node3->getID() << endl;
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
      node2->getSpokeListNC().insertAtNext( le,
                                            spokIter.NodePtr() ); //put edge2 in SPOKELIST
   }
   spokIter.Reset( node1->getSpokeListNC() );
   le = edgIter.PrevP();                     //step backward once in edgeList
     //cout << "node " << node1->getID() << ": ";
   if( node1->getSpokeListNC().isEmpty() )
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " in otherwise empty list" << endl;
      node1->insertFrontSpokeList( le );
      node1->getSpokeListNC().makeCircular();
   }
   else if( spokIter.ReportNextP() == spokIter.DatPtr() )
   {
      node1->insertFrontSpokeList( le );
        //node2->getSpokeListNC().makeCircular();
   }
   else
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " w/ reference to node " << node3->getID() << endl;
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
      node1->getSpokeListNC().insertAtPrev( le,
                                            spokIter.NodePtr() );//put edge1 in SPOKELIST
   }
   
   nedges+=2;
   //cout << "2 edges added" << endl;
   //reset edge id's
   for( ce = edgIter.FirstP(), i = 0; !( edgIter.AtEnd() ); ce = edgIter.NextP(), i++ )
   {
      ce->setID( i );
   }
   return 1;
}

   
template< class tSubNode >
int tGrid< tSubNode >::
AddEdgeAndMakeTriangle( tPtrList< tSubNode > &nbrList,
                        tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) );
     //cout << "AddEdgeAndMakeTriangle" << endl;
     //cout << "aemt 0\n" << flush;
   int flowflag = 1;
     //cout << "aemt 0.1\n" << flush;
   int i, j, newid;
   tSubNode *cn, *cnn, *cnnn;
   cn = nbrIter.DatPtr();
   cnn = nbrIter.NextP();
   cnnn = nbrIter.ReportNextP();
   nbrIter.Prev();
   tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
       p2( cnnn->get2DCoords() );
   if( !PointsCCW( p0, p1, p2 ) )
       cout << "in AEAMT nodes not CCW: " << cn->getID() << ", "
            << cnn->getID() << ", " << cnnn->getID() << endl;
     //cout << "aemt 1/2\n" << flush;
   tEdge tempEdge1, tempEdge2, *ce, *le;
   tGridListIter< tEdge > edgIter( edgeList );
   tTriangle tempTri;
   tTriangle *nbrtriPtr, *ct;
   tListIter< tTriangle > triIter( triList );
   tPtrListIter< tEdge > spokIter;
        //DumpEdges:
     /*cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }*/

     //deal with new edges and new triangle:
     //cout << "aemt 1\n" << flush;
   tempEdge1.setOriginPtr( cn );                 //set edge1 ORG
   tempEdge2.setDestinationPtr( cn );            //set edge2 DEST
   if( cn->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   tempTri.setPPtr(0, cn );                      //set tri VERTEX ptr 0
   tempTri.setPPtr(1, cnn );                     //set tri VERTEX ptr 1
   tempTri.setPPtr(2, cnnn );                    //set tri VERTEX ptr 2
   tempEdge2.setOriginPtr( cnnn );               //set edge2 ORG
   tempEdge1.setDestinationPtr( cnnn );          //set edge1 DEST
   if( cnnn->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   le = edgIter.LastP();
   newid = le->getID() + 1;
   tempEdge1.setID( newid );                     //set edge1 ID
   newid++;
   tempEdge2.setID( newid );                     //set edge2 ID
   tempEdge1.setFlowAllowed( flowflag );         //set edge1 FLOWALLOWED
   tempEdge2.setFlowAllowed( flowflag );         //set edge2 FLOWALLOWED
     //cout << "aemt 2\n" << flush;
   if( flowflag == 1 )
   {
      edgeList.insertAtActiveBack( tempEdge1 );  //put edge1 active in list
      edgeList.insertAtActiveBack( tempEdge2 );  //put edge2 active in list
      le = edgIter.LastActiveP();                //set edgIter to lastactive
   }
   else
   {
      edgeList.insertAtBack( tempEdge1 );        //put edge1 in list
      edgeList.insertAtBack( tempEdge2 );        //put edge2 in list
      le = edgIter.LastP();                      //set edgIter to last
   }
   tempTri.setEPtr(2, le );                      //set tri EDGE ptr 2
     //add pointers to the new edges to nodes' spokeLists:
     //cout << "aemt 3\n" << flush;
   spokIter.Reset( cnnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
   cnnn->getSpokeListNC().insertAtPrev( le, spokIter.NodePtr() );//put edge2 in SPOKELIST
   le = edgIter.PrevP();
   spokIter.Reset( cn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
   cn->getSpokeListNC().insertAtNext( le, spokIter.NodePtr() ); //put edge1 in SPOKELIST
   tempTri.setEPtr(0, ce );                      //set tri EDGE ptr 0
   spokIter.Reset( cnn->getSpokeListNC() );
     //cout << "aemt 4\n" << flush;
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   tempTri.setEPtr(1, ce );                      //set tri EDGE ptr 1
   tempTri.setID( ntri );                        //set tri ID
   for( i=0; i<3; i++ ) tempTri.setTPtr( i, 0 ); //initialize tri TRI ptrs to nul
   triList.insertAtBack( tempTri );              //put tri in list
   ct = triIter.LastP();                         //set triIter to last
   spokIter.Reset( cnnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );
     //cout << "aemt 5\n" << flush;
   nbrtriPtr = TriWithEdgePtr( ce );
   ct->setTPtr(0, nbrtriPtr );                   //set tri TRI ptr 0
   if( nbrtriPtr != 0 )
   {
      for( i=0; nbrtriPtr->ePtr(i) != ce && i<3; i++ );
      assert( i<3 );
      nbrtriPtr->setTPtr( (i+2)%3, ct );         //set NBR TRI ptr to tri
   }
   spokIter.Reset( cnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );
   nbrtriPtr = TriWithEdgePtr( ce );
   ct->setTPtr(2, nbrtriPtr );                   //set tri TRI ptr 2
     //cout << "aemt 6\n" << flush;
   if( nbrtriPtr != 0 )
   {
      for( i=0; nbrtriPtr->ePtr(i) != ce && i<3; i++ );
      assert( i<3 );
      nbrtriPtr->setTPtr( (i+2)%3, ct );         //set NBR TRI ptr to tri
   }
   ntri++;                                       //increment tGrid::ntri by one
   nedges+=2;                                    //increment tGrid::nedges by two
     //cout << "aemt 7\n" << flush;
   //reset edge id's
   for( ce = edgIter.FirstP(), i=0; !( edgIter.AtEnd() ); ce = edgIter.NextP(), i++ )
       ce->setID( i );
     //cout << "aemt 8\n" << flush;
  //reset triangle id's
   for( ct = triIter.FirstP(), i=0; !( triIter.AtEnd() ); ct = triIter.NextP(), i++ )
   {
      ct->setID( i );
      assert( ( ( ct->tPtr(0) != ct->tPtr(1) && ct->tPtr(0) != ct->tPtr(2) ) ||
                ct->tPtr(0) == 0 ) &&
              ( ( ct->tPtr(1) != ct->tPtr(0) && ct->tPtr(1) != ct->tPtr(2) ) ||
                ct->tPtr(1) == 0 ) &&
              ( ( ct->tPtr(2) != ct->tPtr(0) && ct->tPtr(2) != ct->tPtr(1) ) ||
                ct->tPtr(2) == 0 ) );
   }
   
     //cout << "  finished" << endl << flush;
   return 1;
}

template< class tSubNode >
int tGrid< tSubNode >::
MakeTriangle( tPtrList< tSubNode > &nbrList,
              tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   assert( nbrList.getSize() == 3 );
     //cout << "MakeTriangle" << endl;
   int i, j;
   tTriangle tempTri;
   tTriangle *nbrtriPtr;
   tSubNode *cn, *cnn, *cnnn;
   tEdge *ce, *dce;
   tTriangle *ct;
   tListIter< tTriangle > triIter( triList );
   tGridListIter< tEdge > edgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   assert( nbrList.getSize() == 3 );
   cn = nbrIter.DatPtr();
   cnn = nbrIter.NextP();
   cnnn = nbrIter.NextP();
   nbrIter.Next();
   tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
       p2( cnnn->get2DCoords() );
   if( !PointsCCW( p0, p1, p2 ) )
   {
      cout << "in MT nodes not CCW: " << cn->getID() << ", "
           << cnn->getID() << ", " << cnnn->getID();
      if( cn->Meanders() ) p0 = cn->getNew2DCoords();
      else p0 = cn->get2DCoords();
      if( cnn->Meanders() ) p1 = cnn->getNew2DCoords();
      else p1 = cnn->get2DCoords();
      if( cnnn->Meanders() ) p2 = cnnn->getNew2DCoords();
      else p2 = cnnn->get2DCoords();
      if( !PointsCCW( p0, p1, p2 ) )
          cout << "; nor are new coords CCW ";
      cout << endl;
   }
   
     //for debugging:
     //DumpEdges();
     //DumpSpokes:
     //deal with (last) new triangle:
   ct = triIter.LastP();
   int newid = ct->getID() + 1;
   tempTri.setID( newid );                          //set tri ID
     //set edge and vertex ptrs & add to triList:
   for( i=0; i<3; i++ )
   {
      tempTri.setPPtr( i, cn );       //set tri VERTEX ptr i
      spokIter.Reset( cn->getSpokeListNC() );
      cn = nbrIter.NextP();                     //step forward once in nbrList   
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != cn && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );
      tempTri.setEPtr( i, ce );      //set tri EDGE ptr i
      tempTri.setTPtr( i, 0 );                      //initialize tri TRI ptrs to nul
   }
   triList.insertAtBack( tempTri );       //put tri in list
   ct = triIter.LastP();                               //set triIter to last
   assert( cn == ct->pPtr(0) );                           //make sure we're where we
                                                          //think we are
     //set tri ptrs:
   dce = 0;
   nbrtriPtr == 0;
   for( j=0; j<3; j++ )
   {
      nbrIter.Next();                     //step forward once in nbrList   
      cn = nbrIter.NextP();                     //step forward once in nbrList   
      spokIter.Reset( cn->getSpokeListNC() );
      nbrIter.Next();                     //step forward once in nbrList
      cn = nbrIter.NextP();                     //step forward once in nbrList
      if( j>0 ) dce = ce;
      for( ce = spokIter.FirstP(); ce->getDestinationPtrNC() != cn && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );
        //********BUG: following assertion failed; called from FlipEdge,
        //from CheckForFlip, from CheckLocallyDelaunay, from MoveNodes***************
      if( !( TriWithEdgePtr( ce ) != nbrtriPtr || nbrtriPtr == 0 ) )
      {
         p0 = cn->get2DCoords();
         p1 = cnn->get2DCoords();
         p2 = cnnn->get2DCoords();
         
         if( PointsCCW( p0, p1, p2 ) )
             cout << "something FUNNY going on";
         else cout << "tri not CCW: " << nbrtriPtr->getID() << endl;
      }
      nbrtriPtr = TriWithEdgePtr( ce );
      ct->setTPtr( j, nbrtriPtr );      //set tri TRI ptr j
      if( nbrtriPtr != 0 )
      {
         for( i=0; i<3; i++ )
         {
            assert( nbrtriPtr->ePtr(i) != 0 );
            assert( ce != 0 );
            if( nbrtriPtr->ePtr(i) == ce ) break;
         }
         assert( i < 3 );
         nbrtriPtr->setTPtr( (i+2)%3, ct );//set NBR TRI ptr to tri
      }
   }
   ntri++;
   //reset triangle id's
   for( ct = triIter.FirstP(), i=0; !( triIter.AtEnd() ); ct = triIter.NextP(), i++ )
   {
      ct->setID( i );
   }
   
   return 1;
}

/*  AddNode: add a node with value/properties of referenced node   */
template< class tSubNode >
int tGrid< tSubNode >::
AddNode( tSubNode &nodeRef )
{
   assert( &nodeRef != 0 );
   tArray< double > xyz( nodeRef.get3DCoords() );
   cout << "AddNode at " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;
   tTriangle *tri;
   cout << "locate tri" << endl << flush;
   tri = LocateTriangle( xyz[0], xyz[1] );
   assert( tri != 0 );  //if( tri == 0 ) return 0;
   int i, j, k;
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn;
   int newid = nodIter.LastP()->getID() + 1;
   nodeRef.setID( newid );
   nodeList.insertAtActiveBack( nodeRef );
   assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
     //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = (tSubNode *) tri->pPtr(i);
      bndyList.insertAtBack( tmpPtr );
   }
   bndyList.makeCircular();
   //delete triangle
   i = DeleteTriangle( tri );
   assert( i != 0 );  //if ( !DeleteTriangle( tri ) ) return 0;
   //make 3 new triangles
   tPtrListIter< tSubNode > bndyIter( bndyList );
   tSubNode *node3 = bndyIter.FirstP();
   tSubNode *node2 = nodIter.LastActiveP();
   tSubNode *node1 = bndyIter.NextP();
   tSubNode *node4 = bndyIter.NextP();
   tArray< double > p1( node1->get2DCoords() ),
       p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
       p4( node4->get2DCoords() );
   if( xyz.getSize() == 3)
   {
      cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
           << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   else
   {
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) p3 = node3->getNew2DCoords();
      if( node4->Meanders() ) p4 = node4->getNew2DCoords();  
      cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
           << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   AddEdge( node1, node2, node3 );  //add edge between node1 and node2
   tPtrList< tSubNode > tmpList;
   tmpList.insertAtBack( node3 );
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node2 );
   tPtrListIter< tSubNode > tmpIter( tmpList );
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node4 );
   tmpIter.First();
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );
   tmpList.insertAtBack( node4 );
   tmpList.insertAtBack( node3 );
   tmpList.makeCircular();
   tmpIter.First();
   MakeTriangle( tmpList, tmpIter );
   //put 3 resulting triangles in ptr list
   if( xyz.getSize() == 3 )
   {
      cout << "flip checking" << endl;
      tPtrList< tTriangle > triptrList;
      tListIter< tTriangle > triIter( triList );
      tPtrListIter< tTriangle > triptrIter( triptrList );
      tTriangle *ct;
      triptrList.insertAtBack( triIter.LastP() );
      triptrList.insertAtBack( triIter.PrevP() );
      triptrList.insertAtBack( triIter.PrevP() );
        //check list for flips; if flip, put new triangles at end of list
      int flip = 1;
      while( !( triptrList.isEmpty() ) )
      {
         ct = triptrIter.FirstP();
         for( i=0; i<3; i++ )
         {
            if( ct->tPtr(i) != 0 )
            {
               if( CheckForFlip( ct, i, flip ) )
               {
                  triptrList.insertAtBack( triIter.LastP() );
                  triptrList.insertAtBack( triIter.PrevP() );
                  break;
               }
            }
         }
         triptrList.removeFromFront( ct );
      }
   }
   //reset node id's
   for( cn = nodIter.FirstP(), i=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), i++ )
   {
      cn->setID( i );
   }
   cout << "AddNode finished" << endl;
   return 1;
}

/*   AddNodeAt: add a node with referenced coordinates to mesh   */

template< class tSubNode >
tSubNode *tGrid< tSubNode >::
AddNodeAt( tArray< double > &xyz )
{
   assert( &xyz != 0 );
   cout << "AddNodeAt " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;
   tTriangle *tri;
   cout << "locate tri" << endl << flush;
   if( xyz.getSize() == 3 ) tri = LocateTriangle( xyz[0], xyz[1] );
   else tri = LocateNewTriangle( xyz[0], xyz[1] );
   if( tri == 0 ) return 0;
   int i, j, k;
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode tempNode, *cn;
   tempNode.set3DCoords( xyz[0], xyz[1], xyz[2]  );
   if( xyz.getSize() != 3 ) tempNode.setNew2DCoords( xyz[0], xyz[1] );
   tempNode.setBoundaryFlag( 0 );

   int newid = nodIter.LastP()->getID() + 1;
   tempNode.setID( newid );
   nodeList.insertAtActiveBack( tempNode );
   assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
     //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = (tSubNode *) tri->pPtr(i);
      bndyList.insertAtBack( tmpPtr );
   }
   bndyList.makeCircular();
   //delete triangle
   if ( !DeleteTriangle( tri ) ) return 0;
   //make 3 new triangles
   tPtrListIter< tSubNode > bndyIter( bndyList );
   tSubNode *node3 = bndyIter.FirstP();
   tSubNode *node2 = nodIter.LastActiveP();
   tSubNode *node1 = bndyIter.NextP();
   tSubNode *node4 = bndyIter.NextP();
   tArray< double > p1( node1->get2DCoords() ),
       p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
       p4( node4->get2DCoords() );
   if( xyz.getSize() == 3)
   {
      cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
           << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   else
   {
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) p3 = node3->getNew2DCoords();
      if( node4->Meanders() ) p4 = node4->getNew2DCoords();  
      cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
           << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   AddEdge( node1, node2, node3 );  //add edge between node1 and node2
   tPtrList< tSubNode > tmpList;
   tmpList.insertAtBack( node3 );
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node2 );
   tPtrListIter< tSubNode > tmpIter( tmpList );
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node4 );
   tmpIter.First();
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );
   tmpList.insertAtBack( node4 );
   tmpList.insertAtBack( node3 );
   tmpList.makeCircular();
   tmpIter.First();
   MakeTriangle( tmpList, tmpIter );
   //put 3 resulting triangles in ptr list
   if( xyz.getSize() == 3 )
   {
      cout << "flip checking" << endl;
      tPtrList< tTriangle > triptrList;
      tListIter< tTriangle > triIter( triList );
      tPtrListIter< tTriangle > triptrIter( triptrList );
      tTriangle *ct;
      triptrList.insertAtBack( triIter.LastP() );
      triptrList.insertAtBack( triIter.PrevP() );
      triptrList.insertAtBack( triIter.PrevP() );
        //check list for flips; if flip, put new triangles at end of list
      int flip = 1;
      while( !( triptrList.isEmpty() ) )
      {
         ct = triptrIter.FirstP();
         for( i=0; i<3; i++ )
         {
            if( ct->tPtr(i) != 0 )
            {
               if( CheckForFlip( ct, i, flip ) )
               {
                  triptrList.insertAtBack( triIter.LastP() );
                  triptrList.insertAtBack( triIter.PrevP() );
                  break;
               }
            }
         }
         triptrList.removeFromFront( ct );
      }
   }
   //reset node id's
   for( cn = nodIter.FirstP(), i=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), i++ )
   {
      cn->setID( i );
   }
   cout << "AddNodeAt finished" << endl;
   return node2;
}


template <class tSubNode>
tGridList<tEdge> * tGrid<tSubNode>::
GetEdgeList() {return &edgeList;}

template <class tSubNode>
tGridList<tSubNode> * tGrid<tSubNode>::
GetNodeList() {return &nodeList;}

template <class tSubNode>
tList< tTriangle > * tGrid<tSubNode>::
GetTriList() {return &triList;}

template< class tSubNode >
tEdge *tGrid< tSubNode >::
getEdgeCompliment( tEdge *edge )
{
   tGridListIter< tEdge > edgIter( edgeList );
   int edgid = edge->getID();
   assert( edgIter.Get( edgid ) );
   if( edgid%2 == 0 ) return edgIter.GetP( edgid + 1 );
   if( edgid%2 == 1 ) return edgIter.GetP( edgid - 1 );
}

template <class tSubNode>
void tGrid<tSubNode>::
UpdateMesh()
{
   cout << "UpdateMesh()" << endl;
   
   //tListIter<tTriangle> tlist( triList );
   tGridListIter<tEdge> elist( edgeList );
   //tGridListIter<tSubNode> nlist( nodeList );
   tEdge * curedg = 0;
   double len;
   tTriangle * curtri;
   
   // Edge lengths
   curedg = elist.FirstP();
   do
   {
      len = curedg->CalcLength();
      //Xcout << "Edge " << curedg->getID() << " length: " << curedg->getLength() << endl;
      curedg = elist.NextP();
      assert( curedg > 0 ); // failure = complementary edges not consecutive
      curedg->setLength( len );
   } while( curedg=elist.NextP() );
   


// Triangle areas
/*   for( tlist.First(); !tlist.AtEnd(); tlist.Next() )
   {
      curtri = tlist.DatPtr();
      curtri->length_sides();
      curtri->CalcArea();
      curtri = curtri->next;
   }
   */
   
   // Voronoi vertices
   //GetVoronoiVertices();

   

   // Voronoi Areas


}


/*****************************************************************************\
**
**      CheckForFlip: checks whether edge between two triangles should be
**                    flipped; may either check, flip, and report, or just
**                    check and report.
**                    Checks whether the present angle or the possible angle
**                    is greater. Greater angle wins. Also uses flip variable
**                    to determine whether to use newx, newy, or x, y.
**
**      Data members updated: Grid
**      Called by: 
**      Calls:  
**        
**      Created: 8/28/97 SL
**      Modified: 12/16/97 SL                                               
**                                                              
**
\*****************************************************************************/
template< class tSubNode >
int tGrid< tSubNode >::
CheckForFlip( tTriangle * tri, int nv, int flip )
{
   if( tri == 0 ) return 0;
   assert( nv < 3 );
     //cout << "THIS IS CheckForFlip(...) " << tri->getID() << endl;
   tSubNode *node0, *node1, *node2, *node3;
   node0 = ( tSubNode * ) tri->pPtr(nv);
   node1 = ( tSubNode * ) tri->pPtr((nv+1)%3);
   node2 = ( tSubNode * ) tri->pPtr((nv+2)%3);
   tTriangle *triop = tri->tPtr(nv);
   int nvop = triop->nVOp( tri );
   node3 = ( tSubNode * ) triop->pPtr( nvop );
   tArray< double > ptest( node3->get2DCoords() ), p0( node0->get2DCoords() ),
       p1( node1->get2DCoords() ), p2( node2->get2DCoords() );
   if( !flip )
   {
      if( node0->Meanders() ) p0 = node0->getNew2DCoords();
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) ptest = node3->getNew2DCoords();
   }
   if( TriPasses( ptest, p0, p1, p2 ) ) return 0;
   if( flip )                     //and make sure there isn't already an edge?
   {
      if( !PointsCCW( p0, p1, ptest ) || !PointsCCW( p0, ptest, p2 ) )return 0;
      cout << "Flip edge" << endl;
      FlipEdge( tri, triop, nv, nvop );
        /*
      assert( DeleteEdge( tri->ePtr( (nv+1)%3 ) ) );
      tPtrList< tSubNode > nbrList;
      nbrList.insertAtBack( node0 );
      nbrList.insertAtBack( node1 );
      nbrList.insertAtBack( node3 );
      nbrList.insertAtBack( node2 );
      assert( RepairMesh( nbrList ) );*/
   }
     //cout << "finished" << endl;
   return 1;
}

/******************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
FlipEdge( tTriangle * tri, tTriangle * triop ,int nv, int nvop )
{
   //cout << "FlipEdge(...)..." << endl;
   tSubNode *cn;
   tPtrList< tSubNode > nbrList;
   DumpTriangles();
   nbrList.insertAtBack( (tSubNode *) tri->pPtr(nv) );
   nbrList.insertAtBack( (tSubNode *) tri->pPtr((nv+1)%3) );
   nbrList.insertAtBack( (tSubNode *) triop->pPtr( nvop ) );
   nbrList.insertAtBack( (tSubNode *) tri->pPtr((nv+2)%3) );
   nbrList.makeCircular();
   DeleteEdge( tri->ePtr( (nv+1)%3 ) );
   tPtrListIter< tSubNode > nbrIter( nbrList );
   AddEdgeAndMakeTriangle( nbrList, nbrIter );
   nbrIter.First();
   nbrList.removeNext( cn, nbrIter.NodePtr() );
   MakeTriangle( nbrList, nbrIter );
   //cout << "finished" << endl;
}


/*****************************************************************************\
**
**      CheckLocallyDelaunay : updates teh triangulation after moving
**             some points.
**             only uses x and y values, which have already been updated in
**             MoveNodes (frmr PreApply).
**             PREAPPLY SHOULD BE CALLED BEFORE THIS FUNCTION IS CALLED
**      Data members updated: Grid
**      Called by: MoveNodes
**      Calls:  
**        
\*****************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
CheckLocallyDelaunay()
{
   cout << "CheckLocallyDelaunay()" << endl;
   tTriangle *at, * trop[3];
   tPtrList< tTriangle > triPtrList;
   tPtrListIter< tTriangle > triPtrIter( triPtrList );
   tListIter< tTriangle > triIter( triList );
   int i, change, flipped;
   int id0, id1, id2;
   tArray< int > npop(3);
   tSubNode *nodPtr;
   
   int flip = 1;
   //first find a way to go through all the points an the line
   //put each triangle into the stack
     //flipped = TRUE;
     /*do
   {*/
     //flipped = FALSE;
   for( at = triIter.FirstP(); !( triIter.AtEnd() ); at = triIter.NextP() )
   {
      change = FALSE;
      for( i = 0; i < 3; i++ )
      {
         nodPtr = ( tSubNode * ) at->pPtr(i);
         if( nodPtr->Meanders() ) change = TRUE;
      }
      if( change ) triPtrList.insertAtBack( at );
   }

   //check list for flips; if flip, put new triangles at end of list
   tPtrListIter< tTriangle > duptriPtrIter( triPtrList );
   tTriangle *tn, *tp;
   while( !( triPtrList.isEmpty() ) )
   {
      at = triPtrIter.FirstP();
      for( i=0; i<3; i++ )
      {
         if( at->tPtr(i) != 0 )
         {
            tp = at->tPtr(i);
            for( tn = duptriPtrIter.FirstP();
                 duptriPtrIter.ReportNextP() != tp &&
                     !( duptriPtrIter.AtEnd() );
                 tn = duptriPtrIter.NextP() );
            tn = 0;
            if( !( duptriPtrIter.AtEnd() ) )
            {
               tn = duptriPtrIter.ReportNextP();
            }
            if( at->tPtr(0) != 0 ) id0 = at->tPtr(0)->getID();
            else id0 = -1;
            if( at->tPtr(1) != 0 ) id1 = at->tPtr(1)->getID();
            else id1 = -1;
            if( at->tPtr(2) != 0 ) id2 = at->tPtr(2)->getID();
            else id2 = -1;
            cout << "check tri " << at->getID() << " with nbrs "
                 << id0 << ", " << id1
                 << ", and " << id2;
               
            if( tp->tPtr(0) != 0 ) id0 = tp->tPtr(0)->getID();
            else id0 = -1;
            if( tp->tPtr(1) != 0 ) id1 = tp->tPtr(1)->getID();
            else id1 = -1;
            if( tp->tPtr(2) != 0 ) id2 = tp->tPtr(2)->getID();
            else id2 = -1;
            cout << " against tri " << tp->getID() << " with nbrs "
                 << id0 << ", " << id1
                 << ", and " << id2 << endl;
            if( CheckForFlip( at, i, flip ) )
            {
               cout << "flipped tri's, got tri ";
               if( tn != 0 )
                   triPtrList.removeNext( tn, duptriPtrIter.NodePtr() );
               tn = triIter.LastP();
               if( tn->tPtr(0) != 0 ) id0 = tn->tPtr(0)->getID();
               else id0 = -1;
               if( tn->tPtr(1) != 0 ) id1 = tn->tPtr(1)->getID();
               else id1 = -1;
               if( tn->tPtr(2) != 0 ) id2 = tn->tPtr(2)->getID();
               else id2 = -1;
               cout << tn->getID() << " with nbrs "
                    << id0 << ", " << id1
                    << ", and " << id2;
               triPtrList.insertAtBack( tn );
               tn = triIter.PrevP();
               if( tn->tPtr(0) != 0 ) id0 = tn->tPtr(0)->getID();
               else id0 = -1;
               if( tn->tPtr(1) != 0 ) id1 = tn->tPtr(1)->getID();
               else id1 = -1;
               if( tn->tPtr(2) != 0 ) id2 = tn->tPtr(2)->getID();
               else id2 = -1;
               cout << " and tri " << tn->getID() << " with nbrs "
                    << id0 << ", " << id1
                    << ", and " << id2 << endl;
               triPtrList.insertAtBack( tn );
               break;
            }
         }
      }
      triPtrList.removeFromFront( at );
   }
      //for each triangle in the stack
/*      for( at = triPtrIter.FirstP(); !( triPtrIter.AtEnd() );
           at = triPtrIter.NextP() )
      {
         for( i=0; i<3; i++ )
         {
            trop[i] = at->tPtr(i);
            if( trop[i] ) npop[i] = trop[i]->nVOp( at );
            else npop[i] = NULL;
         }
         for( i=0; i<3; i++ )
         {
            if( CheckForFlip( trop[i], npop[i], flip ) )
            {
               flipped = TRUE;
            }
         }
      }
   } while( flipped );*/
   cout << "finished" << endl;
}

/*****************************************************************************\
**
**      IntersectsAnyEdge: returns the first edge in the list which intersects
**                         "edge" or NULL if "edge" intersects no other edges
**      Data members updated: Grid
**      Called by: 
**      Calls:  
**
\*****************************************************************************/
template< class tSubNode >
tEdge *tGrid< tSubNode >::
IntersectsAnyEdge( tEdge * edge )
{
   //cout << "IntersectsAnyEdge( tEdge * edge )..." << endl;
   int i;
   tEdge * ce;
   tGridListIter< tEdge > edgIter( edgeList );
   if( !edge )
   {
      cout<<"IntersectsAnyEdge: Warning: invalid edge"<<endl<<flush;
      return( NULL );
   }
     //cout<<"IAE: nedges "<<nedges<<endl<<flush;
     //cout << "call Intersect for edges " << edge->getID()
     //   << " from nodes " << edge->getOriginPtr()->getID()
     //   << " to " << edge->getDestinationPtr()->getID() << "; " << endl;
   for( ce = edgIter.FirstP(); !(edgIter.AtEnd());
        edgIter.Next(), ce = edgIter.NextP() )
   {
      assert( edgIter.NodePtr()->getNext() != 0 );
      if( edge->getID() != ce->getID() &&
          edge->getID() != getEdgeCompliment( edge )->getID() )
      {
           //cout  << " and " << ce->getID() << " from nodes "
           //    << ce->getOriginPtr()->getID()
           //    << " to " << ce->getDestinationPtr()->getID() << endl;
         if( Intersect( edge, ce ) ) return( ce );
      }
      
   }
   assert( edgIter.AtEnd() );
     /*if( i < nedges - 1 )
       cout<<"IntersectsAnyEdge: Warning: whole list not checked"<<endl<<flush;*/
   return( NULL );
}

/*****************************************************************************\
**
**      CheckTriEdgeIntersect():
**        We want to know if the moving point has passed beyond the polygon
**        defined by its spoke edges; if it has, then we will have edges
**        intersecting one another. In the case where the point has simply
**        passed into one of the 'opposite' triangles, then we can just do a
**        flip operation. In the other case, the remedial action is much more
**        complicated, so we just delete the point and add it again.
**
\*****************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
CheckTriEdgeIntersect()
{
   cout << "CheckTriEdgeIntersect()..." << flush << endl;
   int i, j, nv, nvopp, id0, id1, id2;
   int flipped = TRUE;
   int crossed;
   tSubNode *subnodePtr, tempNode, newNode;  
   tEdge * newedg, * cedg, * ccedg, *fedg, *ce, *cex;
   tTriangle * ct, * ctop, *rmtri, *tri;
   tListIter< tTriangle > triIter( triList );
   tGridListIter< tEdge > edgIter( edgeList );
   tGridListIter< tSubNode > nodIter( nodeList );
   tGridListIter< tEdge > xedgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   tGridList< tSubNode > tmpNodeList;
   tGridListIter< tSubNode > tmpIter( tmpNodeList );
   tArray< double > p0, p1, p2, xy, xyz, xy1, xy2;
   tSubNode *cn, *vtx;
   tPtrList< tTriangle > triptrList;
   tPtrListNode< tTriangle > *tpListNode;
   tPtrListIter< tTriangle > tpIter( triptrList );
     //check for triangles with edges which intersect (an)other edge(s)
     //newedg = new tEdge;
   while( flipped )
   {
      flipped = FALSE;
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      {
         for( i=0; i<3; i++ )
         {
            cn = (tSubNode *) ct->pPtr(i);
            if( cn->Meanders() ) break;
         }
         if( i!=3 ) triptrList.insertAtBack( ct );
      }
        //for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      for( ct = tpIter.FirstP(); !(triptrList.isEmpty());
           triptrList.removeFromFront( ct ), ct = tpIter.FirstP() )
      {
           //cout<<"PA: check triangle "<<ct->id<<", w edges "
           //<<ct->e[0]->id<<", "<<ct->e[1]->id<<", "<<ct->e[2]->id<<endl<<flush;
         if( !NewTriCCW( ct ) )
         {
            flipped = TRUE;
            for( i=0, j=0; i<3; i++ )
            {
               if( ct->pPtr(i)->getBoundaryFlag() != kNonBoundary ) j++;
            }
            if( j > 1 )
            {
               for( i=0, j=0; i<3; i++ )
               {
                  subnodePtr = (tSubNode *) ct->pPtr(i);
                  subnodePtr->RevertToOldCoords();
               }
            }
            else
            {   
               crossed = FALSE;
               for( i=0; i<3; i++ )
               {
                  cn = (tSubNode *) ct->pPtr(i);
                  if( cn->Meanders() )
                  {
                     cedg = ct->ePtr( (i+1)%3 );
                     spokIter.Reset( cn->getSpokeListNC() );
                     for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                          ce = spokIter.NextP() )
                     {
                        if( Intersect( ce, cedg ) )
                        {
                           if( ct->tPtr(i) == 0 ) //boundary has been crossed
                           {
                              subnodePtr = (tSubNode *) ct->pPtr(i);
                              subnodePtr->RevertToOldCoords();
                           }
                           else
                           {
                              crossed = TRUE;
                              ctop = ct->tPtr(i);
                              xy = cn->getNew2DCoords();
                                //check to make sure the opposite tri is still CCW;
                                //if so, check whether the point has moved into it;
                                //otherwise delete node and re-add it
                              if( NewTriCCW( ctop ) && InNewTri( xy, ctop ) )
                              {
                                   //if node has simply moved into 'opposite' triangle;
                                   //remove opposite tri from ptr list, flipedge,
                                   //add two new tri's to ptr list.
                                 for( rmtri = tpIter.FirstP();
                                      tpIter.ReportNextP() != ctop && !(tpIter.AtEnd());
                                      rmtri = tpIter.NextP() );
                                 if( !(tpIter.AtEnd()) ) //ctop is in tri ptrlist
                                 {
                                    tpListNode = tpIter.NodePtr();
                                    triptrList.removeNext( rmtri, tpListNode );
                                 }                           
                                 nv = ct->nVOp( ctop );
                                 nvopp = ctop->nVOp( ct );
                                 cout << "call FlipEdge from CTEI for edge between nodes "
                                      << ct->pPtr( (nv+1)%3 )->getID() << " and "
                                      << ct->pPtr( (nv+2)%3 )->getID() << endl;
                                 FlipEdge( ct, ctop, nv, nvopp );
                                 rmtri = triIter.LastP();
                                 triptrList.insertAtBack( rmtri );
                                 rmtri = triIter.PrevP();
                                 triptrList.insertAtBack( rmtri );
                              }
                              else
                                    //things have gotten complicated and it's probably
                                    //easier to just delete the node and add it again
                                    //at the new location
                              {
                                 if( LocateTriangle( xy[0], xy[1] ) != 0 )
                                 {
                              
                                      //tempNode = *cn;
                                      //find spoke tri's in tri ptr list and remove them
                                    for( ce = spokIter.FirstP(); !(spokIter.AtEnd());
                                         ce = spokIter.NextP() )
                                    {
                                       rmtri = TriWithEdgePtr( ce );
                                       for( tri = tpIter.FirstP();
                                            tpIter.ReportNextP() != rmtri &&
                                                !(tpIter.AtEnd());
                                            tri = tpIter.NextP() );
                                       if( !(tpIter.AtEnd()) ) //rmtri is in tri ptrlist
                                       {
                                          tpListNode = tpIter.NodePtr();
                                          triptrList.removeNext( rmtri, tpListNode );
                                       }
                                    }
                                      //delete the node;
                                    xyz = cn->getNew3DCoords();
                                    cout << "delete node at " << xyz[0] << ", " << xyz[1]
                                         << ", " << xyz[2] << endl << flush;
                                    tmpNodeList.insertAtBack( *cn );
                                    DeleteNode( cn );
                                 }
                                 else
                                 {
                                    subnodePtr = (tSubNode *) ct->pPtr(i);
                                    subnodePtr->RevertToOldCoords();
                                 }
                              }
                           }
                           break;
                        }
                     }
                  }
                  if( crossed ) break;
               }
            }      
         }
      }
   }
   
   for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
       if ( cn->Meanders() ) cn->UpdateCoords();
   for( cn = tmpIter.FirstP(); !(tmpIter.AtEnd()); cn = tmpIter.NextP() )
   {
      if ( cn->Meanders() ) cn->UpdateCoords();
      cout << "add node at " << cn->getX() << ", " << cn->getY() << ", "
           << cn->getZ() << endl << flush;
      cn->getSpokeListNC().Flush();
      i = AddNode( *cn );
      assert( i == 1 );
   }
   
   for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
   {
      if( ct->tPtr(0) != 0 ) id0 = ct->tPtr(0)->getID();
      else id0 = -1;
      if( ct->tPtr(1) != 0 ) id1 = ct->tPtr(1)->getID();
      else id1 = -1;
      if( ct->tPtr(2) != 0 ) id2 = ct->tPtr(2)->getID();
      else id2 = -1;
      cout << "end of CTEI tri " << ct->getID() << " with nbrs "
           << id0 << ", " << id1 << ", and " << id2 << endl;
   }
   cout << "finished" << endl;
}

               
/*****************************************************************************\
**
**      MoveNodes (formerly PreApply) :
**             The function that deleted some points so that 
**             they won't be a problem in ApplyChanges
**               uses newx newy x y.
**             Separate call to ApplyChanges now unnecessary!!!
**      Data members updated: Grid
**      Called by: 
**      Calls:  ApplyChanges
**        
**      Created: SL
**                                                              
**
\*****************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
MoveNodes()
{
   cout << "MoveNodes()..." << flush << endl;
   tSubNode * cn;  
   tGridListIter< tSubNode > nodIter( nodeList );
   //check for triangles with edges which intersect (an)other edge(s)
   CheckTriEdgeIntersect();
   //update coordinates
     /*for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
       if ( cn->Meanders() ) cn->UpdateCoords();*/
   //resolve any remaining problems after points moved
   CheckLocallyDelaunay();
   cout << "MoveNodes() finished" << endl;
}



/*****************************************************************************\
**
**      DumpEdges(), DumpSpokes(), DumpTriangles(), DumpNodes(): debugging
**         routines which simply write out information pertaining to the grid;
**      DumpNodes() calls DumpSpokes for each node;
**      DumpSpokes() takes a pointer to a node as an argument.
**
**      Created: SL 1/98
**
\*****************************************************************************/
template<class tSubNode>
void tGrid<tSubNode>::
DumpEdges()
{
   tGridListIter< tEdge > edgIter( edgeList );
   tEdge *ce;
   tTriangle *ct;
   int tid;
   cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      ct = TriWithEdgePtr( ce );
      tid = ( ct != 0 ) ? ct->getID() : -1;
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << "; in tri "
           << tid << endl;
   }
}

template<class tSubNode>
void tGrid<tSubNode>::
DumpSpokes( tSubNode *cn )
{
   tEdge *ce;
   tPtrListIter< tEdge > spokIter( cn->getSpokeListNC() );
   cout << "node " << cn->getID() << " with spoke edges " << endl;
   for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
   {
      cout << "   " << ce->getID()
          << " from node " << ce->getOriginPtrNC()->getID()
              << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }
}

   
template<class tSubNode>
void tGrid<tSubNode>::
DumpTriangles()
{
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct, *nt;
   int tid0, tid1, tid2;
   cout << "triangles:" << endl;
   for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
   {
      nt = ct->tPtr(0);
      tid0 = ( nt != 0 ) ? nt->getID() : -1;
      nt = ct->tPtr(1);
      tid1 = ( nt != 0 ) ? nt->getID() : -1;
      nt = ct->tPtr(2);
      tid2 = ( nt != 0 ) ? nt->getID() : -1;
      cout << ct->getID() << " with vertex nodes "
           << ct->pPtr(0)->getID() << ", "
           << ct->pPtr(1)->getID() << ", and "
           << ct->pPtr(2)->getID() << "; edges "
           << ct->ePtr(0)->getID() << ", "
           << ct->ePtr(1)->getID() << ", and "
           << ct->ePtr(2)->getID() << "; nbr triangles "
           << tid0 << ", "
           << tid1 << ", and "
           << tid2 << endl;
   }
}

template<class tSubNode>
void tGrid<tSubNode>::
DumpNodes()
{
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn;
   cout << "nodes: " << endl;
   for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
   {
      cout << " at " << cn->getX() << ", " << cn->getY() << ", " << cn->getZ() << "; ";
      DumpSpokes( cn );
   }
}
