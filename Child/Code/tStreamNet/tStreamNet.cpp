/**************************************************************************\
**
**  tStreamNet.cpp
**
**  Functions for class tStreamNet.
**
**  $Id: tStreamNet.cpp,v 1.2.1.10 1998-02-13 22:48:42 stlancas Exp $
\**************************************************************************/

#include <assert.h>
#include "../errors/errors.h"
#include "tStreamNet.h"


/**************************************************************************\
**
**  Constructors
**
**  1) The default constructor sets parameters to zero, including the
**     grid and storm pointers (the default constructor has no grid nor
**     storm associated with it)
**
**  2) The second constructor asks for references to a grid (i.e., _the_
**     grid), a storm, and an input file). It sets the _gridPtr_ and
**     _stormPtr_, initializes the rainfall rate, and reads in various
**     hydrologic options and parameters from the input file.
**     It then initializes flow directions, drainage areas, etc.
**
\**************************************************************************/
/*tStreamNet::tStreamNet()
        : bedErode()
{
   cout << "tStreamNet()...";
   gridPtr = 0;
   stormPtr = 0;
   flowgen = filllakes = 0;
   rainrate = trans = infilt = 0;
   cout << "finished" << endl;	
}*/

tStreamNet::tStreamNet( tGrid< tLNode > &gridRef, tStorm &storm,
                        tInputFile &infile )
        : bedErode( infile )
{
   cout << "tStreamNet(...)...";
   assert( &gridRef != 0 );
   gridPtr = &gridRef;
   assert( gridPtr != 0 );
   stormPtr = &storm;
   assert( stormPtr != 0 );
   flowgen = infile.ReadItem( flowgen, "FLOWGEN" );
   filllakes = infile.ReadItem( filllakes, "LAKEFILL" );
   if( flowgen == kSaturatedFlow )
   {
      trans = infile.ReadItem( trans, "TRANSMISSIVITY" );
      infilt = 0;
   }
   else
   {
      trans = 0;
      infilt = infile.ReadItem( infilt, "INFILTRATION" );
   }
   int typbnd = infile.ReadItem( typbnd, "TYP_BOUND" );
   if( typbnd == kOppositeSidesOpen )
       inDrArea = infile.ReadItem( inDrArea, "INDRAREA" );
   else
       inDrArea = 0;
   rainrate = stormPtr->GetRainrate();
   CalcSlopes();  // TODO: should be in tGrid
   InitFlowDirs(); // TODO: should all be done in call to updatenet
   FlowDirs();
   SetVoronoiVertices(); //TODO tgrid
   CalcVAreas();  // TODO: tgrid
   MakeFlow();
   cout << "finished" << endl;	
}


tStreamNet::~tStreamNet()
{
   gridPtr = 0;
   stormPtr = 0;
   cout << "~tStreamNet()" << endl;
}

/**************************************************************************\
**
**  Get/Set functions
**
\**************************************************************************/
void tStreamNet::ResetGrid( tGrid< tLNode > &gridRef )
{
   assert( &gridRef != 0 );
   gridPtr = &gridRef;
   assert( gridPtr != 0 );
}

const tGrid< tLNode > *
tStreamNet::getGridPtr() const {return gridPtr;}

tGrid< tLNode > *
tStreamNet::getGridPtrNC() {return gridPtr;}

const tStorm *tStreamNet::getStormPtr() const {return stormPtr;}

tStorm *tStreamNet::getStormPtrNC() {return stormPtr;}

int tStreamNet::getFlowGenOpt() const {return flowgen;}

int tStreamNet::getFillLakesOpt() const {return filllakes;}

double tStreamNet::getRainRate() const {return rainrate;}

double tStreamNet::getTransmissivity() const {return trans;}

double tStreamNet::getInfilt() const {return infilt;}

double tStreamNet::getInDrArea() const {return inDrArea;}

// TODO: the value checks are nice, but will hurt performance. Should
// probably be removed.
void tStreamNet::setFlowGenOpt( int val )
{flowgen = ( val == 0 || val == 1 ) ? val : 0;}

void tStreamNet::setFillLakesOpt( int val )
{filllakes = ( val == 0 || val == 1 ) ? val : 0;}

void tStreamNet::setRainRate( double val ) {rainrate = ( val >= 0 ) ? val : 0;}

void tStreamNet::setTransmissivity( double val )
{trans = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInfilt( double val )
{infilt = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInDrArea( double val )
{inDrArea = ( val >= 0 ) ? val : 0;}



/**************************************************************************\
**
**  UpdateNet
**
**  This function updates the network and flow field by calling various
**  helper functions.
**
**  TODO: move mesh-related routines -- slopes, voronoi areas, etc --
**         to tGrid
**
\**************************************************************************/
void tStreamNet::UpdateNet()
{
   cout << "UpdateNet()...";
   CalcSlopes();          // TODO: should be in tGrid
   SetVoronoiVertices();  // TODO: should be in tGrid
   CalcVAreas();          // TODO: should be in tgrid
   FlowDirs();
   MakeFlow();
   cout << "UpdateNet() finished" << endl;	
}

void tStreamNet::UpdateNet( tStorm &storm )
{
   cout << "UpdateNet(...)...";
   stormPtr = &storm;
   assert( stormPtr != 0 );
   rainrate = stormPtr->GetRainrate();
   CalcSlopes();   // TODO: as above
   SetVoronoiVertices();
   CalcVAreas();
   FlowDirs();
   MakeFlow();
   cout << "UpdatNet(...) finished" << endl;	
}



/****************************************************************************\
**
**  CalcSlopes
**
**  Computes the slope of each edge in the mesh. This is done by computing
**  the slope for the first of each pair of complementary edges, then
**  assigning the negative of that value to the second of the pair.
**
**  Memory allocated: none
**  Notes:
**   - if edge length is zero, it is calculated
**   - complementary edges on the list are assumed to be organized pairwise;
**     that is, edges AB and BA are always together, for example.
**
**  TODO: should be a member of tGrid!
**
\****************************************************************************/
void tStreamNet::CalcSlopes()
{
	assert( gridPtr != 0 );
	tEdge *curedg;
	tGridListIter<tEdge> i( gridPtr->GetEdgeList() );
  double slp, length;

  cout << "CalcSlopes()...";

  // Loop through each pair of edges on the list
	for( curedg = i.FirstP(); !( i.AtEnd() ); curedg = i.NextP() )
	{
     // Make sure edge is valid, and length is nonzero
     assert( curedg > 0 );
     assert( curedg->getLength() > 0 );
     //if( curedg->getLength() == 0 )
     //{
     //   length = curedg->CalcLength();
     //}

     // Compute the slope and assign it to the current edge
     slp = ( curedg->getOrgZ() - curedg->getDestZ() )
         / curedg->getLength();
     curedg->setSlope( slp );

     // Advance to the edge's complement, and assign it -slp
     curedg = i.NextP();
     assert( !( i.AtEnd() ) );
     curedg->setSlope( -slp );
     //curedg->setLength( length );
     assert( curedg->getLength() > 0 );
	}
  cout << "CalcSlopes() finished" << endl;	
}


// NB: move to tGrid
void tStreamNet::CalcVAreas()
{
   cout << "CalcVAreas()..." << endl << flush;
   double area;
   tLNode * curnode, *cn;
   tEdge *ce;
   tArray< double > xy;
   tGridListIter<tLNode> nodIter( gridPtr->GetNodeList() ),
       nI( gridPtr->GetNodeList() );
   tPtrListIter< tEdge > sI;
   
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
         //area = VoronoiArea( curnode );
      curnode->ComputeVoronoiArea();
      
/*      for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
      {
         sI.Reset( cn->getSpokeListNC() );
         cout << "node " << cn->getID() << ": ";
         for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
         {
            xy = ce->getRVtx();
            cout << xy[0] << " " << xy[1] << "; ";
         }
         cout << endl;
      }*/
         /*assert( area > 0 );
         curnode->setVArea( area );
         curnode->setVArea_Rcp( 1.0 / area );*/
   }
   cout << "CalcVAreas() finished" << endl;
}



/****************************************************************************\
**
**  InitFlowDirs
**
**  Initialize flow directions such that each active (non-boundary) node
**  flows to another active node (or open boundary node). This initialization
**  process allows the FlowDirs function to assume that the previous flowedg
**  is always valid in the sense that it doesn't point to a closed boundary
**  (which otherwise could happen when the grid is first read in).
**
**    Written 12/1/97 gt.
**
\****************************************************************************/
#define kMaxSpokes 100
void tStreamNet::InitFlowDirs()
{
   tGridListIter<tLNode> i( gridPtr->GetNodeList() );
   tLNode * curnode;
   tEdge * flowedg;
   int ctr;

   // For every active (non-boundary) node, initialize it to flow to a
   // non-boundary node (ie, along a "flowAllowed" edge)
   curnode = i.FirstP();
   while( i.IsActive() )
   {
      // Start with the node's default edge
      assert( curnode>0 );
      flowedg = curnode->GetEdg();
      assert( flowedg>0 );

      // As long as the current edge is a no-flow edge, advance to the next one
      // counter-clockwise
      ctr = 0;
      while( !flowedg->FlowAllowed() )
      {
         flowedg = flowedg->GetCCWEdg();
         assert( flowedg>0 );
         //Xcout << " edg " << flowedg->getID() << " fa: "
         //X     << flowedg->FlowAllowed() << endl;
         ctr++;
         if( ctr>kMaxSpokes ) // Make sure to prevent endless loops
         {
            cerr << "Mesh error: node " << curnode->getID()
                 << " appears to be surrounded by closed boundary nodes"
                 << endl;
            ReportFatalError( "Bailing out of InitFlowDirs()" );
         }
      }
      curnode->SetFlowEdg( flowedg );
      assert( curnode->GetFlowEdg() > 0 );
      curnode = i.NextP();
   }
   if( inDrArea > 0 )
   {
      while( curnode->getBoundaryFlag() == kOpenBoundary )
      {
        if( curnode->getZ() > 0 )
        {
           flowedg = curnode->GetEdg();
           assert( flowedg>0 );
           ctr = 0;
           while( !flowedg->FlowAllowed() )
           {
              flowedg = flowedg->GetCCWEdg();
              assert( flowedg>0 );
              ctr++;
              if( ctr>kMaxSpokes ) // Make sure to prevent endless loops
              {
                 cerr << "Mesh error: node " << curnode->getID()
                      << " appears to be surrounded by closed boundary nodes"
                      << endl;
                 ReportFatalError( "Bailing out of InitFlowDirs()" );
              }
           }
           curnode->SetFlowEdg( flowedg );
           assert( curnode->GetFlowEdg() > 0 );
        }
        curnode = i.NextP();
      }
   }
}
#undef kMaxSpokes


/****************************************************************************\
**
**  FlowDirs:  Computes flow directions for each node using a 
**             steepest-direction routing algorithm. Flow proceeds along the
**             steepest of the directed edges having their origin at a given
**             point. This edge is stored in .flowedg. 
**
**      Parameters:     none
**      Called by:
**      Assumes:
**       - each node points to a valid flow edge (returned by GetFlowEdg())
**       - each edge has a valid counter-clockwise edge
**      Updated: 12/19/97 SL; 12/30/97 GT
**
\****************************************************************************/
#define kLargeNegative -1000
void tStreamNet::FlowDirs()
{
   tGridListIter<tLNode> i( gridPtr->GetNodeList() );  // Gets nodes from the list
   double slp;        // steepest slope found so far
   tLNode *curnode;  // ptr to the current node
   tEdge * firstedg; // ptr to first edg
   tEdge * curedg;   // pointer to current edge
   tEdge * nbredg;   // steepest neighbouring edge so far

//#if TRACKFNS
  cout << "FlowDirs" << endl;
//#endif

  // Find the connected edge with the steepest slope
  curnode = i.FirstP();
  while( i.IsActive() )  // DO for each non-boundary (active) node
  {
      curnode->SetFloodStatus( kNotFlooded );  // Init flood status indicator
      firstedg =  curnode->GetFlowEdg();
      if( firstedg <= 0 )
          curnode->TellAll();
      assert( firstedg > 0 );
      slp = firstedg->getSlope();
      nbredg = firstedg;
      curedg = firstedg->GetCCWEdg();
      
      // Check each of the various "spokes", stopping when we've gotten
      // back to the beginning
      while( curedg!=firstedg )
      {
         assert( curedg > 0 );
         if ( curedg->getSlope() > slp && curedg->FlowAllowed() )
         {
            slp = curedg->getSlope();
            nbredg = curedg;
         }
         curedg = curedg->GetCCWEdg();
      }
      
      curnode->SetFlowEdg( nbredg );
      curnode->SetFloodStatus( ( slp>0 ) ? kNotFlooded : kSink );  // (NB: opt branch pred?)
      //cout << "Node " << curnode->getID() << " flows to "
      //     << curnode->GetDownstrmNbr()->getID() << endl;
      curnode = i.NextP();

  }
  if( inDrArea > 0 )
  {
     while( curnode->getBoundaryFlag() == kOpenBoundary )
     {
        if( curnode->getZ() > 0 )
        {
           firstedg =  curnode->GetFlowEdg();
           slp = firstedg->getSlope();
           nbredg = firstedg;
           curedg = firstedg->GetCCWEdg();
           
             // Check each of the various "spokes", stopping when we've gotten
             // back to the beginning
           while( curedg!=firstedg )
           {
              assert( curedg > 0 );
              if ( curedg->getSlope() > slp && curedg->FlowAllowed() )
              {
                 slp = curedg->getSlope();
                 nbredg = curedg;
              }
              curedg = curedg->GetCCWEdg();
           }
           curnode->SetFlowEdg( nbredg );
        }
        curnode = i.NextP();
     }
  }
  cout << "FlowDirs() finished" << endl << flush;
  
}
#undef kLargeNegative



/*void tStreamNet::FlowDirs()
{
//#if TRACKFNS
   cout << "FlowDirs()...";
//#endif
   tGridListIter<tLNode> i( gridPtr->GetNodeList() );  // Gets nodes from the list
   tPtrListIter< tEdge > spokIter;
   double slp;        // steepest slope found so far
   tLNode *curnode;  // ptr to the current node
   tEdge * curedg;   // pointer to current edge
   tEdge * nbredg;   // steepest neighbouring edge so far
     // Find the connected edge with the steepest slope
   curnode = i.FirstP();
   while( i.IsActive() )  // DO for each non-boundary (active) node
   {
      spokIter.Reset( curnode->getSpokeListNC() );
      curnode->SetFloodStatus( kNotFlooded );  // Init flood status indicator
      curedg = spokIter.DatPtr();
      assert( curedg != 0 );
      slp = curedg->getSlope();
      nbredg = curedg;
        // Check each of the various "spokes", stopping when we've gotten
        // back to the beginning
      while( spokIter.NextIsNotFirst() )
      {
         curedg = spokIter.NextP();
         assert( curedg > 0 );
         if ( curedg->getSlope() > slp && curedg->FlowAllowed() )
         {
            slp = curedg->getSlope();
            nbredg = curedg;
         }
      }
      curnode->SetFlowEdg( nbredg );
      curnode->SetFloodStatus( ( slp>0 ) ? kNotFlooded : kSink );
// (NB: opt branch pred?)
        //cout << "Node " << curnode->getID() << " flows to node "
        //   << curnode->GetDownstrmNbr()->getID() << endl;
      curnode = i.NextP();
   }
   cout << "finished" << endl;
}*/


/*****************************************************************************\
**
**  tStreamNet::DrainAreaVoronoi
**
**  Computes drainage areas for nodes on the grid by 
**  getting the Voronoi areas of the nodesand routing
**  flow down the flowedge of the node.
**
\*****************************************************************************/
void tStreamNet::DrainAreaVoronoi()
{
//#if TRACKFNS
   cout << "DrainAreaVoronoi()..." << endl << flush;
//#endif
   double slope, maxslope;
   tLNode * curnode, *innode;
   tGridListIter<tLNode> nodIter( gridPtr->GetNodeList() );
   
   // Reset drainage areas to zero
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
       curnode->SetDrArea( 0 );
   if( inDrArea > 0 )
   {
      maxslope = -10000;
      for( ; curnode->getBoundaryFlag() == kOpenBoundary; curnode = nodIter.NextP() )
      {
         if( curnode->getZ() > 0 )
         {
            slope = curnode->GetSlope();
            if( slope > maxslope )
            {
               maxslope = slope;
               innode = curnode;
            }
            innode->GetDownstrmNbr()->SetDrArea( inDrArea );
         }
      }
   }
   
   // send voronoi area for each node to the node at the other end of the 
   // flowedge and downstream 
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      RouteFlowArea( curnode, curnode->getVArea() );
   }
   
   cout << "DrainAreaVoronoi() finished" << endl << flush;
}


//MOVED TO tNODE
/******************************************************************************\
**
**  VoronoiArea: modified such that it takes the area of the
**                       closest (counterclockwise) intersection of 
**                       perpendicular bisectors; in some cases,
**                       if the voronoi vertices are used arbitrarily,
**                       errors result when perp. bisectors of non-
**                       adjacent edges form a closer intersection.
**                       Also sets edges' vedglen, in some cases to zero..
**                       
**         TODO: make member of tNode              
\******************************************************************************/
/*
double tStreamNet::VoronoiArea( tLNode * centre )
{
//#if DEBUG
   //cout<<"VoronoiArea(...)...";
//#endif
   int cw;
   double area = 0;
   double a, b, c, dx, dy, dx0, dx1, dy0, dy1, dx2, dy2;
   double x, y, x0, y0, x1, y1, x2, y2, m1, m2;
   tEdge * ce, *edgptr;
   tPtrList< tEdge > vedgList;
   tPtrListIter< tEdge > //XspokIter( centre->getSpokeListNC() ),
       vtxIter( vedgList );
   tArray< double > xy, xyn, xynn, xynnn, xy1, xy2;


   // Create a duplicate list of edges.... why?
   ce = centre->GetEdg();
   do
   {
      assert( ce>0 );
      vedgList.insertAtBack( ce );
      ce = ce->GetCCWEdg();
   } while( ce != centre->GetEdg() );
   vedgList.makeCircular();

   // Check boundary status: Voronoi area only defined for non-boundary nodes
   if( centre->getBoundaryFlag() == kNonBoundary )
   {
      centre->CalcSpokeVEdgLengths();
      cw = TRUE;
      do
      {
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
                  x0 = centre->getX();
                  y0 = centre->getY();
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
                     x = (y2-m2*x2-y1+m1*x1) / (m1-m2);
                     y = m1*(x-x1)+y1;
                  }
                  else
                  {
                     m1 = dy1/dx1;
                     m2 = dy2/dx2;
                     y=(m1*y1+x1-m2*y2-x2)/(m1-m2);
                     x= -y*m1+m1*y1+x1;
                  }
                  assert( vtxIter.Prev() );
                  xyn[0] = x;
                  xyn[1] = y;
                  dx = xy[0] - x;
                  dy = xy[1] - y;
                  cout << "reset vedglen and rvtx for edge "
                       << vtxIter.DatPtr()->getID() << " to len "
                       << sqrt( dx*dx + dy*dy )
                       << ", x, y, " << xyn[0] << ", " << xyn[1] << endl;
                  
                  vtxIter.DatPtr()->setVEdgLen( sqrt( dx*dx + dy*dy ) );
                  vtxIter.DatPtr()->setRVtx( xyn );
                  cout << "reset vedglen and rvtx for edge "
                       << vtxIter.ReportNextP()->getID()
                       << " to len 0.0, x, y, " << xynnn[0] << ", "
                       << xynnn[1] << endl;
                  vtxIter.ReportNextP()->setVEdgLen(0.0);
                  vtxIter.ReportNextP()->setRVtx( xynnn );
                  assert( vedgList.removeNext( edgptr, vtxIter.NodePtr() ) );
               }
            }
            assert( vtxIter.Get( ce->getID() ) );
         }
      } while( cw );
      
      for( ce = vtxIter.FirstP();
           vtxIter.NodePtr()->getNext()->getNext() != vedgList.getFirst();
           ce = vtxIter.NextP() )
      {
         xy = ce->getRVtx();
         xyn = vtxIter.NextP()->getRVtx();
         xynn = vtxIter.ReportNextP()->getRVtx();
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
         area += sqrt( 4*a*a*b*b - (c*c - (b*b + a*a))*(c*c - (b*b + a*a))) / 4;
         assert( vtxIter.Prev() );
      }
   }
   
   else //boundary node
       for( ce = vtxIter.FirstP(); !( vtxIter.AtEnd() ); ce = vtxIter.NextP() )
           ce->setVEdgLen(0.0);
   cout<<"Voronoi_area: node id "<<centre->getID()<<", varea "
       <<area<<endl;
     //cout<<"V_area done"<<endl<<flush;
   //cout << "finished: voronoi area = " << area << endl;
   return area;
}*/


/*****************************************************************************\
**
**  tNode::RouteFlowArea
**
**  Starting with the current node, this routine increments 
**  the drainage area of the node and each node downstream by _addedArea_. 
**  (GT 2/97) 
**
\*****************************************************************************/
void tStreamNet::RouteFlowArea( tLNode *curnode, double addedArea )
{
   //cout << "RouteFlowArea()..." << endl << flush;
//#if DEBUG
   int niterations=0;  // Safety feature: prevents endless loops
//#endif

   // As long as the current node is neither a boundary nor a sink, add
   // _addedArea_ to its total drainage area and advance to the next node
   // downstream
   while( !(curnode->getBoundaryFlag()) && (curnode->GetFloodStatus()!=kSink) )
   {
      curnode->AddDrArea( addedArea );
      curnode = curnode->GetDownstrmNbr();
//#if DEBUG
      niterations++;
      if( niterations>9990 )
          cout << "Node " << curnode->getID()
               <<" boundary status "<< curnode->getBoundaryFlag()
               << " flood status " << curnode->GetFloodStatus() 
               << endl << flush;
      assert( niterations < 10000 );
//#endif
   }
   //cout << "RouteFlowArea() finished" << endl << flush;
}


/*****************************************************************************\
**
**  SetVoronoiVertices
**
**  Each Delaunay triangle is associated with an intersection between
**  three Voronoi cells, called a Voronoi vertex. These Voronoi vertices
**  are used in computing the area of each Voronoi cell. The Voronoi
**  vertex associated with each triangle is the circumcenter of the
**  triangle. This routine finds the Voronoi vertex associated with
**  each triangle by finding the triangle's circumcenter. 
**
**  The vertex coordinates are stored in the three counter-clockwise-
**  oriented tEdge objects associated with each triangle. This is a 
**  space-for-time tradeoff: the coordinates could be stored in the triangles, 
**  saving redundancy (3 copies of each point are stored here), but in
**  that case each tEdge would have to point back to a triangle, and an
**  additional level of indirection would be needed in accessing the Voronoi
**  vertices associated with a particular node.
**
**  Note that computation of the Voronoi vertices can be prone to numerical
**  errors, leading to inconsistent Voronoi polygons with cells that
**  overlap or have loops. See note under tNode::ComputeVoronoiArea().
**
**    Assumes: correct triangulation with valid edge pointers in each tri.
**    Data members modified: none
**    Other objects modified: Voronoi vertices set for each tEdge
**    Modifications:
**     - reverted to earlier triangle-based computation, from an edge-based
**       computation that takes 3x as long because NE = 3NT. In so doing,
**       the definition of the Voronoi vertex stored in a tEdge is changed
**       to "left-hand", meaning the V. vertex associated with the edge's
**       lefthand triangle (the vertex itself may or may not lie to the left
**       of the edge). 1/98 GT
**     - also moved circumcenter computation into a tTriangle mbr fn.
**
\*****************************************************************************/
void tStreamNet::SetVoronoiVertices()
{
   //double x, y, x1, y1, x2, y2, dx1, dy1, dx2, dy2, m1, m2;
   //tArray< double > xyo, xyd1, xyd2, xy(2);
   cout << "SetVoronoiVertices()..." << endl;
   tArray< double > xy;
   tListIter< tTriangle > triIter( gridPtr->GetTriList() );
   tTriangle * ct;

   // Find the Voronoi vertex associated with each Delaunay triangle
   for( ct = triIter.FirstP(); !(triIter.AtEnd()); ct = triIter.NextP() )
   {
      xy = ct->FindCircumcenter();    
      //cout << "SetVoronoiVertices(): " << xy[0] << " " << xy[1];
      // Assign the Voronoi point as the left-hand point of the three edges 
      // associated with the current triangle
      ct->ePtr(0)->setRVtx( xy );
      ct->ePtr(1)->setRVtx( xy );
      ct->ePtr(2)->setRVtx( xy );
      //xy = ct->ePtr(0)->getRVtx();
      //cout << "SetVoronoiVertices(): " << xy[0] << " " << xy[1];
   }
   cout << "SetVoronoiVertices() finished" << endl;
}


/*****************************************************************************\
**
**      MakeFlow : flow routing functions 
**
**      Data members updated: 
**      Called by: 
**      Calls: FillLakes, DrainAreaVoronoi, FlowSaturated, FlowUniform
**
**      Created:  YC
**      Added:   YC
**      Modified:  GT made it a mbr func of tGrid 7/97
**         updated 12/19/97 SL
**
\*****************************************************************************/
void tStreamNet::MakeFlow()
{
   cout << "MakeFlow()..." << endl << flush;
   if( filllakes ) FillLakes();
   DrainAreaVoronoi();
   if( flowgen == kSaturatedFlow ) FlowSaturated();
   else FlowUniform();
   cout << "MakeFlow() finished" << endl;
}


/*****************************************************************************\
**
**  FlowUniform
**
**  Computes discharge as the product of precip and 
**                       drainage area. 3/97 GT
**
**  Parameters:  precip -- precipitation rate
**  Called by:  main
**  Modifications:
**    - added infiltration parameter (default zero) (8/97 GT)
**    - updated 12/19/97 SL
**
\*****************************************************************************/
void tStreamNet::FlowUniform()
{
   cout << "FlowUniform..." << endl << flush;
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   tLNode *curnode;
   double runoff = rainrate - infilt;
   double discharge;
   
   if( runoff<0 ) runoff = 0;  // Make sure runoff isn't negative
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      discharge = curnode->getDrArea() * runoff;
      curnode->setDischarge( discharge );
   }
   cout << "FlowUniform finished" << endl;
}


/*****************************************************************************\
**
**  FlowSaturated
**
**  Computes surface runoff using Topmodel concept. Steady-state 
**  subsurface flow capacity is transmissivity times slope. Total
**  steady-state runoff is rainfall rate times drainage area. Surface
**  runoff is total minus subsurface. (5/97 GT)
**
**  Parameters:  pr -- parameter block
**               precip -- precipitation rate
**  Called by:  main
**  Modifications:
**    - updated 12/19/97 SL
**
\*****************************************************************************/
void tStreamNet::FlowSaturated()
{
   cout << "FlowSaturated...";
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   tLNode *curnode;
   double discharge;
   
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      discharge = curnode->getDrArea() * rainrate -
          curnode->GetFlowEdg()->getSlope() * trans;
      discharge *= ( discharge > 0 );
      curnode->setDischarge( discharge );
   }
   cout << "finished" << endl;
}


/*****************************************************************************\
**
**  FillLakes 
**
**              Finds drainage for closed depressions. The algorithm assumes
**              that sinks (nodes that are lower than any of their neighbors)
**              have already been identified during the flow directions
**              procedure. For each sink, the algorithm creates a list of
**              nodes in the current lake, which initially is just the sink
**              itself. It then iteratively looks for the lowest node on the
**              perimeter of the current lake. That node is checked to see 
**              whether it can be an outlet, meaning that one of its
**              neighbors is both lower than itself and is not already
**              flooded (or is an open boundary). If the low node is not an 
**              outlet, it is added to the current lake and the process 
**              is repeated. If it is an outlet, then all of the nodes on the 
**              current-lake list are identified draining it. The list is then
**              cleared, and the next sink is processed. If during the search
**              for the lowest node on the perimeter a flooded node is found
**              that isn't already part of the current lake (i.e., it was
**              flagged as a lake node when a previous sink was processed),
**              then it is simply added to the current-lake list --- in other
**              words, the "new" lake absorbs any "old" ones that are
**              encountered.
**
**  Parameters: (none)
**  Calls: FindLakeNodeOutlet
**  Called by: main
**  Created: 6/97 GT
**  Modified: fixed memory leak on deletion of lakenodes 8/5/97 GT
**  Updated: 12/19/97 SL
**
\*****************************************************************************/
void tStreamNet::FillLakes()
{
//#if TRACKFNS
   cout << "FillLakes()..." << endl << flush;
//#endif
   tLNode *cn,    // Node on list: if a sink, then process
       *thenode,           // Node on lake perimeter
       *lowestNode,        // Lowest node on perimeter found so far
       *cln, // current lake node
       *node; //placeholder
   tGridListIter< tLNode > nodIter( gridPtr->GetNodeList() );
   tPtrList< tLNode > lakeList;
   tPtrListIter< tLNode > lakeIter( lakeList );
   //tPtrListIter< tEdge > spokIter;
   tEdge *ce;           // Pointer to an edge
   double lowestElev;    // Lowest elevation found so far on lake perimeter
   int done;          // Flag indicating whether outlet has been found
   
   // Check each active node to see whether it is a sink
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      if( cn->GetFloodStatus() == kSink )
      {
         // Create a new lake-list, initially containing just the sink node.
         lakeList.insertAtBack( cn );
         cn->SetFloodStatus( kCurrentLake );
         
         // Iteratively search for an outlet along the perimeter of the lake
         done = FALSE;
         do
         {
            lowestNode = lakeIter.FirstP();
            lowestElev = kVeryHigh; // Initialize lowest elev to a very high val.
            // Check the neighbors of every node on the lake-list
            for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
                 cln = lakeIter.NextP() )
            {
               // Check all the neighbors of the node
                      //XspokIter.Reset( cln->getSpokeListNC() );
                   //Xfor( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                   //Xce = spokIter.NextP() )
               ce = cln->GetEdg();
               do
               {
                  thenode = (tLNode *) ce->getDestinationPtrNC();
                  // Is it a potential outlet (ie, not flooded and not
                  // a boundary)?
                  if( thenode->GetFloodStatus() == kNotFlooded
                      && ce->FlowAllowed() )
                  {
                     // Is it lower than the lowest found so far?
                     if( thenode->getZ() < lowestElev )
                     {
                        lowestNode = thenode;
                        lowestElev = thenode->getZ();
                     }
                  }
                  // If it's a previous lake node or a sink, add it to the list
                  else if( thenode->GetFloodStatus() == kFlooded ||
                           thenode->GetFloodStatus() == kSink )
                  {
                     lakeList.insertAtBack( thenode );
                     thenode->SetFloodStatus( kCurrentLake );
                  }
               } while( ( ce=ce->GetCCWEdg() ) != cln->GetEdg() );// END spokes
               
            } /* END lakeList */
            
            // Now we've found the lowest point on the perimeter. Now test
            // to see whether it's an outlet. If it's an open boundary, it's
            // an outlet...
            if( lowestNode->getBoundaryFlag() == kOpenBoundary ) done = TRUE;
            else // ...it's also an outlet if it can drain to a "dry" location.
            {
               // Can lowestNode drain to a non-flooded location?
               if( FindLakeNodeOutlet( lowestNode ) ) done = TRUE;
               // no, it can't, so add it to the list and continue:
               else
               {
                  lakeList.insertAtBack( lowestNode );
                  lowestNode->SetFloodStatus( kCurrentLake );
               }
            }
            assert( lakeList.getSize() < gridPtr->GetNodeList()->getActiveSize() );
         } while( !done );

         // Now we've found an outlet for the current lake.
         // This next bit of code assigns a flowsTo for each node so there's
         // a complete flow path toward the lake's outlet. This isn't strictly
         // necessary --- the nodes could all point directly to the outlet,
         // skipping anything in between --- but it prevents potential problems
         // in ordering the list by network order. This also works by pointing
         // each node toward the first neighboring node they happen to find
         // that has been flagged as having its flow direction resolved. 
         // Initially, the low node is thus flagged, and the algorithm repeats
         // until all the lake nodes are flagged as having a flow direction.
         // The algorithm isn't unique---there are many paths that could be
         // taken; this simply finds the most convenient one.
         lowestNode->SetFloodStatus( kOutletFlag );
//#if LKDEBUG
//         printf("Outlet found to node %d\n",lowestNode->getID() );
//#endif
         // Test for error in mesh: if the lowestNode is a closed boundary, it
         // means no outlet can be found.
         do
         {
            done = TRUE;  // assume done until proven otherwise
            for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
                 cln = lakeIter.NextP() )
            {
               if( cln->GetFloodStatus() != kOutletFlag )
               {
                  done = FALSE;
                  /*XspokIter.Reset( cln->getSpokeListNC() );
                  for( ce = spokIter.FirstP();
                       cln->GetFloodStatus() != kOutletFlag &&
                           !( spokIter.AtEnd() );
                       ce = spokIter.NextP() )*/
                  // Check each neighbor
                  ce = cln->GetEdg();
                  do
                  {
                     node = (tLNode *) ce->getDestinationPtrNC();
                     if( node->GetFloodStatus() == kOutletFlag )
                     {     // found one!  
                        cln->SetFloodStatus( kOutletPreFlag );
                        cln->SetFlowEdg( ce );
                        //cout << "Node " << cln->getID() << " flows to "
                        //     << cln->GetDownstrmNbr()->getID() << endl;
                        
                     }
                  } while( cln->GetFloodStatus() != kOutletFlag
                           && ( ce=ce->GetCCWEdg() ) != cln->GetEdg() );
               } // END if node not flagged as outlet
            } // END for each lake node

            // Now flag all the "preflagged" lake nodes as outlets
            for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
                 cln = lakeIter.NextP() )
                if( cln->GetFloodStatus() == kOutletPreFlag )
                    cln->SetFloodStatus( kOutletFlag );
            
         } while( !done );
         lowestNode->SetFloodStatus( kNotFlooded );
         
         // Finally, flag all of the 
         // nodes in it as "kFlooded" and clear the list so we can move on to
         // the next sink. (Fixed mem leak here 8/5/97 GT).
         for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
              cln = lakeIter.NextP() )
             cln->SetFloodStatus( kFlooded );
         lakeList.Flush();
      } /* END if Sink */
   } /* END Active Nodes */
   cout << "FillLakes() finished" << endl << flush;
   
} // end of tStreamNet::FillLakes


/*****************************************************************************\
**
**  FindLakeNodeOutlet
**
**  This function is part of the lake-filling algorithm. It checks to see 
**  whether there is a valid outlet for the current node, and if so it 
**  assigns that outlet. An "outlet" essentially means a downhill neighbor 
**  that isn't already flooded to the level of the current node. The function 
**  performs basically the same operation as FlowDirs, but with stricter 
**  criteria. The criteria for a valid outlet are:
**
**  (1) It must be lower than the current node (slope > 0)
**  (2) It must not be part of the current lake (a lake can't outlet to itself)
**  (3) It must not be a closed boundary (_flowAllowed_ must be TRUE)
**  (4) The outlet node must drain to a location lower than the current node 
**      (but not necessarily lower than itself; it could be a separate lake)
**      OR the outlet node is an open boundary.
**
**  Returns: TRUE if a valid outlet is found, FALSE otherwise
**  Calls: (none)
**  Called by: FillLakes
**  Created: 6/97 GT
**  Updated: 12/19/97 SL; 1/15/98 gt bug fix (open boundary condition)
**
\*****************************************************************************/
int tStreamNet::FindLakeNodeOutlet( tLNode *node )
{
   double maxslp = 0;  // Maximum slope found so far
   tEdge * ce;        // Current edge
   //XtPtrListIter< tEdge > spokIter( node->getSpokeListNC() );
   tLNode *dn;        // Potential outlet
   
   // Check all the neighbors
      //Xfor( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
   ce = node->GetEdg();
   do
   {
        // If it passes this test, it's a valid outlet
      dn = (tLNode *) ce->getDestinationPtrNC();
      assert( dn>0 );
      if( ce->getSlope() > maxslp &&
          dn->GetFloodStatus() != kCurrentLake &&
          ce->FlowAllowed() &&
          ( dn->getBoundaryFlag()==kOpenBoundary ||
            dn->GetDownstrmNbr()->getZ() < node->getZ() ) )
          //dn->GetDownstrmNbr()->getZ() < node->getZ() )
      {
         maxslp = ce->getSlope();
         node->SetFlowEdg( ce );
         //cout << "Node " << node->getID() << " flows to "
         //     << node->GetDownstrmNbr()->getID() << endl;
      }
   } while( ( ce=ce->GetCCWEdg() ) != node->GetEdg() );
   
   return( maxslp > 0 );
}



/*****************************************************************************\
**
**  SortNodesByNetOrder:
**
**  This function sorts the list of nodes according to their order in the
**  network (upstream to downstream), in preparation for computing erosion &
**  deposition. (Note that this is only necessary when the sediment output
**  from a given node depends on the input, e.g. for mixed bedrock-alluvial
**  mode in which case the channel type [br or alluvial] depends on the
**  difference between sediment influx and carrying capacity). The sorting
**  algorithm is based on the "cascade" algorithm of Braun and Sambridge
**  (Basin Research, 1997, vol. 9, p. 27).
**    The algorithm works by initially assigning a tracer (like a packet
**  of water) to each node. At each iteration, a tracer from each node is
**  sent downstream. Any nodes that have zero tracers left are moved to the
**  bottom of the list (a FIFO stack), so that for example the very first
**  node moved will be the first node on the list when the sorting is
**  completed. The process continues until no unsorted nodes remain.
**
**  Modifications:
**   - adapted from previous CHILD code by GT, 12/97
**
\*****************************************************************************/
void tStreamNet::SortNodesByNetOrder()
{
   int nThisPass;                      // Number moved in current iteration
   int i;
   int done;
   tLNode * cn;
   tGridList<tLNode> nodeList = gridPtr->GetNodeList();
   int nUnsortedNodes = nodeList.getActiveSize();  // Number not yet sorted
   tGridListIter<tLNode> listIter( nodeList );
   
   //test
   cout << "BEFORE: " << endl;
   for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() ) 
       cout << cn->getID() << endl;
   
#if TRACKFNS
   cout << "SortNodesByNetOrder" << endl;
#endif
   
   // Assign initial tracers: use "qs" field, which contains garbage at
   // this stage.
   for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() ) 
      cn->ActivateSortTracer();

   // Iterate: move tracers downstream and sort until no nodes with tracers
   // are left.
   do
   {
      // Send tracers downstream
      cn = listIter.FirstP();
      for( i=1; i<=nUnsortedNodes; i++ )
      {
         cn->MoveSortTracerDownstream();
         cn = listIter.NextP();
      }

      // Scan for any nodes that have no tracers, and move them to the bottom
      // of the list.
      tListNode< tLNode > * nodeToMove;
      nThisPass = 0;
      done = TRUE;
      cn = listIter.FirstP();
      for( i=1; i<=nUnsortedNodes; i++ )
      {
         if( cn->NoMoreTracers() )  // If no tracers, move to bottom of list
         {
            nodeToMove = listIter.NodePtr();
            cn = listIter.NextP();
            nodeList.moveToBack( nodeToMove );
            nThisPass++;
         }
         else
         {
            cn = listIter.NextP();
            done = FALSE;
         }
      }
      
      nUnsortedNodes -= nThisPass;

    } while( !done );

  /*  cout << "AFTER: " << endl;
  //test
cur = firstnode;
  for (i =1; i<=nActiveNodes; i++ ){
  cout << cur->id << " " << cur->drarea << " " << cur->flood 
       << " " << cur->flowsTo->id << endl;
  cur=cur->next;}
  */  
 
}


/*****************************************************************************\
**
**  ErodeDetachLim
**
**  Solves for erosion and deposition during a time interval dtg, for the
**  case in which any sediment detached from the landscape is assumed to
**  be carried away. This might be appropriate, for example, in a bedrock
**  channel system in which deposition is negligible. This case is handled
**  separately from the more general case in which transport capacity and
**  deposition are taken into account, because the numerical solutions
**  to detachment-limited erosion equations tend to be considerably more
**  stable than the general case.
**
**  The function will solve the erosion equation(s) over one or more time
**  intervals within the total "global" time period dtg. First, the
**  maximum time step size dt is computed. Then erosion is computed
**  iteratively in increments of dt until the total time dtg has been used.
**
\*****************************************************************************/
void tStreamNet::ErodeDetachLim( double dtg )
{
   double dt,
       dtmax = 1000000.0; // time increment: initialize to arbitrary large val
   int i;
   tLNode * cn;
   int nActNodes = gridPtr->GetNodeList()->getActiveSize();
   tGridListIter<tLNode> ni( gridPtr->GetNodeList() );
   tArray<double> dz( nActNodes ); // Erosion depth @ each node

   // First, estimate maximum stable/accurate time step size
   for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
       if( (dt=bedErode.SetTimeStep( cn )) < dtmax )
           dtmax = dt;
   dt = dtmax;
   cout << "In ErodeDetachLim, time step set to " << dtmax;

   // Iterate until total time dtg has been consumed
   do
   {
      if( dt > dtg ) dt = dtg; // dt shouldn't exceed remaining time
        // Compute erosion depths
      for( i=0, cn=ni.FirstP(); i<nActNodes; i++, cn=ni.NextP() )
          dz[i] = -bedErode.DetachCapacity( cn, dt );
        // Apply erosion depths
      for( i=0, cn=ni.FirstP(); i<nActNodes; i++, cn=ni.NextP() )
          cn->EroDep( dz[i] );
      dtg -= dt;  // decrease remaining time
   } while( dtg>0 );
   
}

   
