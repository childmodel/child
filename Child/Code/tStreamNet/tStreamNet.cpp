/**************************************************************************\
**
**  tStreamNet.cpp
**
**  Functions for class tStreamNet and related class tInlet.
**
**  $Id: tStreamNet.cpp,v 1.2.1.42 1998-07-15 22:23:17 gtucker Exp $
\**************************************************************************/

#include <assert.h>
#include "../errors/errors.h"
#include "tStreamNet.h"

/*****************************************************************************\
**
**
**      DistanceToLine: given x,y coords, finds distance to the line 
**              defined by given a, b, and c (ax + by + c = 0)
**      Global function.
**      Data members updated: 
**      Called by: 
**      Calls:  
**        
**
\*****************************************************************************/
double DistanceToLine( double x2, double y2, double a, double b, double c )
{
   double f, g, h, x, y, d;

   f = -b;
   g = a;
   h = b * x2 - a * y2;
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
   d = sqrt( (x2 - x) * (x2 - x) + (y2 - y) * (y2 - y) );
   return d;
}

/*****************************************************************************\
**
**
**      DistanceToLine: given x,y coords, finds distance to the line 
**              formed by points  p0->(x, y) and p1->(x, y)
**      Global function.
**      Data members updated: 
**      Called by: 
**      Calls:  
**        
**
\*****************************************************************************/
double DistanceToLine( double x2, double y2, tNode *p0, tNode *p1 )
{
   double a, b, c, f, g, h, x0, y0, x1, y1, x, y, d;

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


/**************************************************************************\
**  Functions for class tStreamNet here:
\**************************************************************************/

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
    : inlet( &gridRef, infile )
{
   cout << "tStreamNet(...)...";
   assert( &gridRef != 0 );
   gridPtr = &gridRef;
   assert( gridPtr != 0 );
   stormPtr = &storm;
   assert( stormPtr != 0 );

   // Read option for runoff generation and get relevant parameters
   flowgen = infile.ReadItem( flowgen, "FLOWGEN" );
   filllakes = infile.ReadItem( filllakes, "LAKEFILL" );
   infilt = trans = 0;
   if( flowgen == kSaturatedFlow1 || flowgen==kSaturatedFlow2 )
      trans = infile.ReadItem( trans, "TRANSMISSIVITY" );
   if( flowgen==kSaturatedFlow2 || flowgen==kConstSoilStore
       || flowgen==kHortonian )
      infilt = infile.ReadItem( infilt, "INFILTRATION" );
   else infilt = 0.0;
   if( flowgen == kConstSoilStore )
       soilStore = infile.ReadItem( soilStore, "SOILSTORE" );
   else soilStore = 0.0;
   
   rainrate = stormPtr->getRainrate();
   int itMeanders = infile.ReadItem( itMeanders, "OPTMNDR" );
   if( itMeanders )
       mndrDirChngProb = infile.ReadItem( mndrDirChngProb, "CHNGPROB" );
   else mndrDirChngProb = 1.0;
   
   optrainvar = infile.ReadItem( optrainvar, "OPTVAR" );
   kwds = infile.ReadItem( kwds, "HYDR_WID_COEFF_DS" );
   //cout << "kwds: " << kwds << endl;
   assert( kwds > 0 );
   ewds = infile.ReadItem( ewds, "HYDR_WID_EXP_DS" );
   //cout << "ewds: " << ewds << endl;
   ewstn = infile.ReadItem( ewstn, "HYDR_WID_EXP_STN" );
   //cout << "ewstn: " << ewstn << endl;
   knds = infile.ReadItem( knds, "HYDR_ROUGH_COEFF_DS" );
   //cout << "knds: " << knds << endl;
   assert( knds > 0 );
   ends = infile.ReadItem( ends, "HYDR_ROUGH_EXP_DS" );
   //cout << "ends: " << ends << endl;
   enstn = infile.ReadItem( enstn, "HYDR_ROUGH_EXP_STN" );
   //cout << "enstn: " << enstn << endl;
   klambda = infile.ReadItem( klambda, "BANK_ROUGH_COEFF" );
   elambda = infile.ReadItem( elambda, "BANK_ROUGH_EXP" );
   
   CalcSlopes();  // TODO: should be in tGrid
   InitFlowDirs(); // TODO: should all be done in call to updatenet
   FlowDirs();
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
**  get/set functions
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

inline int tStreamNet::getFlowGenOpt() const {return flowgen;}

int tStreamNet::getFillLakesOpt() const {return filllakes;}

double tStreamNet::getRainRate() const {return rainrate;}

double tStreamNet::getTransmissivity() const {return trans;}

double tStreamNet::getInfilt() const {return infilt;}

double tStreamNet::getSoilStore() const {return soilStore;}

double tStreamNet::getInDrArea() const {return inlet.inDrArea;}

double tStreamNet::getInSedLoad() const {return inlet.inSedLoad;}

tArray< double >
tStreamNet::getInSedLoadm( ) const
{
   return inlet.inSedLoadm;
}

tLNode *tStreamNet::getInletNodePtr() const {return inlet.innode;}
tLNode *tStreamNet::getInletNodePtrNC() {return inlet.innode;}

double tStreamNet::getMndrDirChngProb() const {return mndrDirChngProb;}

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
{inlet.inDrArea = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInSedLoad( double val )
{inlet.inSedLoad = ( val >= 0 ) ? val : 0;}

void tStreamNet::setInSedLoadm( int i, double val )
{
   if( i > inlet.inSedLoadm.getSize()-1 )
       ReportFatalError("Tried to index a size that doesn't exist in tstreamnet");
   inlet.inSedLoadm[i]=val;
}

void tStreamNet::setInletNodePtr( tLNode *Ptr )
{inlet.innode = ( Ptr > 0 ) ? Ptr : 0;}

void tStreamNet::setMndrDirChngProb( double val )
{mndrDirChngProb = ( val >= 0.0 && val <= 1.0 ) ? val : 1.0;}



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
   //cout << "UpdateNet()...";
   CalcSlopes();          // TODO: should be in tGrid
   //XsetVoronoiVertices();  // TODO: should be in tGrid
   //XCalcVAreas();          // TODO: should be in tgrid
   //InitFlowDirs();
   FlowDirs();
   MakeFlow();
   //CheckNetConsistency();
   //cout << "UpdateNet() finished" << endl;	
}

void tStreamNet::UpdateNet( tStorm &storm )
{
   //cout << "UpdateNet(...)...";
   stormPtr = &storm;
   assert( stormPtr != 0 );
   rainrate = stormPtr->getRainrate();
   CalcSlopes();   // TODO: as above
   //XsetVoronoiVertices();
   //XCalcVAreas();
   //InitFlowDirs(); // necessary?
   FlowDirs();
   MakeFlow();
   CheckNetConsistency();
   //cout << "UpdatNet(...) finished" << endl;	
}

#define kLargeNumber 1000000
void tStreamNet::CheckNetConsistency()
{
   int ctr = 0;
   tLNode *cn, *dn, *ln;
   tGridListIter< tLNode > nI( gridPtr->getNodeList() ),
       tI( gridPtr->getNodeList() );
   for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
   {
      //make sure each node has a viable downstrm nbr:
      dn = cn->getDownstrmNbr();
      if( dn == 0 )
      {
         cerr << "NODE #" << cn->getID() << " has no downstrm nbr\n";
         goto error;
      }
      if( ln = tI.GetP( dn->getID() ) )
      {
         if( ln != dn )
         {
            cerr << "NODE #" << cn->getID()
                 << " downstrm nbr not id. to node in nodeList with same ID\n";
            goto error;
         }
      }
      else
      {
         cerr << "NODE #" << cn->getID()
              << " downstrm nbr is not in nodeList\n";
         goto error;
      }
      if( cn->Meanders() )
      {
         if( cn->getSlope() < 0.0 )
         {
            cerr << "NODE #" << cn->getID()
                 << " meanders and returns negative getSlope\n";
            goto error;
         }
      }
   }
   
   for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
   {
      //make sure each node has path to outlet:
      dn = cn->getDownstrmNbr();
      while( dn->getBoundaryFlag() == kNonBoundary )
      {
         dn = dn->getDownstrmNbr();
         ctr++;
         if( ctr > kLargeNumber )
         {
            cerr << "NODE #" << cn->getID()
                 << " has infinite loop in path downstream\n";
            goto error;
         }
      }
      if( dn->getBoundaryFlag() != kOpenBoundary )
      {
         cerr << "NODE #" << cn->getID()
              << " does not flow to outlet\n";
         goto error;
      } 
   }
      
   cout << "NETWORK PASSED\n";

   return;
   
  error:
   ReportFatalError( "Error in network consistency." );
}
#undef kLargeNumber


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
	tGridListIter<tEdge> i( gridPtr->getEdgeList() );
  double slp, length;

  //cout << "CalcSlopes()...";

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
  //cout << "CalcSlopes() finished" << endl;	
}


// NB: TODO: moved to tGrid; delete!
// TODO: change L-hand to R-hand orientation of Voronoi vertices and
// create faster Voronoi edge length computation scheme
/*Xvoid tStreamNet::CalcVAreas()
{
   //cout << "CalcVAreas()..." << endl << flush;
   double area;
   tLNode * curnode, *cn;
   tEdge *ce;
   tArray< double > xy;
   tGridListIter<tLNode> nodIter( gridPtr->getNodeList() ),
       nI( gridPtr->getNodeList() );
   tPtrListIter< tEdge > sI;
   
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
         //area = VoronoiArea( curnode );
      curnode->ComputeVoronoiArea();
      
   }
   //cout << "CalcVAreas() finished" << endl;
}
*/


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
   tGridListIter<tLNode> i( gridPtr->getNodeList() );
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
      flowedg = curnode->getEdg();
      assert( flowedg>0 );

      // As long as the current edge is a no-flow edge, advance to the next one
      // counter-clockwise
      ctr = 0;
      while( !flowedg->FlowAllowed() )
      {
         flowedg = flowedg->getCCWEdg();
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
      curnode->setFlowEdg( flowedg );
      assert( curnode->getFlowEdg() > 0 );
      curnode = i.NextP();
   }
}
#undef kMaxSpokes

/****************************************************************************\
**
**  TryDamBypass: See whether we can avoid flowing uphill by deleting nearby
**     points.
**
**      Parameters:     none
**      Called by:
**      Assumes:
**         node points to a valid flow edge which connects the lowest nbr
**
\****************************************************************************/
int tStreamNet::DamBypass( tLNode *snknod )
{
   cout << "DamBypass" << endl;
   int nv, nvopp, nvother, maxcnt = 10, cntr = 0;
   double cz, nz, nnz, dis0, dis1, slp;
   tTriangle *mytri, *opptri;
   int uphill = 1;
   tEdge *fedg, *ce;
   tLNode *nbr;
   tLNode *nxtnbr;
   tLNode *othernbr;
   tLNode *pointtodelete;
   tArray< double > cxy, nxy, nnxy;
   tPtrListIter< tEdge > sI( snknod->getSpokeListNC() );
   cxy = snknod->get2DCoords();
   
   cz = snknod->getZ();
   //while( uphill && cntr <= maxcnt )
   //{
      fedg = snknod->getFlowEdg();
      nbr = snknod->getDownstrmNbr();
      nxtnbr = nbr->getDownstrmNbr();
      if( nxtnbr == 0 ) return 0;
      if ( snknod == nxtnbr ) return 0;
      nnz = nxtnbr->getZ();
      if( nnz > cz ) return 0; //proceed no further if next nbr is not downhill
      nxy = nbr->get2DCoords();
      nnxy = nxtnbr->get2DCoords();
      if( PointsCCW( cxy, nxy, nnxy ) )
      {
         mytri = gridPtr->TriWithEdgePtr( fedg );
         nv = mytri->nVtx( snknod );
         othernbr = (tLNode *) mytri->pPtr( (nv+2)%3 );
      }
      else
      {
         mytri = gridPtr->TriWithEdgePtr( gridPtr->getEdgeComplement( fedg ) );
         nv = mytri->nVtx( snknod );
         othernbr = (tLNode *) mytri->pPtr( (nv+1)%3 );
      }
      opptri = mytri->tPtr(nv);
      dis0 = nbr->Dist( snknod, nxtnbr );
      dis1 = othernbr->Dist( snknod, nxtnbr );
      if( dis0 < dis1 ) pointtodelete = nbr;
      else pointtodelete = othernbr;
      if( pointtodelete->getBoundaryFlag() == kNonBoundary )
          gridPtr->DeleteNode( pointtodelete );
      else return 0;
      cntr++;
      for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
      {
         nz = ce->getDestinationPtrNC()->getZ();
         slp = ( cz - nz ) / ce->getLength();
         ce->setSlope( slp );
         gridPtr->getEdgeComplement( ce )->setSlope( -slp );
         //if( ce->getDestinationPtrNC()->getZ() < cz ) uphill = 0;
      }
      //}
      cout << "DamBypass: deleted " << cntr << " nodes" << endl << flush;
   return 1;
}


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
**       - each node points to a valid flow edge (returned by getFlowEdg())
**       - each edge has a valid counter-clockwise edge
**      Updated: 12/19/97 SL; 12/30/97 GT
**
\****************************************************************************/
#define kLargeNegative -1000
#define kMaxSpokes 100
void tStreamNet::FlowDirs()
{
   tGridListIter<tLNode> i( gridPtr->getNodeList() );  // gets nodes from the list
   double slp;        // steepest slope found so far
   tLNode *curnode;  // ptr to the current node
   tLNode *newnode;  // ptr to new downstream node
   tEdge * firstedg; // ptr to first edg
   tEdge * curedg;   // pointer to current edge
   tEdge * nbredg;   // steepest neighbouring edge so far
   long seed = 91324;
   double chngnum;
   int ctr;
   
//#if TRACKFNS
   //cout << "FlowDirs" << endl;
//#endif

   //int redo = 1;
   //while( redo )
   //{
   //   redo = 0;
      // Find the connected edge with the steepest slope
      curnode = i.FirstP();
      while( i.IsActive() )  // DO for each non-boundary (active) node
      {
         curnode->setFloodStatus( kNotFlooded );  // Init flood status flag
         firstedg =  curnode->getFlowEdg();
         if( firstedg <= 0 )
             curnode->TellAll();
         assert( firstedg > 0 );
         slp = firstedg->getSlope();
         nbredg = firstedg;
         curedg = firstedg->getCCWEdg();
         ctr = 0;
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
            curedg = curedg->getCCWEdg();
            ctr++;
            if( ctr>kMaxSpokes ) // Make sure to prevent endless loops
            {
               cerr << "Mesh error: node " << curnode->getID()
                    << " going round and round"
                    << endl;
               ReportFatalError( "Bailing out of FlowDirs()" );
            }
         }
     //add a wrinkle: if node is a meander node and presently flows
     //to another meander node and the new 'nbredg' does not lead to a
     //meander node, then choose a random number and
     //compare it to the probability that a meander node will change
     //flow direction to a non-meander node
           /*if( mndrDirChngProb != 1.0 )
         {
            newnode = (tLNode *) nbredg->getDestinationPtrNC();
            if( curnode->getDownstrmNbr()->Meanders() &&
                curnode->getDownstrmNbr()->getZ() < curnode->getZ() &&
                !(newnode->Meanders()) )
            {
               chngnum = ran3( &seed );
               if( chngnum <= mndrDirChngProb ) curnode->setFlowEdg( nbredg );
            }
            else curnode->setFlowEdg( nbredg );
         }
         else curnode->setFlowEdg( nbredg );*/
      
         curnode->setFlowEdg( nbredg );

           /*if( slp < 0 )
         {
            if( DamBypass( curnode ) )
            {
               InitFlowDirs();
               redo = 1;
               break;
            }
         }*/
         curnode->setFloodStatus( ( slp>0 ) ? kNotFlooded : kSink );  // (NB: opt branch pred?)
         //cout << "Node " << curnode->getID() << " flows to "
         //     << curnode->getDownstrmNbr()->getID() << endl;
         curnode = i.NextP();

      }
      //}
   //cout << "FlowDirs() finished" << endl << flush;
  
}
#undef kLargeNegative
#undef kMaxSpokes





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
   //cout << "DrainAreaVoronoi()..." << endl << flush;
//#endif
   tLNode * curnode;
   tGridListIter<tLNode> nodIter( gridPtr->getNodeList() );
   
   // Reset drainage areas to zero
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
       curnode->setDrArea( 0 );
   
   // send voronoi area for each node to the node at the other end of the 
   // flowedge and downstream 
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      RouteFlowArea( curnode, curnode->getVArea() );
   }
   if( inlet.innode != 0 )
   {
      inlet.FindNewInlet();
      RouteFlowArea( inlet.innode, inlet.inDrArea );
   }

   //cout << "DrainAreaVoronoi() finished" << endl << flush;
}




/*****************************************************************************\
**
**  tStreamNet::RouteFlowArea
**
**  Starting with the current node, this routine increments 
**  the drainage area of the node and each node downstream by _addedArea_.
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
   while( !(curnode->getBoundaryFlag()) && (curnode->getFloodStatus()!=kSink) )
   {
      curnode->AddDrArea( addedArea );
      curnode = curnode->getDownstrmNbr();
//#if DEBUG
      niterations++;
      if( niterations>9990 )
      {
         curnode->TellAll();
         cout  << flush;
      }
      assert( niterations < 10000 );
//#endif
   }
   //cout << "RouteFlowArea() finished" << endl << flush;
}


/*****************************************************************************\
**
**  tStreamNet::RouteRunoff
**
**  This is a variant of RouteFlowArea that adds addedArea and addedRunoff
**  to the drainage area and discharge, respectively, of each downstream node.
**
**  Parameters: curnode -- node to start with
**              addedArea -- drainage area to add to downstream nodes
**              addedRunoff -- runoff rate to add to downstream nodes
**  Created: 6/98 GT
**
\*****************************************************************************/
void tStreamNet::RouteRunoff( tLNode *curnode, double addedArea,
                              double addedRunoff )
{
   //cout << "RouteFlowArea()..." << endl << flush;
//#if DEBUG
   int niterations=0;  // Safety feature: prevents endless loops
//#endif

   // As long as the current node is neither a boundary nor a sink, add
   // _addedArea_ to its total drainage area and _addedRunoff_ to its total
   // discharge and then advance to the next node downstream.
   while( !(curnode->getBoundaryFlag()) && (curnode->getFloodStatus()!=kSink) )
   {
      curnode->AddDrArea( addedArea );
      curnode->AddDischarge( addedRunoff );
      curnode = curnode->getDownstrmNbr();
//#if DEBUG
      niterations++;
      if( niterations>9990 )
      {
         curnode->TellAll();
         cout  << flush;
      }
      assert( niterations < 10000 );
//#endif
   }
   //cout << "RouteFlowArea() finished" << endl << flush;
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
   //cout << "MakeFlow()..." << endl << flush;
   if( filllakes ) FillLakes();
   DrainAreaVoronoi();
   switch( flowgen )
   {
      case kSaturatedFlow1:   // Topmodel-like runoff model with steady-state
          FlowSaturated1();   //   surface flow and possible return flow.
          break;
      case kSaturatedFlow2:   // Topmodel-type runoff model without return
          FlowSaturated2();   //   flow (infiltrated water does not contribute
          break;              //   to storm flow)
      case kConstSoilStore:   // "Bucket" model: spatially uniform soil storage
          FlowBucket();       //   capacity; any excess contributes to runoff.
          break;
      default:
          FlowUniform();      // Spatially uniform infiltration-excess runoff
   }
   
   //cout << "MakeFlow() finished" << endl;
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
   //cout << "FlowUniform..." << endl << flush;
   tGridListIter< tLNode > nodIter( gridPtr->getNodeList() );
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
   //cout << "FlowUniform finished" << endl;
}


/*****************************************************************************\
**
**  tStreamNet::FlowSaturated1
**
**  Computes surface runoff using Topmodel concept. Steady-state 
**  subsurface flow capacity is transmissivity times slope. Total
**  steady-state runoff is rainfall rate times drainage area. Surface
**  runoff is total minus subsurface. The assumption of steady-state
**  subsurface flow implies that subsurface storm flow can contribute
**  significantly to erosion (in other words, both saturation-excess
**  and return flow are modeled). To model cases in which return flow
**  is insigificant, use FlowSaturated2.
**
**  Transmissivity has units of L^2/T; it is hydraulic conductivity
**  times depth. To get a volume discharge, we must multiply the
**  transmissivity by the Voronoi edge length for the flow edge.
**
**  Parameters:  none
**  Calls:  various utilities of tGridListIter and tLNode
**  Called by:  main
**  Modifies: node discharge
**  Modifications:
**    - updated 12/19/97 SL
**    - 5/2/98 SL
**
\*****************************************************************************/
void tStreamNet::FlowSaturated1()
{
   //cout << "FlowSaturated1...";
   tGridListIter< tLNode > nodIter( gridPtr->getNodeList() );
   tLNode *curnode;
   tEdge *fedg;
   double discharge;
   
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      fedg = curnode->getFlowEdg();
      discharge = curnode->getDrArea() * rainrate -
          fedg->getSlope() * fedg->getVEdgLen() * trans;
      discharge *= ( discharge > 0 );
      curnode->setDischarge( discharge );
   }
   //cout << "finished" << endl;
}


/*****************************************************************************\
**
**  tStreamNet::FlowSaturated2
**
**  Computes surface runoff using a more precise formulation of the
**  Topmodel concept, based on Ijjasz-Vasquez et al (1992). Total
**  runoff from a node is the sum of infiltration-excess (if any)
**  and saturation-excess (if any). Saturation-excess is determined
**  by computing a saturation deficit that depends on slope, area,
**  and a parameter (trans), and subtracting the deficit from the storm
**  depth to obtain the saturation-excess runoff depth, which is
**  converted to a rate by multiplying by the storm duration.
**
**  Parameters:  none
**  Called by:  main
**  Modifies:  node discharges
**
\*****************************************************************************/
void tStreamNet::FlowSaturated2()
{
   tGridListIter< tLNode > nodIter( gridPtr->getNodeList() );
   tLNode *curnode;
   tEdge *fedg;
   double discharge;
   double infiltExRunoff,  // Infiltration-excess runoff rate (L/T)
       sd,                 // Storm depth minus any infilt-excess runoff (L)
       asRatio,            // Slope-area-per-width ratio (ie, A/Sb) (L),
                           //   where b is Voronoi edge length
       satDeficit,         // Saturation deficit (L)
       rsat,               // Saturation-excess runoff depth (L)
       runoff,             // Total runoff from node (L/T)
       stormDur = stormPtr->getStormDuration();  // Storm duration
   int nsat=0,nsr=0,nhort=0,nflat=0; // 4dbg
  
#if TRACKFNS
  cout << "FlowSaturated1" << endl << flush;
#endif

  // Reset drainage areas and discharges to zero
  for ( curnode=nodIter.FirstP(); nodIter.IsActive(); curnode=nodIter.NextP() )
  {
     //if( curnode->drarea>1.15e7 )
     //    cout << "Q(" << curnode->id << ") " << curnode->q << endl;
     curnode->setDischarge( 0 );
  }
  
  for ( curnode=nodIter.FirstP(); nodIter.IsActive(); curnode=nodIter.NextP() )
  {
     //cout<<"Node " <<curnode->getID();
     assert( curnode->getFlowEdg() != 0 );
     satDeficit = 0;
     infiltExRunoff = rainrate - infilt;
     if( infiltExRunoff<0 ) infiltExRunoff = 0;
     else nhort++;
     sd = (rainrate - infiltExRunoff)*stormDur;
     if( curnode->getSlope()>0 )
     {
         asRatio = curnode->getDrArea() /
             ( curnode->getSlope() * curnode->getFlowEdg()->getVEdgLen() );
         assert( asRatio>0 );
         if( asRatio < trans  ) // if true, then NOT saturated => deficit
         {
            assert( trans>0 && asRatio>0 );
            satDeficit = log( trans / asRatio );
         }
         else nsat++;
         //cout << " satdef " << satDeficit;
     }
     else {
        //cout<<" zero slp or ve";
        nflat++;
        nsat++;
     }
     rsat = sd - satDeficit;
     if( rsat<0 ) rsat = 0;
     else nsr++;
     runoff = infiltExRunoff + rsat/stormDur;
     //cout<<" sat excess " << rsat << " total " << runoff << endl;
     RouteRunoff( curnode, curnode->getVArea(), runoff*curnode->getVArea() );
  }

  //cout << nhort << " generate Horton runoff, " << nsat
  //     << " pre-saturated, (" << nflat << "flat) "
  //     << nsr << " saturate during storm.\n";
  
}


/*****************************************************************************\
**
**  tStreamNet::FlowBucket
**
**  Computes runoff rate and discharge assuming a spatially and temporally
**  uniform "bucket" soil storage capacity. Runoff rate is the sum of
**  infiltration excess runoff ( = rainrate - infilt ), if any, and
**  saturation excess runoff, if any. Saturation excess runoff is
**  calculated as ( total water depth of infiltrated storm water
**  minus soil storage capacity ) divided by storm duration, where the
**  infiltration depth is just the rainfall depth minus any infiltration-
**  excess runoff. Discharge is runoff rate times contributing area.
**    We save a multiplication operation by noting that:
**      satEx = ( infiltDepth - soilStore ) / stormDuration
**            = (rainrate - infiltEx) - (soilStore/stormDuration)
**
**  Parameters: none
**  Called by:  MakeFlow
**  Assumes:  stormPtr is valid -- ie, to use this function, StreamNet
**                                  must be init'd with a storm object,
**                                  or UpdateNet must be called w/ a
**                                  storm object.
**  Created:  5/98 GT
**  Modifications:
**
\*****************************************************************************/
void tStreamNet::FlowBucket()
{
   assert( stormPtr!=0 );
   //cout << "FlowBucket..." << endl << flush;
   tGridListIter< tLNode > nodIter( gridPtr->getNodeList() );
   tLNode *curnode;
   double infiltEx=0.0,  /* infiltration excess runoff (L/T) */  
       satEx=0.0,        /* saturation excess runoff (L/T) */
       infiltRate,       /* soil store / storm duration */
       runoff,           /* total runoff */
       discharge;
   
   // Compute total runoff as the sum of infiltration excess (if any)
   // and saturation excess (if any)
   rainrate = stormPtr->getRainrate();
   if( rainrate > infilt )
       infiltEx = rainrate - infilt;
   if( (rainrate-infiltEx) >
       (infiltRate=(soilStore/(stormPtr->getStormDuration()) ) ) )
       satEx = (rainrate - infiltEx) - infiltRate;
   runoff = infiltEx + satEx;
   cout << "  R " << runoff << " = IEx " << infiltEx << " + SEx " << satEx << endl;
   
   // Compute and assign discharge for each node
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      discharge = curnode->getDrArea() * runoff;
      curnode->setDischarge( discharge );
   }
   //cout << "FlowBucket finished" << endl;
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
   //cout << "FillLakes()..." << endl << flush;
//#endif
   tLNode *cn,    // Node on list: if a sink, then process
       *thenode,           // Node on lake perimeter
       *lowestNode,        // Lowest node on perimeter found so far
       *cln, // current lake node
       *node; //placeholder
   tGridListIter< tLNode > nodIter( gridPtr->getNodeList() );
   tPtrList< tLNode > lakeList;
   tPtrListIter< tLNode > lakeIter( lakeList );
   //tPtrListIter< tEdge > spokIter;
   tEdge *ce;           // Pointer to an edge
   double lowestElev;    // Lowest elevation found so far on lake perimeter
   int done;          // Flag indicating whether outlet has been found
   
   // Check each active node to see whether it is a sink
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      if( cn->getFloodStatus() == kSink )
      {
         // Create a new lake-list, initially containing just the sink node.
         lakeList.insertAtBack( cn );
         cn->setFloodStatus( kCurrentLake );
         
         // Iteratively search for an outlet along the perimeter of the lake
         done = FALSE;
         do
         {
            lowestNode = lakeIter.FirstP();
            lowestElev = kVeryHigh; // Initialize lowest elev to very high val.
            // Check the neighbors of every node on the lake-list
            for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
                 cln = lakeIter.NextP() )
            {
               //Xcout << "LAKE LIST:\n";
               //Xcln->TellAll();
                      //XspokIter.Reset( cln->getSpokeListNC() );
                   //Xfor( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                   //Xce = spokIter.NextP() )
               // Check all the neighbors of the node
               ce = cln->getEdg();
               do
               {
                  thenode = (tLNode *) ce->getDestinationPtrNC();
                  // Is it a potential outlet (ie, not flooded and not
                  // a boundary)?
                  if( thenode->getFloodStatus() == kNotFlooded
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
                  else if( thenode->getFloodStatus() == kFlooded ||
                           thenode->getFloodStatus() == kSink )
                  {
                     lakeList.insertAtBack( thenode );
                     thenode->setFloodStatus( kCurrentLake );
                  }
               } while( ( ce=ce->getCCWEdg() ) != cln->getEdg() );// END spokes
               
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
                  lowestNode->setFloodStatus( kCurrentLake );
               }
            }
            if( lakeList.getSize() > gridPtr->getNodeList()->getActiveSize() )
            {
               cout << "LAKE LIST SIZE: " << lakeList.getSize() << endl;
            }
            
            assert( lakeList.getSize() <= gridPtr->getNodeList()->getActiveSize() );
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
         lowestNode->setFloodStatus( kOutletFlag );
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
               if( cln->getFloodStatus() != kOutletFlag )
               {
                  done = FALSE;
                  /*XspokIter.Reset( cln->getSpokeListNC() );
                  for( ce = spokIter.FirstP();
                       cln->getFloodStatus() != kOutletFlag &&
                           !( spokIter.AtEnd() );
                       ce = spokIter.NextP() )*/
                  // Check each neighbor
                  ce = cln->getEdg();
                  do
                  {
                     node = (tLNode *) ce->getDestinationPtrNC();
                     if( node->getFloodStatus() == kOutletFlag )
                     {     // found one!  
                        cln->setFloodStatus( kOutletPreFlag );
                        cln->setFlowEdg( ce );
                        //cout << "Node " << cln->getID() << " flows to "
                        //     << cln->getDownstrmNbr()->getID() << endl;
                        
                     }
                  } while( cln->getFloodStatus() != kOutletFlag
                           && ( ce=ce->getCCWEdg() ) != cln->getEdg() );
               } // END if node not flagged as outlet
            } // END for each lake node

            // Now flag all the "preflagged" lake nodes as outlets
            for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
                 cln = lakeIter.NextP() )
                if( cln->getFloodStatus() == kOutletPreFlag )
                    cln->setFloodStatus( kOutletFlag );
            
         } while( !done );
         lowestNode->setFloodStatus( kNotFlooded );
         
         // Finally, flag all of the 
         // nodes in it as "kFlooded" and clear the list so we can move on to
         // the next sink. (Fixed mem leak here 8/5/97 GT).
         for( cln = lakeIter.FirstP(); !( lakeIter.AtEnd() );
              cln = lakeIter.NextP() )
             cln->setFloodStatus( kFlooded );
         lakeList.Flush();
      } /* END if Sink */
   } /* END Active Nodes */
   //cout << "FillLakes() finished" << endl << flush;
   
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
**  (4) If the outlet is itself part of a different lake, the water surface
**      elevation of that lake must be lower than the current node.
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
   tLNode *dn,        // Potential outlet
       *an;           // Node ptr used to find outlet of a previously
                      // identified lake
   
   // Check all the neighbors
   ce = node->getEdg();
   do
   {
      // If it passes this test, it's a valid outlet
      dn = (tLNode *) ce->getDestinationPtrNC();
      assert( dn>0 );
/*X      if( ce->getSlope() > maxslp &&
          dn->getFloodStatus() != kCurrentLake &&
          ce->FlowAllowed() &&
          ( dn->getBoundaryFlag()==kOpenBoundary ||
            dn->getDownstrmNbr()->getZ() < node->getZ() ) )*/
      if( ce->getSlope() > maxslp &&
          dn->getFloodStatus() != kCurrentLake &&
          ce->FlowAllowed() )
      {
         // Handle a very special and rare case: if the "target" node dn is
         // part of a previous lake, it's still a valid exit as long as its
         // water surface elevation is lower than the current lake (whose
         // wse, assuming an outlet is found, would be equal to _node_'s
         // elevation). It can sometimes happen that the target lake's wse is
         // exactly equal in elevation to _node_, in which case
         // the point is not considered an outlet---if it were, infinite loops
         // could result. (This fix added 4/98)
         if( dn->getFloodStatus()==kFlooded )
         {
            // Iterate "downstream" through the "old" lake until reaching the
            // outlet, then test its elevation. If the elevation is exactly
            // equal to _node_, skip the rest and go on to the next iteration.
            an = dn;
            while( an->getFloodStatus()!=kNotFlooded )
                an = an->getDownstrmNbr();
            if( an->getZ()==node->getZ() ) continue;
         }

         // Assign the new max slope and set the flow edge accordingly
         maxslp = ce->getSlope();
         node->setFlowEdg( ce );
         //cout << "Node " << node->getID() << " flows to "
         //     << node->getDownstrmNbr()->getID() << endl;
      }
   } while( ( ce=ce->getCCWEdg() ) != node->getEdg() );
   
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
   int done=0;
   tLNode * cn;
   tGridList<tLNode> *nodeList = gridPtr->getNodeList();
   int nUnsortedNodes = nodeList->getActiveSize();  // Number not yet sorted
   tGridListIter<tLNode> listIter( nodeList );
   
   //test
   /*Xcout << "BEFORE: " << endl;
   for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() ) 
       cout << cn->getID() << endl;*/
   
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
         assert( cn!=0 );
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
            nodeList->moveToActiveBack( nodeToMove );
            nThisPass++;
         }
         else
         {
            cn = listIter.NextP();
            done = FALSE;
         }
      }
      
      nUnsortedNodes -= nThisPass;

      /*cout << "NO. UNSORTED: " << nUnsortedNodes << endl;
      for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() ) 
          cout << cn->getID() << " " << cn->getQ() << " " << cn->getQs()
               << endl;*/

    } while( !done );

   /*Xcout << "AFTER: " << endl;
  cn = listIter.FirstP();
  cout << "First node:\n";
  cn->TellAll();
  for( cn=listIter.FirstP(); listIter.IsActive(); cn=listIter.NextP() ) 
      cout << cn->getID() << " " << cn->getQ() << endl;
  cout << "Leaving Sort\n" << flush;*/
  
 
}

/*****************************************************************************\
**
**       FindHydrGeom: goes through reach nodes and calculates/assigns
**                     values for hydraulic geometry.
**                
**
**		  Parameters: kwds, ewds, ewstn, knds, ends, enstn, optrainvar
**      Data members updated: tLNode->chan.hydrwidth, hydrnrough, hydrdepth
**      Called by: FindReaches (needs to know how long to make reach "tails"
**      Calls: no "major" functions
**
**      Created: SL 1/98         
**
**
\*****************************************************************************/
void tStreamNet::FindHydrGeom()
{
   int i, j, num;
   double hradius, kwdspow, kndspow, widpow, npow, radfactor, qpsec;
   double width, depth, rough, slope;
   tLNode *cn;

   kwdspow = pow(kwds, ewstn / ewds);
   kndspow = pow(knds, enstn / ends);
   widpow = 1.0 - ewstn / ewds;
   npow = 1.0 - enstn / ends;
   //timeadjust = 86400 * days;  /* 86400 = seconds in a day */
   tGridListIter< tLNode > nIter( gridPtr->getNodeList() );
   for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
   {
      //removed an if cn->Meanders(), so stuff calculated everywhere
      //if rainfall varies, find hydraulic width "at-a-station"
      //based on the channel width "downstream":
      if( optrainvar)
      {
         qpsec = cn->getQ();
         width = pow(cn->getChanWidth(), widpow) * kwdspow * pow(qpsec, ewstn);
         cn->setHydrWidth( width );
         rough = pow(cn->getChanRough(), npow) * kndspow * pow(qpsec, enstn);
         cn->setHydrRough( rough );
         slope = cn->getChanSlope();
         assert( slope > 0 );
         radfactor = qpsec * rough / width / sqrt(slope);
         hradius = pow(radfactor, 0.6);
         depth = width / ( width / hradius - 2.0 );
         cn->setHydrSlope( slope );
         cn->setHydrDepth( depth );
      }
      //if rainfall does not vary, set hydraulic geom. = channel geom.
      else
      {
         width = cn->getChanWidth();
         rough = cn->getChanRough();
         depth = cn->getChanDepth();
         slope = cn->getChanSlope();
         cn->setHydrWidth( width );
         cn->setHydrRough( rough );
         cn->setHydrSlope( slope );
         cn->setHydrDepth( depth );
      }
   }
   //cout << "done FindHydrGeom" << endl << flush;
}

/*****************************************************************************\
**
**       FindChanGeom: goes through reach nodes and calculates/assigns
**                     values for channel geometry.
**                
**
**		  Parameters: kwds, ewds, ewstn, knds, ends, enstn
**      Data members updated: tLNode->chan.hydrwidth, hydrnrough, hydrdepth
**      Called by: FindReaches (needs to know how long to make reach "tails"
**      Calls: no "major" functions
**
**      Created: SL 1/98         
**
**
\*****************************************************************************/
#define kSmallNum 0.0000000001
void tStreamNet::FindChanGeom()
{
   int i, j, num;
   double qbf, hradius, qbffactor=0, radfactor, width, depth, rough, slope;
   double lambda;
   double rlen, cz, nz, critS;
   tLNode *cn, *dsn;
   tGridListIter< tLNode > nIter( gridPtr->getNodeList() );
   tPtrList< tLNode > *plPtr;
   //timeadjust = 86400 * days;  /* 86400 = seconds in a day */
   tStorm *sPtr = getStormPtrNC();
   double isdmn = sPtr->getMeanInterstormDur();
   double pmn = sPtr->getMeanPrecip();
   if (isdmn > 0 )  qbffactor = pmn * log(1.5 / isdmn);
   for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
   {
      //took out an if cn->Meanders() so stuff will be calculated at all nodes
// qbffactor is now in m^3/s
      qbf = cn->getDrArea() * qbffactor;
      if( !qbf ) qbf = cn->getQ();  // q is now in m^3/s
      width = kwds * pow(qbf, ewds);
      rough = knds * pow(qbf, ends);
      lambda = klambda * pow(qbf, elambda);
      cn->setChanWidth( width );
      cn->setChanRough( rough );
      cn->setBankRough( lambda );
      slope = cn->getSlope();
      //make sure slope will produce a positive depth:
      critS = qbf * qbf * rough * rough * 8.0 * pow( 2.0, 0.333 )/
          ( width * width * width * width * width * pow( width, 0.333 ) );
      if( slope > critS ) //should also catch negative slope flag
      {
         //cout << "in FindChanGeom, slope = " << slope << endl << flush;
         cn->setChanSlope( slope );
         radfactor = qbf * rough / width / sqrt(slope);
         hradius = pow(radfactor, 0.6); 
         depth = width / (width / hradius - 2.0);
         cn->setChanDepth( depth );
      }
      else cn->setMeanderStatus( kNonMeanderNode );
      if( slope < 0.0 )
      {
         cout << "negative slope,"
              << " probably from infinite loop in tLNode::GetSlope()" << endl;
         ReportFatalError("negative slope in tStreamMeander::FindChanGeom");
      }
   }
   //cout << "done FindChanGeom" << endl;
}
#undef kSmallNum

/**************************************************************************\
**
**  tStreamNet.cpp
**
**  Functions for class tInlet.
**
\**************************************************************************/
tInlet::tInlet()
{
   innode = 0;
   inDrArea = 0;
   gridPtr = 0;
}

/**************************************************************************\
**
**  tInlet Constructor #2
**
**  ...
**
**  Bug: if a new inlet node is added, using AddNodeAt won't assign
**  correct variables for layers, regolith, etc. Fix: new node is assigned
**  properties of the nearest neighbor (not elevation; that's done by
**  interpolation). GT 7/98
**
\**************************************************************************/
#define LARGE_DISTANCE 1e9
tInlet::tInlet( tGrid< tLNode > *gPtr, tInputFile &infile )
{
   int i, inletbc = infile.ReadItem( inletbc, "OPTINLET" ),
       numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   int add = 1;
   char end, name[20];
   double xin, yin,   // Coords of inlet node
       mindist,       // Minimum distance above which new node will be added
       dist,          // Distance btwn inlet and a nearby node
       x, y,          // Location of a nearby node
       zin = 0,       // Elevation of inlet
       suminvdist = 0,// Sum of 1/dist for all nearby non-boundary nodes
       minDistFound;  // Smallest distance to nearby node found so far
   double help;
   tArray< double > xyz(3);
   tTriangle *intri, *ntri;
   tLNode *cn,
       *closestNode;  // -> to closest nearby node found so far
   tPtrList< tLNode > nPL;            // List of nearby non-boundary nodes
   tPtrListIter< tLNode > itr( nPL ); // Iterator for the above list
   
   gridPtr = gPtr;
   assert( gridPtr != 0 );
   if( inletbc )
   {
      // Read drainage area and sediment load at inlet. If more than one
      // grain size is simulated, read in a sediment load for each size
      // individually
      inDrArea = infile.ReadItem( inDrArea, "INDRAREA" );
      if(numg <= 1)
          inSedLoad = infile.ReadItem( inSedLoad, "INSEDLOAD" );
      else{
         inSedLoadm.setSize(numg+1);
         inSedLoadm[0]=0.0;
         i=1;
         end='1';
         while( i<=numg ){
            strcpy( name, "INSEDLOAD");
            strcat( name, &end ); 
            help = infile.ReadItem( help, name);
            inSedLoadm[i] = help;
            inSedLoadm[0] += help;
            i++;
            end++;
         }
      }

      // Read in the location of the inlet node. If the specified coordinates
      // are "close" to an existing non-boundary node, assign that node as
      // the inlet; otherwise, create a new node. The elevation for the new
      // node is found by interpolation.
      xin = infile.ReadItem( xin, "INLET_X" );
      yin = infile.ReadItem( yin, "INLET_Y" );
      intri = gridPtr->LocateTriangle( xin, yin );
      assert( intri > 0 );
      for( i=0; i<3; i++ )
      {
         cn = (tLNode *) intri->pPtr(i);
         if( cn->getBoundaryFlag() == kNonBoundary ) nPL.insertAtBack( cn );
         ntri = intri->tPtr(i);
         if( ntri != 0 )
         {
            cn = (tLNode *) ntri->tPtr( ntri->nVOp( intri ) );
            if( cn->getBoundaryFlag() == kNonBoundary ) nPL.insertAtBack( cn );
         }
      }
      minDistFound = LARGE_DISTANCE;
      mindist = 0.000001;
      for( cn = itr.FirstP(); !(itr.AtEnd()); cn = itr.NextP() )
      {
         x = cn->getX();
         y = cn->getY();
         dist = sqrt( (xin - x) * (xin - x) + (yin - y) * (yin - y) );
         if( dist < minDistFound )
         {
            minDistFound = dist;
            closestNode = cn;
         }
         if( dist > mindist )
         {
            zin += cn->getZ() / dist;
            suminvdist += 1 / dist;
         }
         else
         {
            innode = cn;
            add = 0;
         }
         
         /*
         if( dist < mindist )
         {
            dist = mindist;
            innode = cn;
         }*/
      }
      if( add ) // fix here:
      {
         // NOTE: this leads to a new bug, as yet unfixed. For unknown
         // reasons, the dgrade array of the top layer of closestNode becomes
         // corrupted, but apparently not until after this routine... so it's
         // an indirect effect. I suspect the problem lies in the copy
         // constructor (or lack of one) in one of component objects that
         // are being copied: tLNode -> tRegolith -> layerlist -> dgrade arr?
         cout << "ADDING INLET: Closest node is:" << endl;
         closestNode->TellAll();
         cout << "DGRADE BEFORE: " << closestNode->getLayerDgrade(0,0);
         tLNode newnode( *closestNode );
         cout << " DGRADE AFTER: " << closestNode->getLayerDgrade(0,0) << endl << flush;
         newnode.setZ( zin / suminvdist );
         newnode.setX( xin );
         newnode.setY( yin );
         //zin = zin / suminvdist;
         //xyz[0] = xin;
         //xyz[1] = yin;
         //xyz[2] = zin;
         //innode = gridPtr->AddNodeAt( xyz );
         innode = gridPtr->AddNode( newnode );
         cout << " DGRADE AFTER AFTER: " << closestNode->getLayerDgrade(0,0) << endl << flush;
         cout << "INLET NODE IS:\n";
         innode->TellAll();
      }
   }
   else
   {
      inDrArea = 0;
      innode = 0;
   }
}

tInlet::~tInlet()
{
   innode = 0;
   gridPtr = 0;
}


/**************************************************************************\
**
**  tInlet::FindNewInlet: search for points 'up-valley' of present inlet;
**     of those points and the present inlet, set new inlet to one with the
**     lowest elevation.
**  5/18/98 SL: Try something less arbitrary. Find active nodes
**
\**************************************************************************/

void tInlet::FindNewInlet()
{
   double xin, yin, zin, x, y, z, zmin, dmn, dmnn, dmin;
   tLNode *cn, *newinnode, *mn;
   tNode *bn0, *bn1, *mnn;
   tEdge *ce, *me;
   tGridListIter< tLNode > nI( gridPtr->getNodeList() );
   tPtrListIter< tEdge > sI, msI;
   int n;
     //tPtrList< tLNode > bList;
     //tPtrListIter< tLNode > bI( bList );
   tArray< double > xyz = innode->get3DCoords();
   yin = xyz[1];
   zmin = xyz[2];
   newinnode = innode;
     //for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
     //go through boundary nodes:
   cn = nI.LastActiveP();
   for( cn = nI.NextP(); !(nI.AtEnd()); cn = nI.NextP() )
   {
        //select for 'northern' bndy nodes:
      if( cn->getY() > yin ) //(cn was originally any active node)
      {
           //go through bndy node's nbrs...
         sI.Reset( cn->getSpokeListNC() );
         for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
         {
            mn = (tLNode *) ce->getDestinationPtrNC();
              //easier to check node's elevation at this point to make
              //sure it's worth going further with all this logic;
              //so, find an active neighbor with elevation lower than
              //present inlet:
            if( mn->getBoundaryFlag() == kNonBoundary && mn->getZ() < zmin )
            {
               msI.Reset( mn->getSpokeListNC() );
               n = 0;
                 //go through the active node's nbrs to find and count
                 //'northern' bndy nodes:
               for(  me = msI.FirstP(); !(msI.AtEnd()); me = msI.NextP() )
               {
                  mnn = me->getDestinationPtrNC();
                  if( mnn->getBoundaryFlag() != kNonBoundary && mnn->getY() > yin )
                  {
                     if( n == 0 ) bn0 = mnn;
                     else if( n == 1 ) bn1 = mnn;
                     n++;
                  }
               }
                 //if active node has >1 northern bndy nbr...
               if( n > 1 )
               {
                    //find node's distance to boundary:
                  dmn = DistanceToLine( mn->getX(), mn->getY(), bn0, bn1 );
                  dmin = dmn;
                    //find it's active nbrs' distances:
                  for(  me = msI.FirstP(); !(msI.AtEnd()); me = msI.NextP() )
                  {
                     mnn = (tLNode *) me->getDestinationPtrNC();
                     if( mnn->getBoundaryFlag() == kNonBoundary )
                     {
                        dmnn = DistanceToLine( mnn->getX(), mnn->getY(),
                                                    bn0, bn1 );
                        if( dmnn < dmin ) dmin = dmnn;
                     }
                  }
                    //if none of active nbrs' distances smaller,
                    //node is the new inlet if no better candidates are found:
                  if( dmin == dmn )
                  {
                     zmin = mn->getZ();
                     newinnode = mn;
                  }
               }
            }
         }           
           /*if( cn->getZ() < zmin )
         {
            zmin = cn->getZ();
            newinnode = cn;
         }*/
      }
   }
     //finally reset inlet node; if no new inlet was found, we're just
     //setting it equal to itself:
   innode = newinnode;
}

double tInlet::getInSedLoad() const {return inSedLoad;}

double tInlet::getInSedLoad( int i )  
{
   if(i>inSedLoadm.getSize()-1)
        ReportFatalError( "Trying to set size in sediment load that doesn't exist");
   return inSedLoadm[i];
}


tArray< double >
tInlet::getInSedLoadm( ) const
{
   return inSedLoadm;
}


void tInlet::setInSedLoad( double val ) {inSedLoad = ( val > 0.0 ) ? val : 0.0;}
void tInlet::setInSedLoad( int i, double val )
{
   if(i>inSedLoadm.getSize()-1)
        ReportFatalError( "Trying to set size in sediment load that doesn't exist");
   inSedLoadm[i]=val;
}
double tInlet::getInDrArea() const {return inDrArea;}
void tInlet::setInDrArea( double val ) {inDrArea = ( val > 0.0 ) ? val : 0.0;}
tLNode *tInlet::getInNodePtr() {return innode;}
void tInlet::setInNodePtr( tLNode *ptr ) {innode = ( ptr > 0 ) ? ptr : 0;}

   
