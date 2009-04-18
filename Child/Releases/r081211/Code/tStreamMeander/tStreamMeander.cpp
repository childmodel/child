/**************************************************************************/
/**
**  @file tStreamMeander.cpp
**  @brief Functions for class tStreamMeander.
**
**  $Id: tStreamMeander.cpp,v 1.108 2005/03/15 17:17:30 childcvs Exp $
*/
/**************************************************************************/

#include "tStreamMeander.h"
#define kBugTime 5000000

#include "meander.h"

#if 0
/*****************************************************************************\
**
**
**      LineRemainder: given x,y coords, finds d = ax + by + c for the line
**              formed by points  p0->(x, y) and p1->(x, y)
**      Global function.
**      Data members updated:
**      Called by:
**      Calls:
**
**
\*****************************************************************************/
static
double LineRemainder( double x, double y, tNode * p0,tNode * p1 )
{
  double a,b,c, x0, y0, x1, y1;

  x0 = p0->getX();
  y0 = p0->getY();
  x1 = p1->getX();
  y1 = p0->getY();
  a = y1 - y0;
  b = x0 - x1;
  c = -( a * x0 + b * y0 );
  return (a * x + b * y + c);
}
#endif

/**************************************************************************\
**
**  Constructors
**
**  1) takes mesh and inputfile references and constructs a tStreamNet
**     object if ptr to tStreamNet is zero
**
**  2) takes mesh ptr and inputfile reference; assumes tStreamNet object
**     has already done its thing. Takes streamnet and inputfile references.
**
**     Modifications:
**       - fixed bug: NB was erroneously read in place of MB (gt 3/99)
**       - removed use of latadjust, which expressed the relationship
**         between KB and bank erodibility. Instead, the bank erodibility
**         is now read into rockerod as a completely separate parameter
**         "BANKERO" (GT 6/99)
**
\**************************************************************************/
tStreamMeander::tStreamMeander():
  meshPtr(0), netPtr(0), infilePtr(0),
  reachList(), rlIter(reachList),
#define FIXCRITFLOWBUG 1
#if FIXCRITFLOWBUG
  critarea(0),
#else
  critflow(0), 
#endif
  optdiamvar(false), optrainvar(false),
  meddiam(0), kwds(0), ewds(0), ewstn(0),
  knds(0), ends(0), enstn(0),
  kdds(0), edds(0), edstn(0),
  klambda(0), elambda(0), dscrtwids(0),
  allowfrac(0), leavefrac(0), /*vegerod(0), */
  rockerod(0), /*latadjust(0), */ Pdz(0),
  rand(0)
{
}

tStreamMeander::tStreamMeander( tStreamNet &netRef, tMesh< tLNode > &mRef,
                                tInputFile &infile, tRand *rand_ ) :
  reachList(), rlIter(reachList),
  rand(rand_) // seed already set in tMesh
{
  // Set pointers to various objects
  netPtr = &netRef;
  assert( netPtr != 0 );
  meshPtr = &mRef;
  assert( meshPtr != 0 );
  infilePtr = &infile;
  assert( infilePtr != 0 );
  
  // Read parameters from input file
#if FIXCRITFLOWBUG
  critarea = infilePtr->ReadItem( critarea, "CRITICAL_AREA" );
#else
  critflow = infilePtr->ReadItem( critflow, "CRITICAL_FLOW" );
#endif
  {
    int tmp_;
    tmp_ = infilePtr->ReadItem( tmp_, "OPT_VAR_SIZE" );
    optdiamvar = tmp_ != 0;
  }
  if( !optdiamvar )
    {
      meddiam = infilePtr->ReadItem( meddiam, "MEDIAN_DIAMETER" );
      assert( meddiam > 0 );
    }
  {
    int tmp_;
    tmp_ = infilePtr->ReadItem( tmp_, "OPTVAR" );
    optrainvar = tmp_ != 0;
  }
  kwds = infilePtr->ReadItem( kwds, "HYDR_WID_COEFF_DS" );
  kdds = infile.ReadItem( kdds, "HYDR_DEP_COEFF_DS" );
  assert( kwds > 0 );
  ewds = infilePtr->ReadItem( ewds, "HYDR_WID_EXP_DS" );
  edds = infile.ReadItem( edds, "HYDR_DEP_EXP_DS" );
  ewstn = infilePtr->ReadItem( ewstn, "HYDR_WID_EXP_STN" );
  edstn = infile.ReadItem( ewstn, "HYDR_DEP_EXP_STN" );
  knds = infilePtr->ReadItem( knds, "HYDR_ROUGH_COEFF_DS" );
  assert( knds > 0 );
  ends = infilePtr->ReadItem( ends, "HYDR_ROUGH_EXP_DS" );
  enstn = infilePtr->ReadItem( enstn, "HYDR_ROUGH_EXP_STN" );
  klambda = infilePtr->ReadItem( klambda, "BANK_ROUGH_COEFF" );
  elambda = infilePtr->ReadItem( elambda, "BANK_ROUGH_EXP" );
  dscrtwids = infilePtr->ReadItem( dscrtwids, "DEF_CHAN_DISCR" );
  assert( dscrtwids > 0 );
  allowfrac = infilePtr->ReadItem( allowfrac, "FRAC_WID_MOVE" );
  assert( allowfrac > 0 );
  leavefrac = infilePtr->ReadItem( leavefrac, "FRAC_WID_ADD" );
  assert( leavefrac > 0 );
  //vegerod = infile.ReadItem( vegerod, "VEG_ERODY" );
  //rockerod = infile.ReadItem( rockerod, "KB" );
  rockerod = infile.ReadItem( rockerod, "BANKERO" );
  //double MB = infile.ReadItem( MB, "MB" ); // bug fix 3/16/99
  //rockerod *= .05*pow(SECPERYEAR,MB);
  //vegerod *= .05*pow(SECPERYEAR,MB);  // added gt 3/15/99
  //latadjust = infile.ReadItem( latadjust, "LATADJUST" );
  //double shrcoeff = 1.0 / ( 1000.0 * GRAV * pow( knds / kwds, 0.6 ) );
  //vegerod *= shrcoeff * latadjust;
  //rockerod *= shrcoeff * latadjust;

  //find dependence of bank erody on bank height, P, 0<=P<=1:
  Pdz = infile.ReadItem( Pdz, "BNKHTDEP" );
  //MakeReaches();
  assert( &reachList != 0 );
}

tStreamMeander::~tStreamMeander()
{
  //if( netPtr != 0 ) delete netPtr;
  meshPtr = 0;
  netPtr = 0;
  infilePtr = 0;
  rand = 0;
}


/*****************************************************************************\
**
**       FindMeander: simply goes through active mesh nodes and assigns
**                    meander status flag based on a threshold discharge
**
**
**		  Parameters: critflow
**      Data members updated: tLNode->chan.migration.meander
**      Called by: MakeReaches
**      Calls: no "major" functions
**
**      Created: SL 1/98
**
**
\*****************************************************************************/
void tStreamMeander::FindMeander()
{
  tLNode * cn;
  tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );

  if (0) //DEBUG
    std::cout << "FindMeander...()";

  for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
    {
      //nmg
      if (0) //DEBUG
	     std::cout<<"FM cn "<<cn->getID()<<" z = "<<cn->getZ()<<std::endl;
      cn->setReachMember( false );
#if FIXCRITFLOWBUG
      if( cn->getDrArea() >= critarea )
#else
      if( cn->getQ() >= critflow )
#endif
#undef FIXCRITFLOWBUG
	  {
	     cn->setMeanderStatus( kMeanderNode );
	     cn->setNew2DCoords( cn->getX(), cn->getY() );
	  }
      else
	  {
	     if (0) {//DEBUG
	        if( cn->Meanders() )
	           std::cout << "FindMeander: node " << cn->getID() << " has Q " << cn->getQ() << " and A " << cn->getDrArea()
		            << " and mndr is being switched off." << std::endl;
	     }
	     cn->setMeanderStatus( kNonMeanderNode );
	   }
    }
  if (0) //DEBUG
    std::cout << "done\n";
}


tLNode* tStreamMeander::BlockShortcut( tLNode* crn, tLNode* bpn, tLNode& nn,
                                       const tArray<double>& ic, double time )
{
   // record a pointer to the present downstream neighbor
   // (don't use flowedge!):
   tLNode* dn = crn->getDownstrmNbr();

   // already assertained that bpn flows to nPtr;
   // called from InterpChannel, so crn doesn't flow to bpn
   // and it needs to be changed
   // make crn flow to bpn
   crn->setDownstrmNbr( bpn );
   //meander-ize the bypassed node:
   bpn->setMeanderStatus( kMeanderNode );
   bpn->setNew2DCoords( bpn->getX(), bpn->getY() );
   bpn->AddDrArea( crn->getDrArea() );
   // find coords for bank node to add
   tArray< double > dumcoords(4);
   dumcoords[0] = ic[0];
   dumcoords[1] = ic[1];
   tArray< double > checkZOld = bpn->getZOld();
   // find coordinates of a new bank node to block
   // the shortcut; coords will be outside channel,
   // so no need to check InChannel.
   if( checkZOld[0] != 0.0 || checkZOld[1] != 0.0 )
   {
      FindBankCoords( bpn, dumcoords );
      nn = *bpn;
      nn.setX( ic[0] );
      nn.setY( ic[1] );
      const tArray< double > zeroArr(4);
      nn.setXYZD( zeroArr );
      // TMP HACK? add a centimeter-SL 9/03
      if( dumcoords[2] < bpn->getZ()+0.01) dumcoords[2] = bpn->getZ()+0.01;
   }
   else{
      FindBankCoords( crn, dumcoords );
      // TMP HACK? add a centimeter-SL 9/03
      if( dumcoords[2] < crn->getZ()+0.01) dumcoords[2] = crn->getZ()+0.01;
   }
   nn.setZ( dumcoords[2] );
   nn.setMeanderStatus( kNonMeanderNode );
   nn.setDrArea( 0.0 );
   tLNode* newnodeP = meshPtr->AddNode( nn, kNoUpdateMesh, time );
   // if addition should fail, unmeanderize bpn
   if( newnodeP == NULL ){
      bpn->setMeanderStatus( kNonMeanderNode );
      bpn->setNew2DCoords( 0.0, 0.0 );
      bpn->AddDrArea( -crn->getDrArea() );
      crn->setDownstrmNbr( dn );
	  if(0) //DEBUG
	     std::cout<<"BlockShortcut: addition failed, un-meanderizing node "<<bpn->getID()<<std::endl;
   }
   
   if(1 && newnodeP!=NULL) //DEBUG
   {
      std::cout<<"BlockShortcut: added node "<<newnodeP->getID()<<" at "<<newnodeP->getX()<<","<<newnodeP->getY()
	     <<","<<newnodeP->getZ()<<" crn="<<crn->getID()<<" bpn="<<bpn->getID()<<std::endl;
   }
   
   return newnodeP;
}

tLNode* tStreamMeander::FindShortcutNode( tLNode* crn, tLNode* nPtr )
{
   tSpkIter crnI( crn );
   bool shortcut = false;
   tLNode* bpn = NULL;
   tEdge* ce = NULL;
   for( ce = crnI.FirstP(); !crnI.AtEnd(); ce = crnI.NextP() )
       if( ce->getDestinationPtr() != nPtr &&
           ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary &&
           ce->getDestinationPtr()->getZ() < crn->getZ() ){
          bpn = static_cast< tLNode* >( ce->getDestinationPtrNC() );
          if( bpn->getFlowEdg()->getDestinationPtrNC() == nPtr )
          {
             shortcut = true;
             break;
          }
       }
   if( shortcut ) return bpn;
   return NULL;
}

//use linear interpolation with noise added
//(for Delaunay mesh, avoid precisely co-linear points):
const tArray< double > tStreamMeander::FindInterpCoords( tLNode* crn, tLNode* nPtr )
{
   const double curseglen = crn->getFlowEdg()->getLength();
   const double defseglen = dscrtwids * crn->getChanWidth();
   // 3.0 originally, any physical basis for this ? (QC) Not really -SL.
   const double bigseglen = 3.0 * defseglen;
   //if we can get the desired spacing by adding a single node
   tArray< double > tic;
   if( curseglen <= bigseglen )
   {
      const double x0 = crn->getX();
      const double y0 = crn->getY();
      const double z0 = crn->getZ();
      const double x1 = nPtr->getX();
      const double y1 = nPtr->getY();
      const double ydif = y1 - y0;
      const double xdif = x1 - x0;
      const double phi = atan2( ydif, xdif );
      const double val1 = rand->ran3() - 0.5;
      const double x = x0 + curseglen * cos(phi) /
          2.0 + 0.001 * val1 * defseglen;//pow(defseglen, 2.0);
      const double val2 = rand->ran3() - 0.5;
      const double y = y0 + curseglen * sin(phi) /
          2.0 + 0.001 * val2 * defseglen;//pow(defseglen, 2.0);
      //const double z = z0 - curseglen / 2.0 * slope; // Removed by quintijn
      const double z = (nPtr->getZ() + z0)/2.0;
      tArray< double > ttic( x, y, z );
      tic = ttic;
   }
   //otherwise, if we need to add more than one point,
   //generate a random walk with uniform spacing in x
   //and exponentially weighted steps in y
   //(weighting to keep last step, back to the orig. downstream node,
   //from being too much of a doozy):
   else
   {
      const int npts = ROUND( curseglen / defseglen );
      tArray< double > xp( npts), yp(npts), zp(npts);
      tArray< double > ttic( 3*(npts-1) );
      const double x0 = crn->getX();
      const double y0 = crn->getY();
      const double z0 = crn->getZ();
      const double x1 = nPtr->getX();
      const double y1 = nPtr->getY();
      const double y = y1 - y0;
      const double x = x1 - x0;
      const double phi = atan2( y, x );
      const double slope = crn->calcSlope();
      for(int i=1; i<npts; i++ )
      {
         xp[i] = static_cast<double>(i) * defseglen;
         const double val = rand->ran3() - 0.5;
         yp[i] = yp[i-1] + 0.001 * val * exp(-yp[i-1] * val)
             * defseglen;//pow(defseglen, 2.0);
         zp[i] = xp[i] * slope;
         const double dd = sqrt(xp[i] * xp[i] + yp[i] * yp[i]);
         const double da = atan2(yp[i], xp[i]);
         const double x = dd * cos(da + phi) + x0;
         const double y = dd * sin(da + phi) + y0;
         const double z = z0 - zp[i];
         ttic[ (i-1)*3 ] = x;
         ttic[ (i-1)*3 + 1 ] = y;
         ttic[ (i-1)*3 + 2 ] = z;
      }
      tic = ttic;
   }
   const tArray< double > ic( tic );
   return ic;
}

/****************************************************************\
**
**	InterpChannel: Fill in points between widely spaced
**			channel nodes with noisy interpolation (i.e., don't make
**      perfectly straight lines). If only one point is added, it's
**      added in between the two channel nodes, colinear except for
**      a small random perturbation. If more than one node is added
**      between a channel pair, then the routine generates a random
**      walk in which the step size is weighted by an exponential
**      which decays as the distance from the line connecting the
**      two channel points increases. This weighting tends to keep
**      the final point of the walk within reasonable bounds.
**
**      Uses reachList.
**
**
**		Parameters: dscrtwids -- default spacing in # hydr widths
**		Called by:  MakeReaches
**    Calls: tMesh::AddNode, tStreamNet::UpdateNet
**		Created:  5/16/97 SL
**    Updated: 1/98 SL
**
\****************************************************************/
int tStreamMeander::InterpChannel( double time )
{
   const double timetrack = time;
   if (0) {//DEBUG
      if( timetrack >= kBugTime )
          std::cout << "InterpChannel()\n";
   }
   const tArray< double > zeroArr(4);
   double slope=0.;
   bool change = false; //haven't added any nodes yet
   bool needmeshupdate = false;
   //loop through reaches:
   int i;
   tPtrList< tLNode > *creach;
   for( creach = rlIter.FirstP(), i=0; !rlIter.AtEnd(); creach = rlIter.NextP(), i++ )
   {
      tPtrListIter< tLNode > rnIter( *creach );
      //loop through reach nodes:
      int j;
      tLNode *crn;
      for( crn = rnIter.FirstP(), j=0; j<nrnodes[i]; crn = rnIter.NextP(), ++j )
      {
         tLNode *nPtr = crn->getDownstrmNbr();
         //mesh changes could make this wrong unless recalculated:
         const double curseglen = crn->getFlowEdg()->CalcLength();
         assert( curseglen > 0 );
         //for negative slope, bail out of
         //these loops through the nodes and call UpdateNet():
         slope = crn->calcSlope();
         if( slope < 0.0 )
         {
            change = true;
            break;
         }
         const double curwidth = crn->getChanWidth();
         const double defseglen = dscrtwids * curwidth;
         const double maxseglen = 2.0 * defseglen;
         //if current flowedg length is greater than
         //a given proportion of hydraulic width
         if( curseglen > maxseglen )
         {
            needmeshupdate = true;//flag so we know that we need to update the mesh
            tLNode nn = *crn; //added nodes are copies of the upstream node except xyz.
            nn.setXYZD( zeroArr ); //and xyzd ('old' coords)
            const tArray< double > ic = FindInterpCoords( crn, nPtr );
            // see if we think we need to interpolate because a single
            // node has been shortcut:
            tLNode* bpn = FindShortcutNode( crn, nPtr );
            // also check whether "ZOld" (elevations at banks) have been set;
            // their being set indicates that meandering has been done
            // already at least once and this is not the first interpolation.
            const tArray< double > checkZOld = crn->getZOld();
            if( ( checkZOld[0] != 0.0 || checkZOld[1] != 0.0 )
                && bpn != NULL ){
               // add node to block shortcut; also changes flow directions.
               tLNode* newnodeP = BlockShortcut( crn, bpn, nn, ic, time );
               if( newnodeP != NULL ){
                  change = true;
                  if (1) //DEBUG
                      std::cout<<"IC BS pt "<<newnodeP->getID()<<" added at "
                          << newnodeP->getX() << "," << newnodeP->getY() << std::endl;
               } else {
		 if (1) //DEBUG
		   std::cout<<"IC BS pt NOT added at "
		       << ic[0] << "," << ic[1] << std::endl;
		 // Process of aborting point addition may have left node(s)
		 // without valid flowedges:
		 netPtr->ReInitFlowDirs();
		 change = true;
               }
            }
            else{
               const int numadd = ic.getSize() / 3;
               tLNode *prevNode = crn;
               tLNode* newnodeP = NULL;
               for( int k=0; k<numadd; ++k ){
                  nn.set3DCoords( ic[k*3+0], ic[k*3+1], ic[k*3+2] );
                  // set flow edge temporarily to zero, so that it is flippable
                  prevNode->setFlowEdgToZero();
                  newnodeP = meshPtr->AddNode( nn, kNoUpdateMesh, time );
                  // did a node get added?
                  // need to reset flow directions here
                  if( newnodeP != NULL )
                  {
                     change = true;
                      if (1) //DEBUG
			std::cout<<"IC pt "<<newnodeP->getID()<<" added at "
			    << ic[k*3+0] << "," << ic[k*3+1] << std::endl;
                    // previous node (prevNode) flows to new node (newnodeP)
                     prevNode->setDownstrmNbr( newnodeP );
                     // (check for rare case when nodes not connected)
                     if( prevNode->getFlowEdg() == NULL )
                         meshPtr->ForceFlow( prevNode, newnodeP, time );
                     // Paranoia
                     assert( prevNode->getFlowEdg()->getDestinationPtr() == newnodeP );
                     prevNode = newnodeP;
                  } else {
		    if (1) //DEBUG
		      std::cout<<"IC pt NOT added at "
			  << ic[k*3+0] << "," << ic[k*3+1] << std::endl;
		    // Process of aborting point addition may have left node(s)
		    // without valid flowedges:
		    netPtr->ReInitFlowDirs();
		    change = true;
                  }
               }
               // after adding all new nodes, make last new node (newnodeP) flow to
               // the next meander node downstream (nPtr)
               prevNode->setDownstrmNbr( nPtr );
               // (check for rare case when nodes not connected)
               if ( prevNode->getFlowEdg() == NULL ) {
		 change = true;
		 meshPtr->ForceFlow( prevNode, nPtr, time );
	       }
               // Paranoia
               assert( newnodeP->getFlowEdg()->getDestinationPtr() == nPtr );
            }
         }
      }
      if( slope < 0.0 ) break; //'change' has already been set to 1
   }
   if( needmeshupdate )
       meshPtr->UpdateMesh();
   if( change )
      netPtr->UpdateNet( time );
   return( ( change ) ? 1 : 0 );
}//End InterpChannel

/*****************************************************************************\
**
**      MakeReaches : Does everything needed to make meandering reaches.
**                    Takes care of the loop of
**                    1) find meandering nodes
**                    2) make reaches out of meandering nodes
**                    3) go through reaches and add points
**                       between nodes which are too far apart
**                    4) Update the areas, stream net, etc.
**                    5) go back to #1 if any points were added
**
**
**      Data members updated: reachList, nrnodes, reachlen, taillen
**      Called by: tStreamMeander(...) constructor and
**                 tStreamMeander::Migrate()
**      Calls: FindMeander, FindReaches,
**             InterpChannel=>calls netPtr->UpdateNet() if points are added
**      Receives: ctime = current time
**
**      Created:  1/98 SL
**      Modified:
**
**
\*****************************************************************************/

void tStreamMeander::MakeReaches( double ctime)
{
  if (0) //DEBUG
    std::cout<<"tStreamMeander::MakeReaches...";
  netPtr->UpdateNet( ctime ); //first update the net
  do
    {
      FindMeander(); //if Q > Qcrit, meander = TRUE
      netPtr->FindChanGeom();//sets meander=false if slope <=0
      netPtr->FindHydrGeom();
      FindReaches(); //find reaches of meandering nodes
    }
  while( InterpChannel( ctime ) ); //updates
  if (0) //DEBUG
    std::cout<<"done MakeReaches"<<std::endl;
}


/*****************************************************************************\
**
**  tStreamMeander::FindReaches
**
**  Constructs the reach objects; does not interpolate channel:
**         1) finds channel head nodes, i.e., furthest upstream extent of
**            channel that meanders, and puts each head at the beginning of
**            a ptr list of meandering nodes; these lists are the reaches;
**            adds each one-member list to a list of reaches;
**         2) goes downstream from heads and adds a node to a reach list if it
**            is not already a reach member, i.e., has not been added to
**            another reach list; when a node is added to the list, its
**            'reachmember' flag is set to keep it from being added to another
**            reach;
**         3) because the meandering "forces" are distal, i.e., they act on
**            points downstream, we tack a "tail" onto each reach; tails have
**            arbitrary length, some number * the channel width at the last
**            reach node;
**         4) during the above, keep track of the number of nodes in each
**            reach, the reach length, and the tail length; these are kept
**            in arrays which are dimensioned after we know how many reaches
**            there are.
**
**
**       **Might be good to put in something that will excise     **
**       **very sharp bends, i.e., bank undulations, from reaches.**
**
**
**      Data members updated: reachList, nrnodes, reachlen, taillen
**      Called by: MakeReaches
**      Calls: FindChanGeom, FindHydrGeom
**
**      Created:  1/98 SL
**      Modified:
**
**
\*****************************************************************************/
void tStreamMeander::FindReaches()
{
   double curwidth, ctaillen;
   int i;
   double rdrop;
   tLNode *cn, *frn, *lrn, *nxn, *ln;
   tMesh< tLNode >::nodeListIter_t nodIter( meshPtr->getNodeList() );
   tPtrListIter< tLNode > rnIter;
   tPtrList< tLNode > rnodList, *plPtr, listtodelete;
   rlListNode_t *tempnode;

   if (0) //DEBUG
       std::cout<<"tStreamMeander::FindReaches()"<<std::endl;

   if( !(reachList.isEmpty()) ) reachList.Flush();
   //loop through active nodes
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      //if node meanders
      if( cn->Meanders() && getUpstreamMeander( cn ) == NULL )
      {
         //if node meanders and has no upstream meandering nbr, then it's a reach "head";
         //put it at the beginning of an otherwise empty reach
         //and add the reach to the list of reaches:
         rnodList.Flush();
         rnodList.insertAtFront( cn );
         if (0){ //DEBUG
            if(cn->calcSlope()<=0) std::cout<<"bad node "<<cn->getID()<<" with slope "<<cn->calcSlope()<<" added to reachlist"<<std::endl;
         }
         reachList.insertAtBack( rnodList );
         cn->setReachMember( true );
      }
   }

   if( reachList.getSize() == 0 ) return;
   {
      tArray< int > *iArrPtr;
      iArrPtr = new tArray< int >( reachList.getSize() );
      nrnodes = *iArrPtr;
      delete iArrPtr;
   }
   {
      tArray< double > *fArrPtr;
      fArrPtr = new tArray< double >( reachList.getSize() );
      reachlen = *fArrPtr;
      taillen = *fArrPtr;
      delete fArrPtr;
   }
   if (0) //DEBUG
       std::cout << "No. reaches: " << reachList.getSize() << std::endl;
   //loop through reaches
   //rlIter is the iterator for reachlist with is a list of ptrs to node lists
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
   {
      assert( reachList.getSize() > 0 );
      if (0) //DEBUG
          std::cout << " on reach " << i << std::endl;
      assert( i<reachList.getSize() );
      //go downstream from reach head
      //and add nodes to the reach
      //if they're not already reach members:
      rnIter.Reset( *plPtr );
      cn = rnIter.FirstP();
      nrnodes[i]++;
      reachlen[i] += cn->getFlowEdg()->getLength();
      cn = cn->getDownstrmNbr();
      while( cn->getBoundaryFlag() == kNonBoundary
             && !cn->getReachMember()
             && cn->Meanders() )
      {
         //assert( cn->GetFlowEdg()->getLength() > 0 );
         if (0){ //DEBUG
            if(cn->calcSlope()<=0) std::cout<<"bad node "<<cn->getID()<<" with slope "<<cn->calcSlope()<<" added to reachlist"<<std::endl;
         }
         nrnodes[i]++;
         reachlen[i] += cn->getFlowEdg()->getLength();
         cn->setReachMember( true );
         plPtr->insertAtBack( cn );
         cn = cn->getDownstrmNbr();
      }

      //make sure reach has more than 4 members:
      if (0) //DEBUG
          std::cout << "reach " << i << " length " << plPtr->getSize() << std::endl;

      if( plPtr->getSize() > 4 )
      {
         //make sure reach has "positive" slope:
         lrn = rnIter.LastP();
         frn = rnIter.FirstP();
         rdrop = frn->getZ() - lrn->getZ();
         if (0) //DEBUG
             std::cout << "reach drop: " << rdrop << std::endl;
         if( rdrop <= 0 )
         {
            //remove reach if it has non-positive slope:
            if (0) //DEBUG
                std::cout << "remove reach for non-positive slope" << std::endl;
            tempnode = rlIter.NodePtr();
            rlIter.Prev();
            reachList.moveToFront( tempnode );
            reachList.removeFromFront( listtodelete );
            reachlen[i] = 0;
            nrnodes[i] = 0;
            i--;
         }
      }
      else
      {
         //remove reach if it has 4 or fewer members:
         if (0) //DEBUG
             std::cout << "remove reach for being too small" << std::endl;
         tempnode = rlIter.NodePtr();
         rlIter.Prev();
         reachList.moveToFront( tempnode );
         reachList.removeFromFront( listtodelete );
         reachlen[i] = 0;
         nrnodes[i] = 0;
         i--;
      }
   }
   if (0){ //DEBUG
      std::cout << "Final no. reaches: " << reachList.getSize() << std::endl;
      for( /*plPtr =*/ rlIter.FirstP(), i=0; !(rlIter.AtEnd());
                       /*plPtr =*/ rlIter.NextP(), i++ )
      {
         std::cout << "reach " << i << " length " << nrnodes[i] << std::endl;
      }
   }
   //construct a tail for the reach
   //which is some number of hydr. widths long:
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
   {
      if (0) //DEBUG
          std::cout << " on reach " << i << std::endl;
      assert( i<reachList.getSize() );
      rnIter.Reset( *plPtr );
      cn = rnIter.LastP();
      curwidth = cn->getChanWidth();
      taillen[i] = 10.0*curwidth;
      ctaillen = 0.0;
      cn = cn->getDownstrmNbr();
      //ng added a check to make sure slope > 0 or else code crashes
      while (ctaillen <= taillen[i]
             && cn->getBoundaryFlag() == kNonBoundary
             && cn->getChanSlope()>0 )
      {
         //assert( cn->GetFlowEdg()->getLength() > 0 );
         ctaillen+= cn->getFlowEdg()->getLength();
         plPtr->insertAtBack( cn );
         cn = cn->getDownstrmNbr();
      }
      taillen[i] = ctaillen;
   }
   //check reaches for correct connectivity:
   for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        plPtr = rlIter.NextP(), i++ )
   {
      rnIter.Reset( *plPtr );
      //ln = rnIter.LastP();
      int j;
      for( cn = rnIter.FirstP(), j=0; !(rnIter.AtEnd()); cn = rnIter.NextP(), j++ )
      {
         if( j < plPtr->getSize() - 1 )  //if( cn != ln )
         {
            nxn = rnIter.ReportNextP();
            if( nxn != cn->getDownstrmNbr() )
            {
               std::cout << "reach with " << plPtr->getSize() << " nodes: " << std::endl;
               for( ln = rnIter.FirstP(); !(rnIter.AtEnd()); ln = rnIter.NextP() )
               {
                  std::cout << ln->getID() << "; ";
               }
               std::cout << std::endl;
               cn->TellAll();
               nxn->TellAll();
               ReportFatalError( "FindReaches: downstream nbr not next reach node" );
            }
         }
      }
   }

   if (0){ //DEBUG
      for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
      {
         std::cout << "end FindReaches, node " << cn->getID() << std::endl;
      }
   }
   if (0) //DEBUG
       std::cout << "done FindReaches" << std::endl;
}


/***********************************************************************\
**
**	tStreamMeander::CalcMigration
**
**  This is the routine that makes the arrays to pass to the fortran
**	routine 'meander_' and sets "new" x and y values;
**  makes a bunch of tArrays, and passes the arrays to
**  fortran by way of the tArray member ptr, gotten
**  with getArrayPtr().
**
**  Also sets "old" x and y to present node coords if they
**  have not been set already for meandering nodes. Leaves
**  old z and channel side flag unset (=0) to indicate that
**  these are not the coords at which a point will be dropped.
**  Rather, these coords are used to determine how far the
**  meandering node has moved since it last dropped a node.
**
**		Parameters:	allowfrac -- fraction of chanwidth
**				          a node is allowed to move in a
**				          given meander iteration
**				        duration -- storm duration,
**                  copy from netPtr->stormPtr->stdur
**                time -- running tab on how long we've meandered
**    NOTE: Nicole made a change so that time is now current time
**          and duration = time at begining of migration call + stdur
**          This shouldn't theoretically affect things, but be aware.
**
**		Called by:	Migrate
**    Calls: FindBankErody, external fortran routine _meander
**		Created: 5/1/97  SL
**    Modified: 11/03 SL, no longer sets "old" x and y.
**
\***************************************************************/
void tStreamMeander::CalcMigration( double &time, double const &duration,
                                    double &cummvmt )
{
   if (0) //DEBUG
       std::cout<<"tStreamMeander::CalcMigration()...";

   //loop through reaches:
   int i;
   tPtrList< tLNode > *creach;
   for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        creach = rlIter.NextP(), i++ )
   {
      tPtrListIter< tLNode > rnIter;
      rnIter.Reset( *creach );
      {
         const int num = nrnodes[i];
         int j;
         tLNode *curnode;
         for( curnode = rnIter.FirstP(), j=0; j<num;
              curnode = rnIter.NextP(), j++ )
             //initialize deltax, deltay, newx, newy:
         {
            curnode->setLatDisplace( 0.0, 0.0 );
            curnode->setNew2DCoords( curnode->getX(), curnode->getY() );
            if (0){ //DEBUG
               const tArray< double > newxy = curnode->getNew2DCoords();
               std::cout << "init. new coords to " << newxy[0] << " " << newxy[1] << std::endl;
            }
         }
      }
      if (0) //DEBUG
          std::cout << "reach " << i << " length " << nrnodes[i] << std::endl;

      const int stations  // number of actual landscape nodes on reach
          = nrnodes[i];
      // total # of reach nodes including 'tail'
      const int nttlnodes = creach->getSize();

      tArray< double >
          xa(nttlnodes), ya(nttlnodes), xsa(nttlnodes), qa(nttlnodes),
          rerodya(nttlnodes), lerodya(nttlnodes), delsa(nttlnodes),
          slopea(nttlnodes), widtha(nttlnodes), deptha(nttlnodes),
          diama(nttlnodes), deltaxa(nttlnodes), deltaya(nttlnodes),
          rdeptha(nttlnodes), ldeptha(nttlnodes), lambdaa(nttlnodes);

      {
         double xs = 0.0;
         tLNode *curnode;
         int j;
         for( curnode = rnIter.FirstP(), j=0; !(rnIter.AtEnd());
              curnode = rnIter.NextP(), j++ )
         {
            // Set up the coordinate, streamwise length & distance, and Q arrays
            tEdge *fedg = curnode->getFlowEdg();
            xa[j] = curnode->getX();
            ya[j] = curnode->getY();
            xsa[j] = xs;
            qa[j] = curnode->getQ()/SECPERYEAR;
            delsa[j] = fedg->getLength();
            xs += delsa[j];

            // For debugging: make sure the next node on the reach is the
            // downstream neighbor of the current node
            if( j < creach->getSize() - 1 )
            {
               tLNode *nxtnode = rnIter.ReportNextP();
               if( nxtnode != curnode->getDownstrmNbr() )
               {
                  curnode->TellAll();
                  nxtnode->TellAll();
                  ReportFatalError( "downstream nbr not next reach node" );
               }
            }

            // Set bank erodibility on left and right banks
            const tArray< double > bankerody = FindBankErody( curnode );
            lerodya[j] = bankerody[0];
            rerodya[j] = bankerody[1];

            // Set slope, width, depth, grainsize, and roughness arrays
            slopea[j] = curnode->getHydrSlope();
            widtha[j] = curnode->getHydrWidth();
            deptha[j] = curnode->getHydrDepth();
            if (0) //DEBUG
                std::cout << "width, depth " << widtha[j] << " " << deptha[j] << std::endl;
            diama[j] = ( optdiamvar ) ? curnode->getDiam() : meddiam;
            lambdaa[j] = curnode->getBankRough();
         }
      }

      // Now we pass all this information to meander.f
      if (0) //DEBUG
          std::cout << "stations, stnserod: " << stations <<" "<< nttlnodes
               << std::endl;
      //this looks horrible, but we need to pass the pointer to the array
      //itself, not the tArray object, which contains the array pointer.
      meander_( &stations,
                &nttlnodes,
                xa.getArrayPtr(),
                ya.getArrayPtr(),
                xsa.getArrayPtr(),
                delsa.getArrayPtr(),
                qa.getArrayPtr(),
                rerodya.getArrayPtr(),
                lerodya.getArrayPtr(),
                slopea.getArrayPtr(),
                widtha.getArrayPtr(),
                deptha.getArrayPtr(),
                diama.getArrayPtr(),
                deltaxa.getArrayPtr(),
                deltaya.getArrayPtr(),
                rdeptha.getArrayPtr(),
                ldeptha.getArrayPtr(),
                lambdaa.getArrayPtr() );

      // Now reset the node values according to the arrays:
      {
         tLNode *curnode;
         int j;
         for( curnode = rnIter.FirstP(), j=0; !(rnIter.AtEnd());
              curnode = rnIter.NextP(), j++ )
         {
            curnode->addLatDisplace( deltaxa[j], deltaya[j] );
         }
      }
      //arbitrary change for simple debugging:
      //( 0.1, 0.01*(rand->ran3()-.5)
      //go through and set the elevations at the right and left banks:
      {
         const int num = nrnodes[i];
         double dbg2=0;
         tLNode *curnode;
         int j;
         for( curnode = rnIter.FirstP(), j=0; j<num;
              curnode = rnIter.NextP(), j++ )
         {
#define POINTBARFIX 1
#if POINTBARFIX
            const double
			    bankfullDepth = curnode->getChanDepth(),
                rz = curnode->getZ() + bankfullDepth * (1.0 - rdeptha[j]/deptha[j]),
                lz = curnode->getZ() + bankfullDepth * (1.0 - ldeptha[j]/deptha[j]);
#else
            const double
                rz = curnode->getZ() + deptha[j] - rdeptha[j],
                lz = curnode->getZ() + deptha[j] - ldeptha[j];
#endif
#undef POINTBARFIX
            //rz = curnode->getZ() + deptha[j];
            //lz = rz + deptha[j];
            dbg2 = dbg2 + (rdeptha[j]-ldeptha[j]);
            curnode->setZOld( rz, lz );
			if(0) //DEBUG
			   std::cout<<"CalcMigration: Flood depth "<<deptha[j]<<" Bank height above chan "<<rz-curnode->getZ()<<" "<<lz-curnode->getZ()<<std::endl;
         }
         if (0) //DEBUG
             std::cout << "MEAN rldepth " << dbg2 << std::endl;
      }
   }
   //calculate ratio of total displacement to length of flow edge,
   //find maximum of ratios and make sure it is less than or equal,
   //to the allowed fraction (e.g., 1/10) and scale displacements
   //as necessary.
   double maxfrac = 0.0;
   for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        creach = rlIter.NextP(), i++ )
   {
      tPtrListIter< tLNode > rnIter;
      rnIter.Reset( *creach );
      const int num = nrnodes[i];
      tLNode *curnode;
      int j;
      for( curnode = rnIter.FirstP(), j=0; j<num;
           curnode = rnIter.NextP(), j++ )
      {
         const tArray< double > delta = curnode->getLatDisplace();
         const double displcmt =
             sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
         const double width = curnode->getChanWidth();
         const double frac = ( width > 0.0 ) ? displcmt / width: 0.;
         if( frac > maxfrac ) maxfrac = frac;
      }
   }
   //maximize the time step:
   //the smaller of the time to move the allowed distance...
   double dtm = ( maxfrac > 0 ) ? allowfrac / maxfrac : 1.0;
   for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        creach = rlIter.NextP(), i++ )
   {
      tPtrListIter< tLNode > rnIter;
      rnIter.Reset( *creach );
      const int num = nrnodes[i];
      tLNode *curnode;
      int j;
      for( curnode = rnIter.FirstP(), j=0; j<num;
           curnode = rnIter.NextP(), j++ )
      {
         const double tmptim = time + dtm;
         //...and the time remaining in the storm
         if( tmptim > duration ) dtm = duration - time;
         tArray< double > delta = curnode->getLatDisplace();
         delta[0] *= dtm;
         delta[1] *= dtm;
         if( curnode == netPtr->getInletNodePtr() )
         {
            delta[0] = 0;
            delta[1] = 0;
         }
         curnode->setLatDisplace( delta[0], delta[1] );
         tArray< double > newxy = curnode->getNew2DCoords();
         newxy[0] += delta[0];
         newxy[1] += delta[1];
         curnode->setNew2DCoords( newxy[0], newxy[1] );
         if (0) //DEBUG
             std::cout << "new coords set to " << newxy[0] << " " << newxy[1] << std::endl;
      }
   }
   time += dtm;
   cummvmt += maxfrac * dtm;

   if (0) //DEBUG
       std::cout<<"done CalcMigration\n";

}


/***********************************************************************\
**
**  Migrate: The "master" meandering routine, it calls all of the
**    routines which "do" the meandering
**		Parameters: storm duration
**		Called by: Main
**    Receives : ctime = current time
**    Calls: MakeReaches--calls all routines necessary to make reaches
**           CalcMigration--calls meander model, sets newx, newy
**           MakeChanBorder--sets "old" coords
**           CheckBanksTooClose--deletes non-meandering nodes
**                               in the channel
**           CheckFlowedgCross--deletes nodes swept over by a
**                              migrating channel
**           CheckBrokenFlowedg--deletes points that are too close
**                               to the channel
**           tMesh::MoveNodes--changes node coords, updates mesh
**           AddChanBorder--"drops" nodes at old coords
**
**		Created: 1/98 SL
**    Modified: 9/98 NG made so that time is passed to it.
**             The time passed is now the time used in the while loop,
**             and duration is set to duration + ctime.
**             This should be OK, but NG didn't test that it was.
**      - fixed "infinite loop when no reaches found" bug, GT 3/99
**
\**********************************************************************/
void tStreamMeander::Migrate( double ctime )
{
   if (0) //DEBUG
       std::cout<<"Migrate time=" << ctime <<std::endl;
   const double duration = netPtr->getStormPtrNC()->getStormDuration()
       + ctime;
   double cummvmt = 0.0;
   //timeadjust = 86400. * pr->days;

   //NOTE: ctime and duration involve close subtraction of increasingly
   // large #'s -- shouldn't we just use the "true" duration? TODO

   while( ctime < duration )
   {
      MakeReaches( ctime ); //updates net, makes reachList, interpolates
      if( !(reachList.isEmpty()) )
      {
         if (0) //DEBUG
             std::cout<<"in loop "<<ctime<<" duration is "<<duration<<std::endl;
         CalcMigration( ctime, duration, cummvmt ); //incr time; uses reachList
         MakeChanBorder( ); //uses reachList
         CheckBndyTooClose();  //uses tMesh::nodeList
         CheckBanksTooClose(); //uses reachList
         CheckFlowedgCross(); //uses reachList
         meshPtr->MoveNodes( ctime ); //uses tMesh::nodeList
         AddChanBorder( ctime ); //uses reachList
      }
      else ctime=duration; // If no reaches, end here (GT added 3/12/99)
   }
   if (0) //DEBUG
       std::cout<<"end migrate"<<std::endl;
}

/****************************************************************\
**  tStreamMeander::FindBankCoords: Function cut from
**    MakeChanBorder.
**
**		Parameters:
**		Called by: MakeChanBorder
**		Created: 9/2003 SL
**
\***************************************************************/
void tStreamMeander::FindBankCoords( tLNode* cn, tArray< double >& posRef )
{
   assert( posRef.getSize() >= 4 );
   const tArray< double > cnpos = cn->get2DCoords();
   
   //find coordinates of banks, i.e., coords at n = + and -width/2
   tLNode* cnbr = cn->getDownstrmNbr();
   const tArray< double > dsnpos = cnbr->get2DCoords();
   double x0 = cn->getX();
   double y0 = cn->getY();
   double x1 = cnbr->getX();
   double y1 = cnbr->getY();
   double delx = x1 - x0;
   double dely = y1 - y0;
   double phi = atan2( dely, delx );
   // add a centimeter so it won't be in the channel:
   double width = cn->getChanWidth();
   double xdisp = ( 0.5 * width + 0.01 ) * sin(phi);
   double ydisp = ( 0.5 * width + 0.01 ) * cos(phi);
   tArray< double > rl = cn->getZOld();
   double z = cn->getZ();
   //find which side of the channel the old coords are at:
   if( posRef[3] == 1.0 || PointsCCW( cnpos, dsnpos, posRef ) )  // Left side
   {
      posRef[0] = x0 - xdisp;
      posRef[1] = y0 + ydisp;
      posRef[2] = ( rl[1] > z ) ? rl[1] : z;
      posRef[3] = 1.0;
      if (0) //DEBUG
          std::cout << "node " << cn->getID()
               << " old pos set, on left side of channel" << std::endl;
   }
   else  // Right side
   {
      posRef[0] = x0 + xdisp;
      posRef[1] = y0 - ydisp;
      posRef[2] = ( rl[0] > z ) ? rl[0] : z;
      posRef[3] = -1.0;
      if (0) //DEBUG
          std::cout << "node " << cn->getID()
               << " old pos set, on right side of channel" << std::endl;
   }
}

/****************************************************************************\
**
**  tStreamMeander::MakeChanBorder()
**
**  This function sets the "old coordinates" of the
**  migrating node; i.e., the coordinates at which a new
**  node will be "dropped" by a meander node. Calls ResetEffNbrCoords
**  for each reach node; that function finds the effective nbr coordinates
**  on the receding bank and sets "old" x-y to those coordinates. Then,
**  check distance from those coordinates.
**  If more than (1/2 + leavefrac) * (channel width), find xyz on channel
**  bank and set old coordinates to those xyz.
**  This should have the advantage of basing new node locations on
**  present geometry rather than on memory, as in older versions.
**
**   This is done after CalcMigration has been called
**   but before tMesh::MoveNodes is called.
**
**              Parameters:     leavefrac, chanwidth
**              Called by:      Migrate
**              Created:        8/18/97 SL
**              Updated: 11/03 SL (major change)
**
\*****************************************************************************/
void tStreamMeander::MakeChanBorder( )
{
   if (0) //DEBUG
       std::cout << "MakeChanBorder()" << std::endl;
   tPtrList< tLNode > *cr;
   int i;
   for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        cr = rlIter.NextP(), i++ )
   {
      tPtrListIter< tLNode > rnIter( cr );
      const int num = nrnodes[i];
      tLNode *cn = NULL;
      int j;
      for( cn = rnIter.FirstP(), j=0; j<num;
           cn = rnIter.NextP(), j++ )
      {
         const double width = cn->getChanWidth();
         // sets xyzd according to direction moving and
         // where nodes on opposite side are
         ResetEffNbrCoords( cn );
         //if more than a channel width from nodes on that side,
         // set bank coords for dropping a new node
         if( cn->DistFromOldXY() >= width * ( 0.5 + leavefrac ) )
         {
            if (0) //DEBUG
                std::cout << "node " << cn->getID()
                     << " >= width from receding bank" << std::endl;
            tArray< double > oldpos = cn->getXYZD();
            //find coordinates of banks, i.e., coords at n = + and -width/2
            FindBankCoords( cn, oldpos );
            cn->setXYZD( oldpos );  //coords found, update node data member
         }
      }
   }
}

/****************************************************************************\
**
**	tStreamMeander::AddChanBorder
**
**  For meandering nodes with placement coords set, check whether
**  a new node should be dropped. Checks to see whether 3D bank coords
**  were set in MakeChanBorder. If so, adds new node if it is within bounds
**  (i.e., can be found by LocateTriangle). That's it.
**
**  Parameters:
**  Called by:      Migrate
**  Created:        8/18/97 SL
**  Updated:        1/98 SL; 2/98 SL; 11/03 SL (major change)
**
\***************************************************************************/
void tStreamMeander::AddChanBorder(double time)
{
   if (0) //DEBUG
       std::cout << "AddChanBorder()" << std::endl;
   bool change = false; //haven't added any nodes yet
   const tArray< double > zeroArr(4);
   tLNode channode;

   int i;
   tPtrList< tLNode > *cr;
   for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd()); cr = rlIter.NextP(), ++i )
   {
      tPtrListIter< tLNode > rnIter( *cr );
      int j;
      tLNode* cn;
      for( cn = rnIter.FirstP(), j=0; j<nrnodes[i]; cn = rnIter.NextP(), ++j )
      {

         tArray< double > oldpos = cn->getXYZD();
         //select for nodes with old coords set:
         if( oldpos[3] != 0.0 )
         {
            if (0) //DEBUG
                std::cout << "node " << cn->getID()
                     << " ready to drop new node" << std::endl;
            //just make sure new node will be in a triangle
            tTriangle* ct = meshPtr->LocateTriangle( oldpos[0], oldpos[1] );
            if( ct != NULL )
            {
               //***NG: HERE IS WHERE YOU CAN FIND A DEPOSIT THICKNESS
               //TO ADD TO THE NEW NODE***
               channode.set3DCoords( oldpos[0], oldpos[1], oldpos[2] );
               tArray< double > xyz(3);
               for( int k=0; k<3; ++k ) xyz[k] = oldpos[k];
               // Make sure the banknode is not lower than the node it
               // originates from, bug fix 8/2003 QC. Causes ponds if the
               // meander path is redirected over the newly added banknode in FlowDir
               if( xyz[2] < cn->getZ()) xyz[2] = cn->getZ();
               channode = *cn;//added node is copy of "mother" except
               channode.set3DCoords( xyz[0], xyz[1], xyz[2] );//xyz
               channode.setXYZD( zeroArr );//and xyzd and meander and drarea
               channode.setMeanderStatus( kNonMeanderNode );
               channode.setDrArea( 0.0 );
               //TODO: NG Need to take care of deposit depth here
               //I was thinking to leave a deposit of depth
               //xyz[2]-cn->getZ() if this depth is positive
               //The texture of this deposit would be
               //the surface texture of cn.  Use erodep.
               tLNode* nnPtr = meshPtr->AddNode( channode, kNoUpdateMesh, time );
               if( nnPtr != NULL ){
		          if(0)//DEBUG
		             std::cout<<"ACB pt "<<nnPtr->getID()<<" added at "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<std::endl;
		          change = true; //flag to update mesh
               }
               cn->setXYZD( zeroArr );
            }
         }
      }
   }
   if( change )
   {
      meshPtr->UpdateMesh();
      if (0){ //DEBUG
         if( time >= kBugTime ) std::cout << "added nodes(s), AddChanBorder finished"
                                     << std::endl;
      }
   }
   
} //End AddChanBorder


/****************************************************************************\
**
**	tStreamMeander::ResetEffNbrCoords
**
**  This function finds the weighted average of appropriate bank coordinates.
**  Check all nbrs; find distances to perp. line. Find 2 pairs of
**  consecutive nbrs which fall on either side of line (going ccw from
**  flowedge, 1st pair is on left bank, 2nd pair is on right).
**  For pair on side opposite migration direction, effective bank
**  coordinates are weighted avg. of two:
**        for d1, d2 = the two distances to the perp. line; D = d1 + d2;
**        x1, x2 = x-coordinates of two pts; nbrx = effective x-coordinate
**              nbrx = [x1 * (D - d1) + x2 * (D - d2)] / D
**
**		Called by: MakeChanBorder
**		Created: 11/03 SL (modified from FindBankErody)
**
\****************************************************************************/
void tStreamMeander::ResetEffNbrCoords( tLNode *nPtr )
{
   assert( nPtr->getBoundaryFlag() == kNonBoundary );

   const tArray< double > xy1 = nPtr->get2DCoords();
   tLNode* nNPtr = nPtr->getDownstrmNbr();
   const tArray< double > xy2 = nNPtr->get2DCoords();
   const tArray< double > xy1n = nPtr->getNew2DCoords();
   const int pccw = PointsCCW( xy1, xy2, xy1n )?1:0;
   // find "effective" coordinates of bank nodes, i.e.,
   // weighted average of coordinates of channel neighbors,
   // choose the appropriate side,
   // set those coordinates in xyzd.
   tSpkIter sI(nPtr);
   tArray< double > spD( sI.getNumSpokes() * 2 ),
       spR( sI.getNumSpokes() * 2 );
   double a = xy2[0] - xy1[0];
   double b = xy2[1] - xy1[1];
   double c = -b * xy1[1] - a * xy1[0];
   const int n = sI.getNumSpokes();
   //find distance and remainders of pts wrt line perpendicular to downstream direction
   tEdge *fe = nPtr->getFlowEdg();
   tEdge *ce = fe;
   int i = 0;
   do
   {
      const tArray< double > xy = ce->getDestinationPtrNC()->get2DCoords();
      const double d = a * xy[0] + b * xy[1] + c;
      spD[i] = spD[n + i] = d;
      //if d=0 then point is on line
      spR[i] = spR[n + i] = DistanceToLine( xy[0], xy[1], a, b, c );
      ce = ce->getCCWEdg();
      ++i;
   } while( ce != fe );
   //find point pairs:
   i=0;
   int j=0;
   ce = fe;
   do
   {
      //find signs of 'this' and 'next' remainders
      double s1;
      double s2;
      if( spD[i] != 0.0 ) s1 = spD[i] / fabs( spD[i] );
      else s1 = 0.0;
      if( spD[i + 1] != 0.0 ) s2 = spD[i + 1] / fabs( spD[i + 1] );
      else
      {
         s2 = 0.0;
         spD[i + 1] = spD[i+2];
         //NG added in the above line. If spD=0, the below if will be
         //entered twice with the same node. This new line should
         //prevent that from happening.
         //Problem is that one of the surrounding nodes lies
         //exactly on perpendicular line.
      }
      if( s1 != s2 ) //points are on opposite sides of the perp.
      {
         // pccw is 1 if new coords are on left side, so want "old"
         // coords on right, or j = 1;
         // pccw is 0 if new coords are on right side, so want "old"
         // coords on left, or j = 0.
         if( j == pccw ){
            //find distances to line and their sum:
            const double d1 = spR[i];
            const double d2 = spR[i + 1];
            const double D = d1 + d2;
            assert( D > 0.0 );
            //find nbr node locations
            tEdge* ne = ce->getCCWEdg();
            tLNode* node1 = static_cast<tLNode *>( ce->getDestinationPtrNC() );
            tLNode* node2 = static_cast<tLNode *>( ne->getDestinationPtrNC() );
            // find weighted average coordinates
            const double nbrx = ( node1->getX() * d2 + node2->getX() * d1 ) / D;
            const double nbry = ( node1->getY() * d2 + node2->getY() * d1 ) / D;
            // set "old" coords
            tArray< double > oldpos(4);
            oldpos[0] = nbrx;
            oldpos[1] = nbry;
            nPtr->setXYZD( oldpos );
            break;
         }
         j++;
      }
      i++;
      ce = ce->getCCWEdg();
   } while( ce != fe );
}



/****************************************************************************\
**
**	tStreamMeander::FindBankErody
**
**  This is the routine that finds the effective erodibility of each bank
**	based on a reach node's neighbor's erodibility
**  and relative height above the channel.
**  Check all nbrs; find distances to perp. line. Find 2 pairs of
**  consecutive nbrs which fall on either side of line (going ccw from
**  flowedge, 1st pair is on left bank, 2nd pair is on right). For each
**  pair, find each pt's erod'y w.r.t. the channel pt (e.g. elev.
**  difference could be a factor) (for points lower than channel pt's
**  water surface elevation, assume level w/ surface.); erod'y for chan.
**  pt is weighted avg. of two:
**        for d1, d2 = the two distances to the perp. line; D = d1 + d2;
**        E1, E2 = erod'y of two pts; E = effective erod'y;
**              E = [E1 * (D - d1) + E2 * (D - d2)] / D
**        Bank erody may be totally, partially, or not at all dependent on
**        bank height, according to whether Pdz=1, 0<Pdz<1, or Pdz = 0:
**        E1 = E1 * (Pdz * H / dz + (1 - Pdz)), dz >= H.
**
**
**		Parameters:	tSurface::vegerody; tBedrock::erodibility
**		Called by: CalcMigration
**		Created: 5/1/97 SL
**    Modifications:
**      - removed previous call to now-obsolete function getAlluvThickness.
**        Erodibility now assumed constant and equal to rockerod. Should be
**        based on layer exposed at channel position -- TODO. (GT 3/16/99)
**
\****************************************************************************/
tArray< double >
tStreamMeander::FindBankErody( tLNode *nPtr ) const
{
   if (0) //DEBUG
       std::cout << "FBE\n";

   tArray< double > lrerody(2);
   if( nPtr->getBoundaryFlag() != kNonBoundary ) return lrerody;

   tSpkIter sI(nPtr);
   tArray< double > spD( sI.getNumSpokes() * 2 ),
       spR( sI.getNumSpokes() * 2 );

   tLNode *nNPtr = nPtr->getDownstrmNbr();
   const tArray< double > xyz1 = nPtr->get3DCoords();
   const tArray< double > xy2 = nNPtr->get2DCoords();
   tArray< double > dxy(2);
   dxy[0] = xy2[0] - xyz1[0];
   dxy[1] = xy2[1] - xyz1[1];
   double a = dxy[0];
   double b = dxy[1];
   double c = -dxy[1] * xyz1[1] - dxy[0] * xyz1[0];
   const int n = sI.getNumSpokes();
   //find distance and remainders of pts wrt line perpendicular to downstream direction
   tEdge *fe = nPtr->getFlowEdg();
   tEdge *ce = fe;
   int i = 0;
   do
   {
      tLNode *cn = static_cast<tLNode *>( ce->getDestinationPtrNC() );
      const tArray< double > xy = cn->get2DCoords();
      const double d = a * xy[0] + b * xy[1] + c;
      spD[i] = spD[n + i] = d;
      //if d=0 then point is on line
      spR[i] = spR[n + i] = DistanceToLine( xy[0], xy[1], a, b, c );
      ce = ce->getCCWEdg();
      ++i;
   } while( ce != fe );
   const double H = nPtr->getHydrDepth();
   //find point pairs:
   i=0;
   int j=0;
   ce = fe;
   do
   {
      //find signs of 'this' and 'next' remainders
      double s1;
      double s2;
      if( spD[i] != 0.0 ) s1 = spD[i] / fabs( spD[i] );
      else s1 = 0.0;
      if( spD[i + 1] != 0.0 ) s2 = spD[i + 1] / fabs( spD[i + 1] );
      else
      {
         s2 = 0.0;
         spD[i + 1] = spD[i+2];
         //NG added in the above line. If spD=0, the below if will be
         //entered twice with the same node. This new line should
         //prevent that from happening.
         //Problem is that one of the surrounding nodes lies
         //exactly on perpendicular line.
      }


      if( s1 != s2 ) //points are on opposite sides of the perp.
      {
         assert( j<2 );
         //find distances to line and their sum:
         const double d1 = spR[i];
         const double d2 = spR[i + 1];
         const double D = d1 + d2;
         //find erodibilities of nbr nodes wrt nPtr
         tEdge* ne = ce->getCCWEdg();
         tLNode* node1 = static_cast<tLNode *>( ce->getDestinationPtrNC() );
         tLNode* node2 = static_cast<tLNode *>( ne->getDestinationPtrNC() );
         //find elev. diff's:
         const double dz1 = node1->getZ() - xyz1[2];
         const double dz2 = node2->getZ() - xyz1[2];

         //find whether bedrock or alluvial bank:
         // (modified by GT 3/99: getAlluvThickness is obsolete. Layer
         // info should be used. For now assume constant erodibility. TODO)
         /*if( dz1 > node1->getAlluvThickness() ) E1 = rockerod;//node1->getBedErody();
           else E1 = vegerod;//node1->getVegErody();
           if( dz2 > node2->getAlluvThickness() ) E2 = rockerod;//node2->getBedErody();
           else E2 = vegerod;//node2->getVegErody();*/
         double E1 = rockerod, E2 = rockerod; // added 3/99

         if (0) //DEBUG
             std::cout << "E1 " << E1 << "  E2 " << E2 << std::endl;
         //find height dependence:
         //if elev diff > hydraulic depth, ratio of depth to bank height;
         //o.w., keep nominal erody:
         if (0) //DEBUG
             std::cout << "FBE 4" << H << " " << dz1 << std::endl;
         if( dz1 > H ) E1 *= (Pdz * H / dz1 + (1 - Pdz));
         if (0) //DEBUG
             std::cout << "FBE 5" << std::endl;
         if( dz2 > H ) E2 *= (Pdz * H / dz2 + (1 - Pdz));
         //now we've found erod'ies at ea. node, find weighted avg:
         assert( D > 0.0 );
         if (0) //DEBUG
             std::cout << "FBE 6" << std::endl;
         lrerody[j] = (E1 * d2 + E2 * d1) / D;
         j++;
         if (0) //DEBUG
             std::cout << "FBE 7" << std::endl;
      }
      i++;
      ce = ce->getCCWEdg();
   } while( ce != fe );
   return lrerody;
}


/*****************************************************************************\
**
**  CheckBndyTooClose():
**   Go through boundary part of node list and check for meandering nbrs. Find
**   new distance of mndr nbrs from boundary edges; if it's closer than 1/2 the
**   hydraulic width of the mndr node, do a RevertToOldCoords() on that MF IFF
**   migration would move it closer to the boundary (i.e., don't want to pin it
**   next to the boundary if it's already there).
**
**		Parameters:
**		Called by: Migrate
**		Created: 5/98 SL
**
**
\*****************************************************************************/
void tStreamMeander::CheckBndyTooClose()
{
   tMesh< tLNode >::nodeListIter_t nI( meshPtr->getNodeList() );
   tLNode *cn, *mn;
   tNode *nn, *bn0(0), *bn1(0);
   tEdge *ce;
   int n;
   double width, mindist, d0, d1, d2, d3, xp, yp;

   if (0) //DEBUG
       std::cout << "CBTC\n";

   /*cn =*/ nI.LastActiveP();
   //go through boundary nodes
   for( cn = nI.NextP(); !(nI.AtEnd()); cn = nI.NextP() )
   {
      tSpkIter sI( cn );
      n = 0;
      //count number of meandering nbrs:
      for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
      {
         mn = static_cast<tLNode *>( ce->getDestinationPtrNC() );
         if( mn->Meanders() )
         {
            n++;
         }
      }
      //proceed only if bndy node has meandering nbr(s):
      if( n > 0 )
      {
         n = 0;
         //count number of and find boundary nbrs:
         for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
         {
            nn = ce->getDestinationPtrNC();
            if( nn->getBoundaryFlag() != kNonBoundary )
            {
               if( n == 0 ) bn0 = nn;
               else if( n == 1 ) bn1 = nn;
               //else std::cout << 'Warning: >2 bndy nbrs found in CheckBndyTooClose\n';
               n++;
            }
         }
         //find meandering nbrs:
         for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
         {
            mn = static_cast<tLNode *>( ce->getDestinationPtrNC() );
            if( mn->Meanders() )
            {
               width = mn->getHydrWidth();
               mindist = width / 2.0;
               //find distances to boundary edges:
               const tArray< double > xyn = mn->getNew2DCoords();
               xp = xyn[0];
               yp = xyn[1];
               d0 = DistanceToLine( xp, yp, cn, bn0 );
               d1 = DistanceToLine( xp, yp, cn, bn1 );
               //if too close, RevertToOldCoords() on the meandering node:
               if( d0 < mindist || d1 < mindist )
               {
                  const tArray< double > xy = mn->get2DCoords();
                  xp = xy[0];
                  yp = xy[1];
                  d2 = DistanceToLine( xp, yp, cn, bn0 );
                  d3 = DistanceToLine( xp, yp, cn, bn1 );
                  if( d0 < d2 || d1 < d3 ) mn->RevertToOldCoords();
               }
            }
         }
      }
   }
   return;
}


/*****************************************************************************\
**
**  CheckBanksTooClose
**
**  This simply checks the non-meandering neighbors of all
**  meandering nodes to make sure they're not in the
**  channel segment defined by the meandering node and its
**  downstream neighbor.
**
**		Parameters:
**		Called by:
**		Created: 1/98 SL
**
\*****************************************************************************/
void tStreamMeander::CheckBanksTooClose()
{
   if (0) //DEBUG
       std::cout << "CheckBanksTooClose()..." << std::endl;


   // For each reach
   int i;
   tPtrList< tLNode > *cr;  // ptr to current reach
   for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
        cr = rlIter.NextP(), i++ )
   {
      tPtrList< tLNode > delPtrList;
      tPtrListIter< tLNode > dIter( delPtrList );
      tPtrList< tLNode > chanPtrList;
      // For each node on reach
      tPtrListIter< tLNode > rnIter( cr );   // iterator for nodes on reachrnIter.Reset( *cr );
      const int num = nrnodes[i];
      int j;
      tLNode* cn;
      for( cn = rnIter.FirstP(), j=0; j<num;
           cn = rnIter.NextP(), j++ )
      {
         // Check neighboring nodes
         tSpkIter spokIter( cn );
         for( tEdge *ce = spokIter.FirstP(); !( spokIter.AtEnd() );
              ce = spokIter.NextP() )
         {
            tLNode* sn = static_cast<tLNode *>( ce->getDestinationPtrNC() );
            //check for proximity to channel:
            if( !(sn->Meanders()) && InChannel( cn, sn ) )
            {
               // If node isn't a boundary and isn't already on the
               // deletion list, put it on the deletion list now
               if (0) //DEBUG
                   std::cout<<"too close: cn, cn->hydrwidth "<<cn->getID()<<" "
                       <<cn->getHydrWidth()<<std::endl;
               tLNode* pointtodelete = static_cast<tLNode *>( ce->getDestinationPtrNC() );
               if( pointtodelete->getBoundaryFlag() == kNonBoundary )
               {
                  if ( pointtodelete->getDrArea() < cn->getDrArea() )
                  {
                     bool onlist = false;
                     for( tLNode* dn = dIter.FirstP(); !(dIter.AtEnd());
                          dn = dIter.NextP() )
                         if( pointtodelete == dn ) onlist = true;
                     if( !onlist )
                     {
                        if (0) //DEBUG
                            std::cout << "add to delete list: "
                                 << pointtodelete->getID() << std::endl;
                        delPtrList.insertAtBack( pointtodelete );
                        // also add associated channel node to list
                        chanPtrList.insertAtBack( cn );
                     }
                  }
               }
               else cn->RevertToOldCoords();
            }
         }
      }
      // Having found all the nodes that have been "swept away" by the
      // channel and placed them on the delPtrList, we now delete them
      if( !delPtrList.isEmpty() )
      {
         for( tLNode *dn = dIter.FirstP(); !(dIter.AtEnd()); dn = dIter.FirstP())
         {
            if (0) //DEBUG
                std::cout << "CBTC: delete node " << dn->getID() << " at "
                     << dn->getX() << ", " << dn->getY() << std::endl;
            meshPtr->DeleteNode( dn, kRepairMesh, kNoUpdateMesh );
            delPtrList.removeFromFront();
            chanPtrList.removeFromFront();
         }
         meshPtr->UpdateMesh();
      }
   }
   if (0) //DEBUG
       std::cout << "finished" << std::endl;
}

/*****************************************************************************\
**
**       CheckFlowedgCross(): checks to see whether any nodes have been
**          crossed/swept over by flowedges; i.e., whether the channel has
**          eroded the node away.
**          -finds tri's with !NewTriCCW;
**          -are two vtcs. connected by a flowedg?
**          -if so, is >1 nbr tri !NewTriCCW, or does a spoke of the 3rd vtx.
**           cross the flowedg?
**          -if so, mark the 3rd node for deletion;
**          -if a point is to be deleted, is it a meandering node with greater
**           drainage area than the upstream node of the found flowedg?
**          -if so, mark the upstream node of the flowedg for deletion instead.
**          -if a point is to be deleted, is it a boundary node?
**          -if so, if either/both of the other nodes have moved outside the
**           boundaries, revert it/them to the old coords., and assert that
**           one or the other is indeed outside the boundary, o.w., exit w/err.
**          -if the point to be deleted is not a boundary node, delete it.
**
**		Parameters:
**		Called by:	Migrate
**		Created: 1/98 SL; moved/modified from tMesh routine
**    Modified, 6/2003 SL: added list of nodes to delete, as in
**      CheckBanksTooClose--was having problem with ending up "nowhere" in
**      triList.
**    Also: Did check repeatedly through list of triangles. Now check only
**    those triangles associated with meandering reach flow edges and their
**    complements.
**    TODO: Is repetition necessary?  -SL
\*****************************************************************************/
void tStreamMeander::CheckFlowedgCross()
{
   if (0) //DEBUG
       std::cout << "CheckFlowedgCross()..." << std::endl;
   tPtrList< tLNode > delPtrList;
   tPtrListIter< tLNode > dIter( delPtrList );
   //delete node crossed by flowedg:
   //if a new triangle is !CCW and two vtcs. are connected by a flowedg AND
   //  a spoke of the third vtx. intersects the flowedg OR
   //  more than one nbr. tri. is also !CCW
   //then delete third vtx. because it has been "crossed" by the flowedg
   bool crossed = true;
   do
   {
      crossed = false;
      //For each reach
      tPtrList< tLNode > *cr;
      int k;
      for( cr = rlIter.FirstP(), k=0; !(rlIter.AtEnd());
           cr = rlIter.NextP(), ++k )
      {
         // For each node on reach
         tPtrListIter< tLNode > rnIter( cr );
         tLNode* cn;
         int n;
         for( cn = rnIter.FirstP(), n=0; n<nrnodes[k]; cn = rnIter.NextP(), ++n )
         {
            int ft = 0;
            tLNode *pointtodelete = NULL;
            tLNode* nod = NULL;
            // check the triangles associated with the flowedge and its
            // complement:
            tEdge *fedg = cn->getFlowEdg();
            tTriangle *ct = fedg->TriWithEdgePtr();
            tLNode* dscn = cn->getDownstrmNbr();
            if( ct && !NewTriCCW( ct ) )
            {

               ft = 1;
               //set 'nod' to third node
               const int j = ct->nVtx( cn );
               nod = static_cast< tLNode* >( ct->pPtr( (j+2)%3 ) );
            }
            else
            {
               fedg = fedg->getComplementEdge();
               ct = fedg->TriWithEdgePtr();
               if( ct && !NewTriCCW( ct ) )
               {

                  ft = 1;
                  //set 'nod' to third node
                  const int j = ct->nVtx( cn );
                  nod = static_cast< tLNode* >( ct->pPtr( (j+1)%3 ) );
               }
            }
            //if 'flow triangle', i.e., meandering node and dnstrm nbr
            if( ft )
            {
               //check nbr tri's for new pos'n CCW
               int j = 0;
               for( int i=0; i<3; i++ )
               {
                  tTriangle* nt = ct->tPtr(i);
                  if( nt != 0 )
                      if( !NewTriCCW( nt ) )
                          ++j;
               }
               //if >1 nbr tri !CCW
               if( j>1 )
               {
                  pointtodelete = nod;
               }
               else
               {
                  tSpkIter spokIter( nod );
                  fedg = cn->getFlowEdg();
                  for( tEdge* ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                       ce = spokIter.NextP() )
                  {
                     //if flowedg crossed by spoke(s) of the third vtx.
                     if( Intersect( ce, fedg ) )
                     {
                        pointtodelete = nod;
                        break;
                     }
                  }
               }
               if( pointtodelete != 0 )
               {
                  // SL, 7/2003: wait to set "crossed" until ready to delete
                  // nodes, otherwise boundary nodes will lead to infinite loops
                  //crossed = true;
                  //if the node to delete is a meandering node
                  //with greater flow, delete the node we started with
                  // (SL, 7/2003) and remove from reach while we're here;
                  // if it doesn't have greater flow we'll put it on the
                  // list for deletion but won't remove it from reachList
                  // now; it can't get added to the deletion list twice
                  if( pointtodelete->Meanders() &&
                      cn->getDrArea() < pointtodelete->getDrArea() )
                  {
                     pointtodelete = cn;
                     if( rnIter.Prev() ){
                        --n;
                        cr->removeNext( rnIter.NodePtr() );
                     }
                     else{
                        cr->removeFromFront();
                        rnIter.First();
                     }

                     // decrement number of reach nodes:
                     --nrnodes[k];
                     // unset "reachmember"
                     pointtodelete->setReachMember(false);
                  }
                  //don't delete boundary node
                  else if( pointtodelete->getBoundaryFlag() )
                  {
                     pointtodelete = 0;
                  }
                  if( pointtodelete != 0 )
                  {
                     //meshPtr->DeleteNode( pointtodelete );
                     // 6/2003 SL: add to list of nodes to be deleted
                     // rather than delete right away:
                     bool onlist = false;

                     for( tLNode* dn = dIter.FirstP(); !(dIter.AtEnd()); dn = dIter.NextP() )
                         if( pointtodelete == dn ) onlist = true;
                     if( !onlist )
                     {
                        if (0) //DEBUG
                            std::cout << "add to delete list: "
                                 << pointtodelete->getID() << std::endl;
                        delPtrList.insertAtBack( pointtodelete );
                     }
                  }
                  else
                  {
                     //we had found a crossed node but it was on the boundary
                     //so we revert the other two to their original coords
                     //if they fall outside the mesh
                     const tArray< double > xy0 = cn->getNew2DCoords();
                     if( meshPtr->LocateTriangle( xy0[0], xy0[1] ) == 0 )
                         cn->RevertToOldCoords();
                     const tArray< double > xy1 = dscn->getNew2DCoords();
                     if( meshPtr->LocateTriangle( xy1[0], xy1[1] ) == 0 )
                         dscn->RevertToOldCoords();
                  }
               }
            }
         }
      }

      // Having found all the nodes that have been "swept away" by the
      // channel and placed them on the delPtrList, we now delete them;
      if( !delPtrList.isEmpty() )
      {
         // SL, 7/2003: set "crossed" here; otherwise crossed boundary nodes
         // can lead to infinite loops:
         crossed = true;
         for( tLNode* dn = dIter.FirstP(); !(dIter.AtEnd()); dn = dIter.FirstP() )
         {
            // rare case in which deleted nodes are not removed from reaches:
            if( dn->getReachMember() )
            {
               // didn't get deleted from reach, so we need to do it here:
               tLNode *cn = NULL;
               tPtrList< tLNode > *cr;
               int k;
               for( cr = rlIter.FirstP(), k=0; !(rlIter.AtEnd());
                    cr = rlIter.NextP(), ++k )
               {
                  tPtrListIter< tLNode > rnIter( cr );
                  cn = rnIter.GetP( dn );
                  if( cn != NULL ){
                     if( rnIter.Prev() )
                         cr->removeNext( rnIter.NodePtr() );
                     else
                     {
                        cr->removeFromFront();
                        rnIter.First();
                     }
                     // decrement number of reach nodes:
                     --nrnodes[k];
                     // unset "reachmember"
                     cn->setReachMember(false);
                     break;
                  }
               }
            }
            if (0) //DEBUG
                std::cout << "CFC: delete node " << dn->getID() << " at "
                     << dn->getX() << ", " << dn->getY() << std::endl;
            meshPtr->DeleteNode( dn, kRepairMesh, kNoUpdateMesh, true );
            delPtrList.removeFromFront();
         }
         meshPtr->UpdateMesh();
      }
   } while( crossed );
   if (0) //DEBUG
       std::cout << "finished" << std::endl;
}

/****************************************************************\
**
**   tStreamMeander::InChannel
**
**   Find whether a node (bnode) falls within an oval with
**   focii at up- and downstream nodes of channel segment,
**   minor axis of one channel width.
**
**		Parameters:	mnode -- upstream meandering channel node
**                bnode -- the potential new bank node
**		Called by:
**		Created: 1/98 SL
**    Modified:
**     - 6/99 GT: changed from hydraulic width to channel,
**       since we're interested in the morphological channel
**       not the current flood width.
**
\***************************************************************/
int tStreamMeander::InChannel( tLNode *mnode, tLNode const *bnode )
{
   //b = mnode->getHydrWidth();
   const double b = mnode->getChanWidth();  // GT changed to ChanWidth, 6/99
   if( b == 0. ) return 0;

   tMesh< tLNode >::nodeListIter_t dI( meshPtr->getNodeList() );
   tMesh< tLNode >::edgeListIter_t eI( meshPtr->getEdgeList() );
   tLNode *dnode = mnode->getDownstrmNbr();
   tEdge *fe = mnode->getFlowEdg();

   //if downstrm nbr exists and is still in nodeList; also flowedge:
   if( dnode != 0 && fe != 0 )
   {
      const double L = mnode->getFlowEdg()->getLength();
      //mindist = sqrt( L * L + b * b );
      //for perp. distance from mnode > b/2:
      const double mindist = sqrt( L * L + 0.5 * b * b
                                   + b * sqrt( L * L + 0.25 * b * b ) );
      const tArray< double >
          up = mnode->getNew2DCoords(),
          dn = dnode->getNew2DCoords(),
          bnk = bnode->get2DCoords();
      const double a = sqrt( (bnk[0] - up[0]) * (bnk[0] - up[0]) +
                             (bnk[1] - up[1]) * (bnk[1] - up[1]) );
      const double d = sqrt( (bnk[0] - dn[0]) * (bnk[0] - dn[0]) +
                             (bnk[1] - dn[1]) * (bnk[1] - dn[1]) );
      return ( a + d < mindist ) ? 1 : 0;
   }
   else return 0;
}

#undef kBugTime

/****************************************************************\
**  getUpstreamMeander: calls tLNode::getUpstrmNbr(), checks
**    whether node returned by that function is meandering
**
**		Parameters:
**		Called by:
**		Created: 9/03 SL
**
\***************************************************************/
tLNode *tStreamMeander::getUpstreamMeander(tLNode *cn)
{

   if( cn!=NULL ){
      tLNode *upstreamnode = cn->getUpstrmNbr();
      if( upstreamnode != NULL ){
         if( !upstreamnode->Meanders() ){
            if (0) //DEBUG
	      std::cout<<"upstream node non-meandering for node "<<cn->getID()<<std::endl;
            upstreamnode = NULL;
         }
      }
      else if (0) //DEBUG
	std::cout<<"no upstream node for node "<<cn->getID()<<std::endl;
      return upstreamnode;
   }
   else return NULL;
}


/****************************************************************\
 **
**
**		Parameters:
**		Called by:
**		Created: 1/98 SL
**
\***************************************************************/


