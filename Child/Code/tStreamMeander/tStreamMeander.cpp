/**************************************************************************/
/**
**  @file tStreamMeander.cpp
**  @brief Functions for class tStreamMeander.
**
**  $Id: tStreamMeander.cpp,v 1.98 2003-11-03 12:42:08 childcvs Exp $
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
  critflow(0), optdiamvar(0), optrainvar(0),
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
  //if( netPtr != 0 ) netPtr = new tStreamNet( gRef );
  netPtr = &netRef;
  assert( netPtr != 0 );
  meshPtr = &mRef;
  assert( meshPtr != 0 );
  infilePtr = &infile;
  assert( infilePtr != 0 );
  critflow = infilePtr->ReadItem( critflow, "CRITICAL_FLOW" );
  optdiamvar = infilePtr->ReadItem( optdiamvar, "OPT_VAR_SIZE" );
  if( !optdiamvar )
    {
      meddiam = infilePtr->ReadItem( meddiam, "MEDIAN_DIAMETER" );
      assert( meddiam > 0 );
    }
  optrainvar = infilePtr->ReadItem( optrainvar, "OPTVAR" );
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
  tMeshListIter< tLNode > nodIter( meshPtr->getNodeList() );

  if (0) //DEBUG
    cout << "FindMeander...()";

  for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
    {
      //nmg
      if (0) //DEBUG
	cout<<"FM cn "<<cn->getID()<<" z = "<<cn->getZ()<<endl;
      cn->setReachMember( false );
      if( cn->getQ() >= critflow )
	{
	  cn->setMeanderStatus( kMeanderNode );
	  cn->setNew2DCoords( cn->getX(), cn->getY() );
	  if (0) //DEBUG
	    cout << "FM node " << cn->getID() << " meanders discharge " << cn->getQ()<<" crit Q "<<critflow<<endl;
	}
      else
	{
	  cn->setMeanderStatus( kNonMeanderNode );
	  if (0) {//DEBUG
	    if( cn->getQ() >= critflow )
	      cout << "node " << cn->getID() << " has Q " << cn->getQ()
		   << " but is flooded" << endl;
	  }
	}
    }
  if (0) //DEBUG
    cout << "done\n";
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
      cout << "InterpChannel()\n";
  }
  const tArray< double > zeroArr(4);
  double slope=0.;
  bool change = false; //haven't added any nodes yet
  //loop through reaches:
  if (0) {//DEBUG
    if( timetrack >= kBugTime ) cout << "loop through reaches" << endl;
  }

  int i;
  tPtrList< tLNode > *creach;
  for( creach = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
       creach = rlIter.NextP(), i++ )
    {
      if (0) //DEBUG
	cout<<"reach " << i << endl;
      tPtrListIter< tLNode > rnIter( *creach );
      const int num = nrnodes[i];
      //loop through reach nodes:
      if (0) {//DEBUG
	if( timetrack >= kBugTime ) cout << "loop through reach nodes" << endl;
      }
      int j;
      tLNode *crn;
      for( crn = rnIter.FirstP(), j=0; j<num;
           crn = rnIter.NextP(), j++ )
	{
	  if (0) {//DEBUG
	    if( timetrack >= kBugTime ) cout << "node " << crn->getID() << endl;
	  }
	  // GT changed from hydr to chan width, 6/99:
	  const double curwidth = crn->getChanWidth();
	  //curwidth = crn->getHydrWidth();
	  if (0){ //DEBUG
	    if( timetrack >= kBugTime ) cout << "found width" << endl;
	  }
	  tLNode *nPtr = crn->getDownstrmNbr();
	  if (0){ //DEBUG
	    if( timetrack >= kBugTime ) cout << "found dnstrm nbr" << endl;
	  }
	  const double curseglen = crn->getFlowEdg()->getLength();
	  if (0){ //DEBUG
	    if( timetrack >= kBugTime ) cout << "found flowedge length" << endl;
	  }
	  assert( curseglen > 0 );
	  //ran into an infinite loop in tLNode::GetSlope; now, if it runs into
	  //an infinite loop, it returns a negative number (-1) as an alarm
	  //flag instead of generating a fatal error; in such case, bail out of
	  //these loops through the nodes and call UpdateNet():
	  slope = crn->calcSlope();
	  if( slope < 0.0 )
	    {
	      change = true;
	      break;
	    }
	  if (0){ //DEBUG
	    if( timetrack >= kBugTime ) cout << "found slope" << endl;
	  }
	  const double defseglen = dscrtwids * curwidth;
	  const double maxseglen = 2.0 * defseglen;
	  //if current flowedg length is greater than
	  //a given proportion of hydraulic width
	  if( curseglen > maxseglen )
	    {
	      change = true;//flag so we know that we need to update the network
	      tLNode nn = *crn; //added nodes are copies of the upstream node except xyz.
	      nn.setXYZD( zeroArr ); //and xyzd ('old' coords)

	      const double x0 = crn->getX();
	      const double y0 = crn->getY();
	      const double z0 = crn->getZ();
	      const double x1 = nPtr->getX();
	      const double y1 = nPtr->getY();
	      const double y = y1 - y0;
	      const double x = x1 - x0;
	      const double phi = atan2( y, x );
	      const double bigseglen = 3.0 * defseglen;
	      //if we can get the desired spacing by adding a single node
	      //use linear interpolation with noise added
	      //(for Delaunay mesh, avoid precisely co-linear points):
	      if( curseglen <= bigseglen )
		{
		  const double val1 = rand->ran3() - 0.5;
		  const double x = x0 + curseglen * cos(phi) /
		    2.0 + 0.01 * val1 * defseglen;//pow(defseglen, 2.0);
		  const double val2 = rand->ran3() - 0.5;
		  const double y = y0 + curseglen * sin(phi) /
		    2.0 + 0.01 * val2 * defseglen;//pow(defseglen, 2.0);
		  const double z = z0 - curseglen / 2.0 * slope;
		  nn.set3DCoords( x, y, z );
		  if (0){ //DEBUG
		    if( timetrack >= kBugTime ) cout << "add a node" << endl;
		  }
		  // set flow edge temporarily to zero, so that it is flippable
		  crn->setFlowEdgToZero();
		  // Add new node
		  tLNode * newnodeP = meshPtr->AddNode( nn, kUpdateMesh, time );
		  if (0) //DEBUG
		    cout<<"IC pt added at " << x << "," << y << endl;

		  // Reconnect flow edges
		  newnodeP->flowTo( nPtr );
		  crn->flowTo( newnodeP );
		  // Paranoia
		  assert( newnodeP->getFlowEdg()->getDestinationPtr() == nPtr );
		  assert( crn->getFlowEdg()->getDestinationPtr() == newnodeP );
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
		  tLNode *prevNode = crn;
		  for(int i=1; i<npts; i++ )
		    {
		      xp[i] = static_cast<double>(i) * defseglen;
		      const double val = rand->ran3() - 0.5;
		      yp[i] = yp[i-1] + 0.01 * val * exp(-yp[i-1] * val)
			* defseglen;//pow(defseglen, 2.0);
		      zp[i] = xp[i] * slope;
		      const double dd = sqrt(xp[i] * xp[i] + yp[i] * yp[i]);
		      const double da = atan2(yp[i], xp[i]);
		      const double x = dd * cos(da + phi) + x0;
		      const double y = dd * sin(da + phi) + y0;
		      const double z = z0 - zp[i];
		      if (0) //DEBUG
			cout<<"InterpChannel: call AddNode"<<endl;
		      nn.set3DCoords( x, y, z );
		      if (0){ //DEBUG
			if( timetrack >= kBugTime ) cout << "add a node" << endl;
		      }
		      // set flow edge temporarily to zero, so that it is flippable
		      prevNode->setFlowEdgToZero();
		      tLNode * newnodeP = meshPtr->AddNode( nn, kUpdateMesh, time );
		      if (0) //DEBUG
			cout<<"IC pt added at " << x << "," << y << endl;

		      // Reconnect flow edges
		      newnodeP->flowTo( nPtr );
		      prevNode->flowTo( newnodeP );
		      // Update
		      prevNode = newnodeP;
		    }
		}
	    }
	}
      if( slope < 0.0 ) break; //'change' has already been set to 1
    }
  if( change )
    {
      if (0){ //DEBUG
	if( timetrack >= kBugTime ) cout << "update network" << endl;
      }
      //gridPtr->UpdateMesh();
      //netPtr->InitFlowDirs();
      //NOTE****!!! the zero param below should be replaced with current time,
      // which needs to be passed to Migrate, etc....TODO
      netPtr->UpdateNet( time );
      if (0){ //DEBUG
	if( timetrack >= kBugTime ) cout << "added nodes(s), InterpChannel finished"
					 << endl;
      }
      return 1;
    }
  else
    {
      if (0){ //DEBUG
	if( timetrack >= kBugTime ) cout << "InterpChannel finished" << endl;
      }
      return 0;
    }

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
    cout<<"tStreamMeander::MakeReaches...";

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
    cout<<"done MakeReaches"<<endl;
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
  tLNode *cn, *lnPtr, *frn, *lrn, *nxn, *ln;
  tEdge *ce;
  tMeshListIter< tLNode > nodIter( meshPtr->getNodeList() );
  tPtrListIter< tLNode > rnIter;
  tPtrList< tLNode > rnodList, *plPtr, listtodelete;
  tListNode< tPtrList< tLNode > > *tempnode;

  if (0) //DEBUG
    cout<<"tStreamMeander::FindReaches()"<<endl;

  if( !(reachList.isEmpty()) ) reachList.Flush();
  //loop through active nodes
  for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
    {
      //if node meanders
      if( cn->Meanders() )
	{
	  int nmndrnbrs = 0;

	  if (0) //DEBUG
	    cout<<"FR node "<<cn->getID()<<" meanders"<<endl;
	  tSpkIter spokIter( cn );
	  //loop through spokes to find upstream meandering nbr
	  for( ce = spokIter.FirstP(); !(spokIter.AtEnd());
	       ce = spokIter.NextP() )
	    {
	      lnPtr = static_cast<tLNode *>( ce->getDestinationPtrNC() );
	      //lnPtr points to the downstream neighbor of the current edge
	      if( lnPtr->getBoundaryFlag() == kNonBoundary
		  && lnPtr->getDownstrmNbr() == cn
		  && lnPtr->Meanders() )
		  {
		    //If you enter here, the cn is downstream of one
		    //of the nodes that one of the spokes on the spoke
		    //list is pointing to.
		    nmndrnbrs++;
		    break;
		  }
	    }
	  //if no upstream meandering nbr, then it's a reach "head";
	  //put it at the beginning of an otherwise empty reach
	  //and add the reach to the list of reaches:
	  if (nmndrnbrs==0)
	    {
	      rnodList.Flush();
	      rnodList.insertAtFront( cn );
	      if (0){ //DEBUG
		if(cn->calcSlope()<=0) cout<<"bad node "<<cn->getID()<<" with slope "<<cn->calcSlope()<<" added to reachlist"<<endl;
	      }
	      reachList.insertAtBack( rnodList );
	      cn->setReachMember( true );
	    }
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
    cout << "No. reaches: " << reachList.getSize() << endl;
  //loop through reaches
  //rlIter is the iterator for reachlist with is a list of ptrs to node lists
  for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
       plPtr = rlIter.NextP(), i++ )
    {
      assert( reachList.getSize() > 0 );
      if (0) //DEBUG
	cout << " on reach " << i << endl;
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
	    if(cn->calcSlope()<=0) cout<<"bad node "<<cn->getID()<<" with slope "<<cn->calcSlope()<<" added to reachlist"<<endl;
	  }
	  nrnodes[i]++;
	  reachlen[i] += cn->getFlowEdg()->getLength();
	  cn->setReachMember( true );
	  plPtr->insertAtBack( cn );
	  cn = cn->getDownstrmNbr();
	}

      //make sure reach has more than 4 members:
      if (0) //DEBUG
	cout << "reach " << i << " length " << plPtr->getSize() << endl;

      if( plPtr->getSize() > 4 )
	{
	  //make sure reach has "positive" slope:
	  lrn = rnIter.LastP();
	  frn = rnIter.FirstP();
	  rdrop = frn->getZ() - lrn->getZ();
	  if (0) //DEBUG
	    cout << "reach drop: " << rdrop << endl;
	  if( rdrop <= 0 )
	    {
	      //remove reach if it has non-positive slope:
	      if (0) //DEBUG
		cout << "remove reach for non-positive slope" << endl;
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
	    cout << "remove reach for being too small" << endl;
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
    cout << "Final no. reaches: " << reachList.getSize() << endl;
    for( /*plPtr =*/ rlIter.FirstP(), i=0; !(rlIter.AtEnd());
	 /*plPtr =*/ rlIter.NextP(), i++ )
      {
	cout << "reach " << i << " length " << nrnodes[i] << endl;
      }
  }
  //construct a tail for the reach
  //which is some number of hydr. widths long:
  for( plPtr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
       plPtr = rlIter.NextP(), i++ )
    {
      if (0) //DEBUG
	cout << " on reach " << i << endl;
      assert( i<reachList.getSize() );
      rnIter.Reset( *plPtr );
      cn = rnIter.LastP();
      // GT changed this to chan width, 6/99:
      curwidth = cn->getChanWidth();
      //curwidth = cn->getHydrWidth();
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
		  cout << "reach with " << plPtr->getSize() << " nodes: " << endl;
		  for( ln = rnIter.FirstP(); !(rnIter.AtEnd()); ln = rnIter.NextP() )
		    {
		      cout << ln->getID() << "; ";
		    }
		  cout << endl;
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
	cout << "end FindReaches, node " << cn->getID() << endl;
      }
  }
  if (0) //DEBUG
    cout << "done FindReaches" << endl;
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
**
\***************************************************************/
void tStreamMeander::CalcMigration( double &time, double const &duration,
                                    double &cummvmt )
{
  if (0) //DEBUG
    cout<<"tStreamMeander::CalcMigration()...";

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
	    //set old x,y if necessary
	    tArray< double > oldpos = curnode->getXYZD();
	    if( oldpos[0] == 0.0 && oldpos[1] == 0.0 )
	      {
		oldpos[0] = curnode->getX();
		oldpos[1] = curnode->getY();
		curnode->setXYZD( oldpos );
		if (0) //DEBUG
		  cout << "set old x,y to current coords for node "
		       << curnode->getID() << endl;
	      }
	    if (0){ //DEBUG
	      const tArray< double > newxy = curnode->getNew2DCoords();
	      cout << "init. new coords to " << newxy[0] << " " << newxy[1] << endl;
	    }
	  }
      }
      if (0) //DEBUG
	cout << "reach " << i << " length " << nrnodes[i] << endl;

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
	    //xs += fedg->getLength(); BUG?? GT 3/16/99 seems to dup line below!
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
	    //curnode->GetErodibility( rerody, lerody, pr->kf);
	    lerodya[j] = bankerody[0];
	    rerodya[j] = bankerody[1];

	    // Set slope, width, depth, grainsize, and roughness arrays
	    slopea[j] = curnode->getHydrSlope();
	    widtha[j] = curnode->getHydrWidth();
	    deptha[j] = curnode->getHydrDepth();
	    if (0) //DEBUG
	      cout << "width, depth " << widtha[j] << " " << deptha[j] << endl;
	    /*diama[j] = curnode->diam;*/
	    diama[j] = ( optdiamvar ) ? curnode->getDiam() : meddiam;
	    lambdaa[j] = curnode->getBankRough();
	  }
      }

      // Now we pass all this information to meander.f
      if (0) //DEBUG
	cout << "stations, stnserod: " << stations <<" "<< nttlnodes
	     << endl;
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
      double dbg=0.;
      static double cumdbg=0.0; //debug
      {
	tLNode *curnode;
	int j;
	for( curnode = rnIter.FirstP(), j=0; !(rnIter.AtEnd());
	     curnode = rnIter.NextP(), j++ )
	  {
	    dbg+=deltaxa[j];
	    cumdbg+=deltaxa[j];
	    curnode->addLatDisplace( deltaxa[j], deltaya[j] );
	  }
	if (0) //DEBUG
	  cout << "MEAN X CHG: " << dbg/static_cast<double>(j)
	       << "  CUM: " << cumdbg << endl;
	//assert( cumdbg>-1e6 );
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
	    const double
	      rz = curnode->getZ() + deptha[j] - rdeptha[j],
	      lz = curnode->getZ() + deptha[j] - ldeptha[j];
	    //rz = curnode->getZ() + deptha[j];
	    //lz = rz + deptha[j];
	    dbg2 = dbg2 + (rdeptha[j]-ldeptha[j]);
	    curnode->setZOld( rz, lz );
	  }
	if (0) //DEBUG
	  cout << "MEAN rldepth " << dbg2 << endl;
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
  //if( maxfrac < allowfrac )
  //    dtm = ( maxfrac > 0 ) ? allowfrac / maxfrac : 1.0;
  //else dtm = allowfrac / maxfrac;

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
	  //newxy = curnode->getNew2DCoords();
	  if (0) //DEBUG
	    cout << "new coords set to " << newxy[0] << " " << newxy[1] << endl;
	  //XXXdo 'uplift' on old z-values:
	  //uplift creates spikes if points are added after a long time
	  /*oldpos = curnode->getXYZD();
	    if( oldpos[3] != 0.0 )
	    {
            oldpos[2] += curnode->getUplift() * dtm;
            curnode->setXYZD( oldpos );
	    }*/
	}
    }
  time += dtm;
  cummvmt += maxfrac * dtm;

  if (0) //DEBUG
    cout<<"done CalcMigration\n";

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
    cout<<"Migrate time=" << ctime <<endl;
  const double duration = netPtr->getStormPtrNC()->getStormDuration()
    + ctime;
  double cummvmt = 0.0;
  //timeadjust = 86400. * pr->days;

  //NOTE: ctime and duration involve close subtraction of increasingly
  // large #'s -- shouldn't we just use the "true" duration? TODO

  while( ctime < duration )
    {
      MakeReaches( ctime ); //updates net, makes reachList
      if( !(reachList.isEmpty()) )
	{
	  if (0) //DEBUG
	    cout<<"in loop "<<ctime<<" duration is "<<duration<<endl;
	  CalcMigration( ctime, duration, cummvmt ); //incr time; uses reachList
	  MakeChanBorder( ); //uses reachList
	  CheckBndyTooClose();  //uses tMesh::nodeList
	  CheckBanksTooClose(); //uses reachList
	  CheckFlowedgCross(); //uses tMesh::nodeList
	  CheckBrokenFlowedg(); //uses tMesh::nodeList
	  meshPtr->MoveNodes( ctime ); //uses tMesh::nodeList
	  AddChanBorder( ctime ); //uses tMesh::nodeList
	}
      else ctime=duration; // If no reaches, end here (GT added 3/12/99)
      //MakeReaches(); had called from main routine and here
    }
  if (0) //DEBUG
    cout<<"end migrate"<<endl;
}


/****************************************************************************\
**
**  tStreamMeander::MakeChanBorder()
**
**  This function sets the "old coordinates" of the
**  migrating node; i.e., the coordinates at which a new
**  node will be "dropped" by a meander node. Checks present
**  coords against the coords set in CalcMigration. If the
**  present coords are more than half a hydraulic width
**  away from the old coords, then the old x, y, and z
**  values are set to correspond to the present location
**  of the bank away from which the node is moving and
**  the bed elevation at that bank, respectively. Also,
**  "d" is set to indicate which side of the channel that
**  bank is on, 1 for the left side, -1 for the right side.
**  That way, now and when it's time for the node to be dropped,
**  we can make sure it's still on the same side of the
**  channel as when the coords were set. If the channel
**  side changes, reset z and d to zero.
**
**   All of this is done after CalcMigration has been called
**   but before tMesh::MoveNodes is called.
**
**              Parameters:
**              Called by:      Migrate
**              Created:        8/18/97 SL
**              Updated: 1/98 SL
**
\*****************************************************************************/
void tStreamMeander::MakeChanBorder( )
{
  if (0) //DEBUG
    cout << "MakeChanBorder()" << endl;
  int i, j, num, pccw;
  double x0, y0, x1, y1, z, delx, dely, phi, width, xdisp, ydisp;
  tPtrList< tLNode > *cr;
  tPtrListIter< tLNode > rnIter;
  tLNode *cn, *cnbr;
  tArray< double > cnpos, dsnpos, oldpos, rl;

  for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
       cr = rlIter.NextP(), i++ )
    {
      rnIter.Reset( *cr );
      num = nrnodes[i];
      for( cn = rnIter.FirstP(), j=0; j<num;
           cn = rnIter.NextP(), j++ )
	{
	  // GT changed from hydr to chan width, 6/99
	  width = cn->getChanWidth();
	  //width = cn->getHydrWidth();
	  oldpos = cn->getXYZD();
	  if( oldpos[3] == 0.0 ) //if channel side is still unset,
	    {
	      if( cn->DistFromOldXY() >= 0.5 * width ) //if moved more than 1/2-width
		{
		  if (0) //DEBUG
		    cout << "node " << cn->getID()
			 << " >= width/2 from old coords" << endl;
		  //find coordinates of banks, i.e., coords at n = + and -width
		  cnpos = cn->get2DCoords();
		  cnbr = cn->getDownstrmNbr();
		  dsnpos = cnbr->get2DCoords();
		  x0 = cn->getX();
		  y0 = cn->getY();
		  x1 = cnbr->getX();
		  y1 = cnbr->getY();
		  delx = x1 - x0;
		  dely = y1 - y0;
		  phi = atan2( dely, delx );
		  xdisp = 0.5 * width * sin(phi);
		  ydisp = 0.5 * width * cos(phi);
		  rl = cn->getZOld();
		  z = cn->getZ();
		  //find which side of the channel the old coords are at:
		  if( PointsCCW( cnpos, dsnpos, oldpos ) )
		    {
		      oldpos[0] = x0 - xdisp;
		      oldpos[1] = y0 + ydisp;
		      oldpos[2] = ( rl[1] > z ) ? rl[1] : z;
		      oldpos[3] = 1.0;
		      if (0) //DEBUG
			cout << "node " << cn->getID()
			     << " old pos set, on left side of channel" << endl;
		    }
		  else
		    {
		      oldpos[0] = x0 + xdisp;
		      oldpos[1] = y0 - ydisp;
		      oldpos[2] = ( rl[0] > z ) ? rl[0] : z;
		      oldpos[3] = -1.0;
		      if (0) //DEBUG
			cout << "node " << cn->getID()
			     << " old pos set, on right side of channel" << endl;
		    }
		  //if( oldpos[2]>2.0 ) cout <<"**** OLDPOS z = " << oldpos[2] << endl;
		  cn->setXYZD( oldpos );  //coords found, update node data member
		}
	    }
	  else
	    {
	      //make sure old coords are on same side of channel
	      //they were last time:
	      cnpos = cn->get2DCoords();
	      cnbr = cn->getDownstrmNbr();
	      dsnpos = cnbr->get2DCoords();
	      pccw = PointsCCW( cnpos, dsnpos, oldpos );
	      //if not, erase old z (i.e., set to zero):
	      if( !( ( pccw && oldpos[3] == 1.0 ) ||
		     ( !pccw && oldpos[3] == -1.0 ) ) )
		{
		  oldpos[2] = 0.0;
		  oldpos[3] = 0.0;
		  cn->setXYZD( oldpos );
		  if (0) //DEBUG
		    cout << "node " << cn->getID()
			 << " switched sides of channel: reinitialize z coord."
			 << endl;
		}
	    }
	}
    }
  if (0){ //DEBUG
    tMeshListIter< tLNode > nIter( meshPtr->getNodeList() );
    for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
      {
	cout << "node " << cn->getID() << endl;
	if( cn->Meanders() )
	  {
	    cout << "   meanders" << endl;
	  }
      }
  }
}


/****************************************************************************\
**
**	tStreamMeander::AddChanBorder
**
**  For meandering nodes with placement coords set, check whether
**  a new node should be dropped. First, check the distance from
**  the old coords vs. the distance specified in the input file.
**  If far enough, make sure the old coords are
**  not presently in the channel and that they are on the same
**  side of the channel as when they were set. If either of
**  these conditions are not met, reset old coords to zero.
**  If the conditions are met, add a new node.
**
**  Parameters:     uses node->hydrwidth and, indirectly,
**                  parameters which determine hydraulic geometry
**  Called by:      Migrate
**  Created:        8/18/97 SL
**  Updated:        1/98 SL; 2/98 SL
**   - 6/99 GT: fixed bug in which channode was not initialized before
**     being passed to InChannel. The effect of this was probably to
**     allow new nodes to land inside channels.
**
\***************************************************************************/
void tStreamMeander::AddChanBorder(double time)
{
  if (0) //DEBUG
    cout << "AddChanBorder()" << endl;
  int i, inchan, pccw, sameside;
  double lvdist, width;
  tArray< double > oldpos, zeroArr(4);
  tTriangle *ct;
  tLNode *cn, *tn, *channodePtr, channode;
  tMeshListIter< tLNode > nIter( meshPtr->getNodeList() ),
    tI( meshPtr->getNodeList() );

  //go through active nodes:
  for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
    {
      //select for meandering nodes:
      if( cn->Meanders() )
	{
	  oldpos = cn->getXYZD();
	  //select for nodes with old coords set:
	  if( oldpos[3] != 0.0 )
	    {
	      // GT changed hydr width to chan width, 6/99:
	      width = cn->getChanWidth();
	      //width = cn->getHydrWidth();
	      lvdist = leavefrac * width;
	      //select for nodes far enough away from the old coords:
	      if( cn->DistFromOldXY() > lvdist )
		{
		  if (0) //DEBUG
		    cout << "node " << cn->getID()
			 << " ready to drop new node" << endl;
		  //just make sure new node will be
		  //(a) not in a channel and
		  //(b) on the same side of the channel:
		  if( (ct = meshPtr->LocateTriangle( oldpos[0], oldpos[1] ))
		      != NULL)
		    {
		      // WHY COMMENTED OUT??
		      //channodePtr = cn;
		      //channode = *channodePtr;
		      //***NG: HERE IS WHERE YOU CAN FIND A DEPOSIT THICKNESS
		      //TO ADD TO THE NEW NODE***
		      // GT uncommented this line to fix bug, 6/99:
		      channode.set3DCoords( oldpos[0], oldpos[1], oldpos[2] );
		      //channode.setXYZD( zeroArr ); //initialize xyzd
		      //channode.SetMeanderStatus(0); //zero meander
		      inchan = 0;
		      for( i=0; i<3; i++ )
			{
			  tn = static_cast<tLNode *>( ct->pPtr(i) );
			  if( tn->Meanders() )
			    {
			      // ! This passes uninitialized channode ptr ! TODO
			      if( (inchan = InChannel( tn, &channode )) != 0)
				{
				  if (0) //DEBUG
				    cout << "old coord's in channel" << endl;
				  break;
				}
			    }
			}
		      if( !inchan )
			{
			  sameside = 1;
			  for( i=0; i<3; i++ )
			    {
			      tn = static_cast<tLNode *>( ct->pPtr(i) );
			      if( tn->Meanders() )
				{
				  const tArray< double >
				    xy = tn->get2DCoords(),
				    xyd = tn->getDownstrmNbr()->get2DCoords();
				  pccw = PointsCCW( xy, xyd, oldpos );
				  if( !( ( pccw && oldpos[3] == 1.0 ) ||
					 ( !pccw && oldpos[3] == -1.0 ) ) )
				    {
				      sameside = 0;
				      if (0) //DEBUG
					cout << "old coord's switched sides" << endl;
				      break;
				    }
				}
			    }
			  if( sameside )
			    {
			      if (0){ //DEBUG
				cout << "node " << cn->getID()
				     << "'s old coords pass: add new node" << endl;
				cout << "** ELEV is " << oldpos[2] << endl;
			      }
			      tArray< double > xyz(3);
			      for( i=0; i<3; i++ ) xyz[i] = oldpos[i];
			      // Make sure the banknode is not lower than the node it
			      // originates from, bug fix 8/2003 QC. Causes ponds if the
			      // meander path is redirected over the newly added banknode
			      // in FlowDir
			      if( xyz[2] < cn->getZ())
				xyz[2] = cn->getZ();
			      channodePtr = meshPtr->AddNodeAt( xyz, time );
			      channodePtr->setRock( cn->getRock() );
			      //XchannodePtr->setSurf( cn->getSurf() );
			      channodePtr->setVegCover( cn );  // GT 1/2000
			      channodePtr->setReg( cn->getReg() );
			      //TODO: NG Need to take care of deposit depth here
			      //I was thinking to leave a deposit of depth
			      //xyz[2]-cn->getZ() if this depth is positive
			      //The texture of this deposit would be
			      //the surface texture of cn.  Use erodep.
			      //if(xyz[2]-cn->getZ()>0){
			      //}
			      //gridPtr->AddNode( channode );
			      for( i=0; i<4; i++ ) oldpos[i] = 0.0;
			      cn->setXYZD( oldpos );
			      if (0) //DEBUG
				cout << "node " << cn->getID()
				     << ": reinitialize old coords" << endl;
			    }
			  else
			    {

			      for( i=0; i<4; i++ ) oldpos[i] = 0.0;
			      cn->setXYZD( oldpos );
			      if (0) //DEBUG
				cout << "node " << cn->getID()
				     << ": reinitialize old coords" << endl;
			    }
			}
		      else
			{
			  if (0) //DEBUG
			    cout << "if old coord's not in channel, no vtx was meandering"
				 << endl;
			  if( i == 3 ) //no vtx was meandering
			    {
			      for( i=0; i<4; i++ ) oldpos[i] = 0.0;
			      cn->setXYZD( oldpos );
			      if (0) //DEBUG
				cout << "node " << cn->getID()
				     << ": reinitialize old coords" << endl;
			    }
			}
		    }
		  else
		    {
		      if (0) //DEBUG
			cout << "old coords not in any triangle" << endl;
		      for( i=0; i<4; i++ ) oldpos[i] = 0.0;
		      cn->setXYZD( oldpos );
		      if (0) //DEBUG
			cout << "node " << cn->getID()
			     << ": reinitialize old coords" << endl;
		    }
		}
	    }
	}
    }

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
    cout << "FBE\n";

  tArray< double > lrerody(2);
  if( nPtr->getBoundaryFlag() != kNonBoundary ) return lrerody;

  tSpkIter sI(nPtr);
  tArray< double > spD( sI.getNumSpokes() * 2 ),
    spR( sI.getNumSpokes() * 2 );
  int i, j, n;
  tLNode *nNPtr, *cn, *node1, *node2;
  tEdge *ce, *fe, *ne;
  double a, b, c, d, s1, s2, D, d1, d2, E1, E2, dz1, dz2, H;

  nNPtr = nPtr->getDownstrmNbr();
  const tArray< double > xyz1 = nPtr->get3DCoords();
  const tArray< double > xy2 = nNPtr->get2DCoords();
  tArray< double > dxy(2);
  dxy[0] = xy2[0] - xyz1[0];
  dxy[1] = xy2[1] - xyz1[1];
  a = dxy[0];
  b = dxy[1];
  c = -dxy[1] * xyz1[1] - dxy[0] * xyz1[0];
  n = sI.getNumSpokes();
  //find distance and remainders of pts wrt line perpendicular to downstream direction
  fe = nPtr->getFlowEdg();
  ce = fe;
  i = 0;
  do
    {
      cn = static_cast<tLNode *>( ce->getDestinationPtrNC() );
      const tArray< double > xy = cn->get2DCoords();
      d = a * xy[0] + b * xy[1] + c;
      spD[i] = spD[n + i] = d;
      //if d=0 then point is on line
      spR[i] = spR[n + i] = DistanceToLine( xy[0], xy[1], a, b, c );
      ce = ce->getCCWEdg();
      i++;
    } while( ce != fe );
  H = nPtr->getHydrDepth();
  //find point pairs:
  i=0;
  j=0;
  ce = fe;
  do
    {
      //find signs of 'this' and 'next' remainders
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
	  d1 = spR[i];
	  d2 = spR[i + 1];
	  D = d1 + d2;
	  //find erodibilities of nbr nodes wrt nPtr
	  ne = ce->getCCWEdg();
	  node1 = static_cast<tLNode *>( ce->getDestinationPtrNC() );
	  node2 = static_cast<tLNode *>( ne->getDestinationPtrNC() );
	  //find elev. diff's:
	  dz1 = node1->getZ() - xyz1[2];
	  dz2 = node2->getZ() - xyz1[2];

	  //find whether bedrock or alluvial bank:
	  // (modified by GT 3/99: getAlluvThickness is obsolete. Layer
	  // info should be used. For now assume constant erodibility. TODO)
	  /*if( dz1 > node1->getAlluvThickness() ) E1 = rockerod;//node1->getBedErody();
	    else E1 = vegerod;//node1->getVegErody();
	    if( dz2 > node2->getAlluvThickness() ) E2 = rockerod;//node2->getBedErody();
	    else E2 = vegerod;//node2->getVegErody();*/
	  E1 = E2 = rockerod; // added 3/99

	  if (0) //DEBUG
	    cout << "E1 " << E1 << "  E2 " << E2 << endl;
	  //find height dependence:
	  //if elev diff > hydraulic depth, ratio of depth to bank height;
	  //o.w., keep nominal erody:
	  if (0) //DEBUG
	    cout << "FBE 4" << H << " " << dz1 << endl;
	  if( dz1 > H ) E1 *= (Pdz * H / dz1 + (1 - Pdz));
	  if (0) //DEBUG
	    cout << "FBE 5" << endl;
	  if( dz2 > H ) E2 *= (Pdz * H / dz2 + (1 - Pdz));
	  //now we've found erod'ies at ea. node, find weighted avg:
	  assert( D > 0.0 );
	  if (0) //DEBUG
	    cout << "FBE 6" << endl;
	  lrerody[j] = (E1 * d2 + E2 * d1) / D;
	  j++;
	  if (0) //DEBUG
	    cout << "FBE 7" << endl;
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
  tMeshListIter< tLNode > nI( meshPtr->getNodeList() );
  tLNode *cn, *mn;
  tNode *nn, *bn0(0), *bn1(0);
  tEdge *ce;
  int n;
  double width, mindist, d0, d1, d2, d3, xp, yp;

  if (0) //DEBUG
    cout << "CBTC\n";

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
		  //else cout << 'Warning: >2 bndy nbrs found in CheckBndyTooClose\n';
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
    cout << "CheckBanksTooClose()..." << endl;
  int i, j, num, onlist;
  tPtrList< tLNode > delPtrList;
  tPtrListIter< tLNode > dIter( delPtrList );
  tLNode * cn, *pointtodelete, *dn, *sn;
  tEdge *ce;
  tPtrList< tLNode > *cr;          // ptr to current reach
  tPtrListIter< tLNode > rnIter;   // iterator for nodes on reach

  // For each reach
  for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
       cr = rlIter.NextP(), i++ )
    {
      // For each node on reach
      rnIter.Reset( *cr );
      num = nrnodes[i];
      for( cn = rnIter.FirstP(), j=0; j<num;
           cn = rnIter.NextP(), j++ )
	{
	  // Check neighboring nodes
	  tSpkIter spokIter( cn );
	  for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
	       ce = spokIter.NextP() )
	    {
	      sn = static_cast<tLNode *>( ce->getDestinationPtrNC() );
              //check for proximity to channel:
	      if( !(sn->Meanders()) && InChannel( cn, sn ) )
		{
		  // If node isn't a boundary and isn't already on the
		  // deletion list, put it on the deletion list now
		  if (0) //DEBUG
		    cout<<"too close: cn, cn->hydrwidth "<<cn->getID()<<" "
			<<cn->getHydrWidth()<<endl;
		  pointtodelete = static_cast<tLNode *>( ce->getDestinationPtrNC() );
		  if( pointtodelete->getBoundaryFlag() == kNonBoundary )
		    {
		      if ( pointtodelete->getDrArea() < cn->getDrArea() )
			{
			  onlist = 0;
			  for( dn = dIter.FirstP(); !(dIter.AtEnd());
			       dn = dIter.NextP() )
			    if( pointtodelete == dn ) onlist = 1;
			  if( !onlist )
			    {
			      if (0) //DEBUG
				cout << "add to delete list: "
				     << pointtodelete->getID() << endl;
			      delPtrList.insertAtBack( pointtodelete );
			    }
			}
		    }
		  else cn->RevertToOldCoords();
		}
	    }
	}
    }

  // Having found all the nodes that have been "swept away" by the
  // channel and placed them on the delPtrList, we now delete them
  for( dn = dIter.FirstP(); !(dIter.AtEnd()); dn = dIter.FirstP() )
    {
      if (0) //DEBUG
	cout << "CBTC: delete node " << dn->getID() << endl;
      meshPtr->DeleteNode( dn );
      /*cn =*/ delPtrList.removeFromFront();
    }

  if (0) //DEBUG
    cout << "finished" << endl;
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
**
\*****************************************************************************/
void tStreamMeander::CheckFlowedgCross()
{
  if (0) //DEBUG
    cout << "CheckFlowedgCross()..." << endl;
  int i, j;
  int ft;
  bool crossed;
  tLNode *pointtodelete, *nod(0), *cn, *dscn(0);
  tEdge * fedg, *ce;
  tTriangle * ct, *nt;
  tListIter< tTriangle > triIter( meshPtr->getTriList() );
  //delete node crossed by flowedg:
  //if a new triangle is !CCW and two vtcs. are connected by a flowedg AND
  //  a spoke of the third vtx. intersects the flowedg OR
  //  more than one nbr. tri. is also !CCW
  //then delete third vtx. because it has been "crossed" by the flowedg
  do
    {
      crossed = false;
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
	{
	  if( !NewTriCCW( ct ) )
	    {
	      ft = 0;
	      pointtodelete = 0;
	      for( i=0; i<3; i++ )
		{
		  cn = static_cast<tLNode *>( ct->pPtr(i) );
		  if( cn->Meanders() )
		    {
		      dscn = cn->getDownstrmNbr();
		      for( j=1; j<3; j++ )
			{
			  nod = static_cast<tLNode *>( ct->pPtr( (i+j)%3 ) );
			  //if another vtx is cn's downstream nbr
			  if( nod == dscn )
			    {
			      ft = 1;
			      //set 'nod' to third node
			      if( j == 1 ) nod = static_cast<tLNode *>( ct->pPtr( (i+2)%3 ) );
			      else nod = static_cast<tLNode *>( ct->pPtr( (i+1)%3 ) );
			    }
			  if( ft ) break;
			}
		    }
		  if( ft ) break;
		}
              //if 'flow triangle', i.e., meandering node and dnstrm nbr
	      if( ft )
		{
		  //check nbr tri's for new pos'n CCW
		  for( i=0, j=0; i<3; i++ )
		    {
		      nt = ct->tPtr(i);
		      if( nt != 0 )
			if( !NewTriCCW( nt ) )
                          j++;
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
		      for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
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
		      crossed = true;
		      //if the node to delete is a meandering node
		      //with greater flow, delete the node we started with
		      if( pointtodelete->Meanders() &&
			  cn->getDrArea() < pointtodelete->getDrArea() )
			{
			  pointtodelete = cn;
			}
		      //don't delete boundary node
		      else if( pointtodelete->getBoundaryFlag() )
			{
			  pointtodelete = 0;
			}
		      if( pointtodelete != 0 )
			{
			  meshPtr->DeleteNode( pointtodelete );
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
			  //assert( gridPtr->LocateTriangle( xy0[0], xy0[1] ) == 0 ||
			  //        gridPtr->LocateTriangle( xy1[0], xy1[1] ) == 0 );
			}
		    }
		}
	    }
	}
    } while( crossed );
  /*OLD CODE BLOCK:
    int num;
    tPtrList< tLNode > *cr;
    tPtrListIter< tLNode > rnIter;
    tArray< double > xyz;
    for( cr = rlIter.FirstP(), i=0; !(rlIter.AtEnd());
    cr = rlIter.NextP(), i++ )
    {
    rnIter.Reset( *cr );
    num = nrnodes[i];
    for( cn = rnIter.FirstP(), j=0; j<num;
    cn = rnIter.NextP(), j++ )
    {
    xyz = cn->getNew3DCoords();
    cout << "CFC: new coords " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
    }
    }
  */
  if (0) //DEBUG
    cout << "finished" << endl;
}


/*****************************************************************************\
**
**         CheckBrokenFlowedg(): checks to see whether the flowedge of a
**            meandering node will need to be flipped; if so, check the two
**            nodes that are the "third" vertices of the triangles on either
**            side of the flowedge: if both drainage areas are less than that
**            of the meandering node, delete the closer "third" vtx. node
**            (but make sure it's not a boundary node; if it is, revert node
**            to old coord's.
**
**		Parameters:
**		Called by: Migrate
**		Created: 2/98 SL
**
\*****************************************************************************/
#define MAXLOOPS 10
void tStreamMeander::CheckBrokenFlowedg()
{
  if (0) //DEBUG
    cout << "CheckBrokenFlowedg()..." << endl;
  int nrn, nln;
  bool breakedge;
  int nloops = 0;
  double area;
  const bool flip = false;
  double dis0, dis1;
  tLNode *cn, *dn, *rn, *ln;
  tEdge *fedg, *cedg;
  tTriangle *rtri, *ltri;
  tMeshListIter< tLNode > nIter( meshPtr->getNodeList() ),
    dI( meshPtr->getNodeList() );
  tMeshListIter< tEdge > eIter( meshPtr->getEdgeList() );
  tPtrListIter< tEdge > sIter;

  /*for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
    {
    cout << "node " << cn->getID() << endl;
    if( cn->Meanders() )
    {
    cout << "   meanders" << endl;
    }
    }*/

  do
    {
      nloops++;
      breakedge = false;
      if (0) //DEBUG
	cout << "checking..." << endl;
      //look through meandering nodes:
      for( cn = nIter.FirstP(); nIter.IsActive(); cn = nIter.NextP() )
	{
	  if (0) //DEBUG
	    cout << "node " << cn->getID() << endl;
	  //for( ccn = dI.FirstP(); dI.IsActive(); ccn = dI.NextP() );
	  if( cn->Meanders() )
	    {
	      if (0) //DEBUG
		cout << "   meanders" << endl;
	      dn = cn->getDownstrmNbr();
              //if downstrm nbr exists and is still in nodeList (formerly asserted):
	      if( dn != 0 )
		{
		  if (1) //DEBUG
		    assert(dI.Get( dn->getID() ));
		  fedg = cn->getFlowEdg();
		  assert( fedg != 0 );
		  //now getting a bug: assertion in getEdgeComplement that fedg is in
		  //the edgeList failed; so I guess I need to add another firewall here
		  //to make sure it's in the list; if not, just bail 'cos we don't need
		  //to worry about breaking it if it doesn't really exist:
		  if( eIter.Get( fedg->getID() ) )
		    {
		      cedg = meshPtr->getEdgeComplement( fedg );
		      assert( cedg != 0 );
		      ltri = meshPtr->TriWithEdgePtr( fedg );
		      rtri = meshPtr->TriWithEdgePtr( cedg );
		      assert( ltri != 0 );
		      assert( rtri != 0 );
		      //if( ltri != 0 && rtri != 0 )
		      //{
		      nln = ltri->nVOp( rtri );
		      //nrn = rtri->nVOp( ltri );
		      //check for flip of flowedge (tri's on either side of flowedge:
		      if( meshPtr->CheckForFlip( ltri, nln, flip ) != FLIP_NOT_NEEDED )
			{
			  //if flowedge to be flipped, delete closer of two "third" nodes
			  //if their drareas are smaller:
			  nrn = rtri->nVOp( ltri );
			  ln = static_cast<tLNode *>( ltri->pPtr( nln ) );
			  rn = static_cast<tLNode *>( rtri->pPtr( nrn ) );
			  area = cn->getDrArea();
			  if( ln->getDrArea() < area && rn->getDrArea() < area )
			    {
			      breakedge = true;
			      dis0 = ln->DistNew( cn, dn );
			      dis1 = rn->DistNew( cn, dn );
			      //delete closer node if it's not a boundary node;
			      //if it is, revert to old coords.
			      if( dis0 < dis1 )
				{
				  if( ln->getBoundaryFlag() == kNonBoundary && !ln->isMobile())
				    {
				      meshPtr->DeleteNode( ln );
				      cn->setFlowEdg( cn->EdgToNod( dn ) );
				      if (0) //DEBUG
					cout<<"Node "<<cn->getID()<<" flows to "<<dn->getID()
					    <<" in tStreamMeander::CheckBrokenFlowedg"<<endl;
				    }
				  else cn->RevertToOldCoords();
				}
			      else
				{
				  if( rn->getBoundaryFlag() == kNonBoundary && !rn->isMobile())
				    {
				      meshPtr->DeleteNode( rn );
				      cn->setFlowEdg( cn->EdgToNod( dn ) );
				      if (0) //DEBUG
					cout<<"Node "<<cn->getID()<<" flows to "<<dn->getID()
					    <<" in tStreamMeander::CheckBrokenFlowedg"<<endl;

				    }
				  else cn->RevertToOldCoords();
				}
			    }
			}
		    }
		}
	    }
	}
      //repeat
    } while( breakedge && nloops < MAXLOOPS );
  if (0) //DEBUG
    cout << "finished CheckBrokenFlowedg" << endl;
}
#undef MAXLOOPS

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

  tMeshListIter< tLNode > dI( meshPtr->getNodeList() );
  tMeshListIter< tEdge > eI( meshPtr->getEdgeList() );
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
 **
**
**		Parameters:
**		Called by:
**		Created: 1/98 SL
**
\***************************************************************/
