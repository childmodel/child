/**************************************************************************/
/**
**  @file tFloodplain.cpp
**  @brief functions for class tFloodplain and its helper class tFloodNode
**
**  Class tFloodplain simulates overbank deposition in the CHILD model,
**  using a modified version of Howard's (1996) diffusive overbank
**  deposition model.
**
**  The class provides a constructor, which reads in and initializes
**  all necessary parameters, and a DepositOverbank function that
**  implements floodplain deposition.
**
**  The overbank deposition model works as follows. Nodes with a
**  drainage area greater than a given threshold size are considered
**  large enough to initiate an overbank flood with a significant amount
**  of deposition. The minimum drainage area for a "flood node" is
**  read from the input file and recorded in drarea_min. The rate of
**  overbank deposition at any node depends on (1) its distance from the
**  nearest flood node (which is zero if it _is_ a flood node), and
**  (2) its elevation. The distance dependence is meant to simulate the
**  presumed dropoff in average sedimentation rate of suspended sand or
**  silt with distance from the main channel, and is modeled as an
**  exponential as in the original Howard model. The elevation dependence
**  is meant to simulate the increased likelihood of overbank flooding
**  in low spots versus high spots, as well as a presumed increase in
**  the total amount of sediment raining out of a deeper water column in
**  lower spots. In the original Howard model, the elevation dependence
**  was expressed as the difference between local elevation and a
**  maximum floodplain height above the channel. This was later modified
**  to use an exponential elevation weighting that accounts for a
**  distribution of event sizes. In the implementation here, we can take
**  advantage of the fact that we have an explicit distribution of event
**  sizes (i.e., the storm events that drive the model). Thus, the
**  elevation dependence is modeled as the difference between the local
**  elevation and the water surface height at the nearest flood node
**  during the current event. Water surface height is computed as the
**  sum of the flood node elevation plus water depth as computed from a
**  simple hydraulic geometry relation. Only events larger than a
**  specified precipitation rate event_min are assumed to generate
**  overbank deposition. Taken together, the deposition rate at a point
**  for a given event is (assuming P>Pe and WSH>Z):
**
**    F = MU * (WSH-Z) * exp( - Df / LAMDA )
**
**    WSH = Zf + Kdb * Pe^(mqb-mqs) * A^(mqb-mqs) * Q^mqs
**
**  where (with code variable names in parentheses):
**
**     MU (fpmu) = deposition rate per unit depth at distance = 0
**     WSH = water surface height
**     Z = elevation of node
**     Df = distance to nearest flood node
**     LAMDA (fplamda) = e-folding distance for overbank deposition
**     Zf = elevation of nearest flood node
**     Kdb = bankfull depth-discharge coefficient (ie, Db = Kdb Qb ^ mqb)
**     Pe (event_min) = assumed bankfull event precipitation rate
**     mqb = bankfull depth-discharge exponent
**     mqs = at-a-station depth-discharge exponent
**     A = drainage area of flood node
**     Q = flood node discharge
**
**  (Created 1/99 by GT)
**
**  $Id: tFloodplain.cpp,v 1.22 2004-01-29 17:09:22 childcvs Exp $
*/
/**************************************************************************/

#include "tFloodplain.h"

/**************************************************************************\
**
**  @fn tFloodplain( tInputFile &infile, tMesh<tLNode *mp );
**  @brief Main constructor for tFloodplain
**
**  @param infile Input file from which parameters are read
**  @param mp Pointer to the mesh
**
**  Takes a tInputFile as an argument and reads from it the various
**  necessary parameters. For convenience and speed, mqbmqs = mqb - mqs
**  and kdb = Kdb * Pe^(mqb-mqs).
**
\**************************************************************************/
tFloodplain::tFloodplain( const tInputFile &infile, tMesh<tLNode> *mp )
  :
  chanDriver(0),
  meshPtr(mp)
{
   int numg;

   // Keep a pointer to the mesh in order to access list of nodes
   assert( meshPtr!=0 );

   // Read in parameters
   drarea_min = infile.ReadItem( drarea_min, "FP_DRAREAMIN" );
   kdb = infile.ReadItem( kdb, "HYDR_DEP_COEFF_DS" );
   mqs = infile.ReadItem( mqs, "HYDR_DEP_EXP_STN" );
   mqbmqs = infile.ReadItem( mqbmqs, "HYDR_DEP_EXP_DS" ) - mqs;
   event_min = infile.ReadItem( event_min, "FP_BANKFULLEVENT" );
   fpmu = infile.ReadItem( fpmu, "FP_MU" );
   fplamda = infile.ReadItem( fplamda, "FP_LAMBDA" );

   kdb = kdb*pow( event_min, mqbmqs );

   //cout << "kdb: " << kdb << "  mqbmqs " << mqbmqs << endl;

   numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   deparr.setSize(numg); // dimension & init to 0 the deparr array

   {
     int tmp_;
     optControlMainChan = BOOL(infile.ReadItem( tmp_,
						"FP_OPTCONTROLCHAN" ));
   }
   if( optControlMainChan )
     chanDriver = new tMainChannelDriver( infile );

}

/**************************************************************************\
**
**  @fn ~tFloodplain
**  @brief destructor
**
\**************************************************************************/
tFloodplain::~tFloodplain()
{
  meshPtr = 0;
  delete chanDriver;
}

/**************************************************************************\
**
**  tFloodplain::DepositOverbank
**
**  Implements the overbank deposition model using the following
**  algorithm:
**
**  IF precip > event_min (ie, overbank flood event)
**    Create "flood list" of all flood nodes, computing WSH for each and
**         recording the maximum WSH
**    FOR each landscape node
**      IF node elevation is below maximum WSH
**        Find the closest flood node
**        IF node elevation < WSH at closest flood node
**          Calculate total deposition depth and update node elevation
**
**  Modifications:
**   - hydraulic geometry parameters should be interpreted for discharge
**     units of meters per second, not meters per year. Corrected. GT 3/99
**   - now uses call to EroDep with layer info. deparr is passed to
**     EroDep; all overbank deposits are assumed to be finest fraction,
**     which is assumed to be the 1st size in the array. GT 3/13/99
**   - Deposition rate WITHIN main channel now reduced to 20% of
**     nominal rate, based on the argument that (a) fine material is
**     much less likely to settle out under the swift in-channel current,
**     and (b) the coarse bedload material, which will, is typically
**     ~20% of the total load. This is not ideal, but it allows
**     in-channel deposition even under "detachment limited" conditions,
**     and is a useful approximation for large-scale floodplain sim's.
**     GT 6/99.
**
**    Parameters:
**      precip -- precipitation rate for current storm event
**      delt -- duration of current storm event
**      ctime -- current time
**    Notes:
**
\**************************************************************************/
void tFloodplain::DepositOverbank( double precip, double delt, double ctime )
{
   if( precip < event_min ) return;

   tMeshListIter<tLNode> ni( meshPtr->getNodeList() ); // iterator for nodes
   tList<tFloodNode> floodList;    // list of "flood nodes"
   tFloodNode *fn;       // ptr to current flood node
   tLNode *cn,           // current landscape node
       *closestNode;     // closest flood node
   double maxWSH = 0.0,  // maximum water surface height at any flood node
       minDist,          // minimum distance to a flood node
       dist,             // distance to flood node
       floodDepth,       // local flood depth
       dx,dy,            // x and y distance to flood node
       wsh=0.0,          // water surface height
       drarea;           // drainage area at flood node

   //cout << "tFloodplain\n";

   // Make a list of all nodes with a drainage area large enough to be
   // considered "flood generators". Record the highest water surface height.
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
      if( ( drarea = cn->getDrArea() ) >= drarea_min )
      {
	 tFloodNode floodNode(cn,
			      kdb*pow( drarea, mqbmqs )
			      *pow( cn->getQ()/SECPERYEAR, mqs )
			      + cn->getZ() );
         //cout << "flood depth " << cn->getID() << " = " << floodNode.wsh-cn->getZ() << endl;
         if( floodNode.wsh > maxWSH )
             maxWSH = floodNode.wsh;
         floodList.insertAtBack( floodNode );
      }
   }

   // Just in case there are no flood nodes, stop here
   if( floodList.isEmpty() ) return;

   // For each node, find the nearest flood node and if it's WSH is above
   // the local elevation, deposit stuff
   for( cn=ni.FirstP(); ni.IsActive(); cn=ni.NextP() )
   {
      closestNode = 0;
      minDist = kVeryFar;
      if( cn->getZ() < maxWSH ) // don't bother if node is above max flood ht
      {
         // Sweep through the flood node list to find the closest flood node,
         // recording its distance and wsh
         for( fn=floodList.FirstP(); fn!=0; fn=floodList.NextP() )
         {
            //cout << "  check fn " << fn->nodePtr->getID() << endl << flush;
            dx = cn->getX() - fn->nodePtr->getX();
            dy = cn->getY() - fn->nodePtr->getY();
            dist = sqrt( dx*dx + dy*dy );
            if( dist<minDist )
            {
               minDist = dist;
               closestNode = fn->nodePtr;
               wsh = fn->wsh;
            }

         }
         assert( closestNode!=0 ); // (should always find one)
	 assert( wsh>0 );   // should always find a closest node & set wsh
         /*cout << " got it: " << closestNode->getID() << " dist=" << minDist
           << endl << flush;*/

         // If current node is below flood level, do some deposition
         // (note: deposit is assumed to consist of 100% of the finest
         // grain size fraction, which is assumed to be the first fraction.
         // All other entries in deparr are zero)
         //   NOTE: when delt > 1/mu exp( minDist/fplamda ), the solution
         // is invalid because deposition occurs above the water surface!
         // An analytical solution znew = z + D0(1-exp(-fx delt)) could
         // be used in these cases (ie, when delt is large to speed up
         // computation)
         if( ( floodDepth = wsh - cn->getZ() ) > 0.0 )
         {
            // Depth of deposition: if in a main channel, use 20% of
            // dep rate as the bedload contribution, assuming that susp load
            // is not deposited in main channel (obviously an approximation!);
            // otherwise, use Howard formula
            if( minDist==0.0 ) deparr[0] = 0.2*floodDepth*fpmu*delt;
            else deparr[0] = floodDepth*fpmu*exp( -minDist/fplamda )*delt;
            if( deparr[0]>floodDepth)
                cout << " *WARNING, deposit thicker than flood depth\n";
            //Xcn->ChangeZ( depo ); // (note: use layering-TODO)
            cn->EroDep( 0, deparr, ctime );
            //cout << " OBDep " << deparr[0] << " at (" << cn->getX()
            //  << "," << cn->getY() << ")\n";
         }
      }
   }

}


/**************************************************************************\
**
**  @fn OptControlMainChan()
**  @brief Returns TRUE if user has requested control of main channel elev.
**
\**************************************************************************/
bool tFloodplain::OptControlMainChan() const { return optControlMainChan; }


/**************************************************************************\
**
**  @fn UpdateMainChannelHeight()
**  @brief Calls function in tMainChannelDriver to set main channel height.
**
\**************************************************************************/
void tFloodplain::UpdateMainChannelHeight( double tm, tLNode * inletNode )
{
  chanDriver->UpdateMainChannelElevation( tm, inletNode );
}



/**************************************************************************\
**
**  @fn tMainChannelDriver( tInputFile &infile )
**  @brief tMainChannelDriver constructor
**
**  @param infile Reference to the input file for the run
**
\**************************************************************************/
tMainChannelDriver::tMainChannelDriver( const tInputFile &infile )
{
  drop = infile.ReadItem( drop, "FP_VALDROP" );
  num_grnsize_fractions = infile.ReadItem( num_grnsize_fractions,
					   "NUMGRNSIZE" );
  if( num_grnsize_fractions <= 0 )
    ReportFatalError( "NUMGRNSIZE must be >= 1" );

  infile.WarnObsoleteKeyword("FP_OPTCHANALITUDEVARIATION", "FP_INLET_ELEVATION");
  infile.WarnObsoleteKeyword("FP_PERIOD", "FP_INLET_ELEVATION");
  infile.WarnObsoleteKeyword("FP_AMPLITUDE", "FP_INLET_ELEVATION");

  infile.ReadItem( InletElevationVariation, "FP_INLET_ELEVATION");
}


/**************************************************************************\
**
**  @fn UpdateMainChannelElevation( double tm )
**  @brief Updates elevations along the main stream channel.
**
**  @param tm Current time in simulation
**
**  This function sets elevations along the main channel as a function of
**  (a) time (via sinusoidal, step-wave, or user-specified variation), and
**  (b) a user-specified valley gradient (contained in the parameter
**      "drop", which is the height drop between the top and bottom of
**      the main valley.
**
**  Created: GT, May 2003
**
**  Assumes: (1) delz array elements initialized to zero;
**           (2) highest grain size fraction is the coarsest
**
\**************************************************************************/
void tMainChannelDriver::UpdateMainChannelElevation( double tm,
						     tLNode * inletNode )
{
  double newInletElevation,
    elev,
    chanslp;
  tEdge * fe;    // Ptr to flow edge of current node
  tLNode * cn;   // Current node along main channel

  assert( inletNode != 0 );

  // Update elevation at head of channel reach (inlet)
  newInletElevation = drop + InletElevationVariation.calc( tm );

  // Compute length and slope of channel
  cn = inletNode;
  double totlen = 0.0;  // Total length of channel found so far
  do
    {
      fe = cn->getFlowEdg();
      assert( fe );
      totlen += fe->getLength();
      cn = cn->getDownstrmNbr();
      assert( cn );
    }
  while( cn->getBoundaryFlag()==kNonBoundary );
  chanslp = drop / totlen;

  // Update elevations of all channel points below inlet
  tArray<double> delz( num_grnsize_fractions );
  elev = newInletElevation;
  cn = inletNode;
  do
    {
      // Set elevation of current node
      delz[num_grnsize_fractions-1] = elev - cn->getZ();
      cn->EroDep( 0, delz, tm );

      // Calculate new elevation of downstream node
      fe = cn->getFlowEdg();
      assert( fe );
      elev = elev - chanslp * fe->getLength();

      // Move to downstream node
      cn = cn->getDownstrmNbr();
      assert( cn );
    }
  while( cn->getBoundaryFlag()==kNonBoundary );


}
