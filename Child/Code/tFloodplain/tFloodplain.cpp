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
**  $Id: tFloodplain.cpp,v 1.11 2003-01-17 17:30:27 childcvs Exp $
*/
/**************************************************************************/

#include "tFloodplain.h"

/**************************************************************************\
**
**  Constructor:
**
**  Takes a tInputFile as an argument and reads from it the various
**  necessary parameters. For convenience and speed, mqbmqs = mqb - mqs
**  and kdb = Kdb * Pe^(mqb-mqs).
**
\**************************************************************************/
tFloodplain::tFloodplain( tInputFile &infile, tMesh<tLNode> *mp )
  :
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
#define kYearpersec 3.171e-8 // 1/SecondsPerYear
void tFloodplain::DepositOverbank( double precip, double delt, double ctime )
{
   if( precip < event_min ) return;
   
   tMeshListIter<tLNode> ni( meshPtr->getNodeList() ); // iterator for nodes
   tList<tFloodNode> floodList;    // list of "flood nodes"
   tFloodNode floodNode,           // flood node to be added to list
       *fn;                        // ptr to current flood node
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
         floodNode.nodePtr = cn;
         floodNode.wsh = kdb*pow( drarea, mqbmqs )
             *pow( cn->getQ()*kYearpersec, mqs )
             + cn->getZ();
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


