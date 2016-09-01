/**************************************************************************/
/**
**  @file tVegetation.cpp
**  @brief Functions for tVegetation and tVegCover classes
**
**  Classes tVegetation and tVegCover represents the vegetation cover
**  across a terrain. Class tVegCover represents the properties of
**  vegetation at a point on the landscape (e.g., the %cover), while
**  class tVegetation represents the landscape-wide ("global")
**  properties (e.g., regrowth coefficient). Class tVegCover is designed
**  to be embedded in a node (in other words, each node "has a"
**  vegetation cover).
**
**  Created January, 2000, GT
**  State read back from file if needed - November 2003, AD
**
**  $Id: tVegetation.cpp,v 1.16 2004/05/10 10:52:52 childcvs Exp $
*/
/**************************************************************************/

#include <assert.h>
#include "tVegetation.h"
#include "../tMesh/tMesh.h"

/**************************************************************************\
**                 FUNCTIONS FOR CLASS tVegetation
\**************************************************************************/

/**************************************************************************\
**
**  tVegetation constructors:
**    1. Default
**    2. Input file
**
\**************************************************************************/

tVegetation::tVegetation()
  :
  mdKvd(0),
  mdTVeg(1),
  mdTauCritBare(0), mdTauCritVeg(0)
{}

tVegetation::tVegetation( tMesh<class tLNode> * meshPtr, const tInputFile &infile )
  :
  mdKvd(0),
  mdTVeg(1),
  mdTauCritBare(0), mdTauCritVeg(0)
{
   mdKvd = infile.ReadItem( mdKvd, "VEG_KVD" );
   mdTVeg = infile.ReadItem( mdTVeg, "VEG_TV" );
   mdTauCritBare = infile.ReadItem( mdTauCritBare, "TAUC" );
   mdTauCritVeg = infile.ReadItem( mdTVeg, "VEG_TAUCVEG" );

   // Loop through nodes and set initial vegetation cover & threshold
   // (for now, assume constant; later need to add restart capability)
   // Note: assumes initially 100% cover.
   tMesh<tLNode>::nodeListIter_t niter( meshPtr->getNodeList() );

   // unused
   // intlVegCover = infile.ReadItem( intlVegCover, "INTLVEGCOV" );

   // Initialise nodes
   int opt;
   if ( (opt = infile.ReadItem( opt, "OPTREADINPUT" ))
	== OPTREADINPUT_PREVIOUS) {
     // Start from a previous computation
     tListInputDataVegetation inputVegData( infile );
     const int nnodes = meshPtr->getNodeList()->getSize();
     if (inputVegData.vegCov.getSize() != static_cast<size_t>(nnodes))
       ReportFatalError( "tVegetation(): invalid number of records"
			 " in input file." );
     // Rely on the fact that the IDs have not been re-numbered.
     // for fast lookup per ID
     const tMesh< tLNode >::tIdArrayNode_t NodeTable(*(meshPtr->getNodeList()));
     for( int id=0; id < nnodes; ++id )
       NodeTable[id]->getVegCover().mdVeg = inputVegData.vegCov[id];
     for( tLNode *cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
       cn->setTauCrit( mdTauCritBare + mdTauCritVeg );
   } else {
     // Start from scratch
     for( tLNode *cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
       {
	 cn->getVegCover().mdVeg = 1.0;
	 cn->setTauCrit( mdTauCritBare + mdTauCritVeg );
       }
   }
}

/**************************************************************************\
**
**  tVegetation::UpdateVegetation
**
**  This routine updates the % vegetation cover at each node in response
**  to a storm event and its ensuing interstorm period. It also updates
**  the critical shear stress at each node according to its present
**  vegetation cover.
**
**  Erosion of vegetation during a storm is computed as:
**    dV/dt = -Kvd V ( tau - tauc )
**  Where V represents the proportional cover (0 to 1), tau is shear
**  stress, tauc is critical shear stress, and Kvd is a vegetation
**  erodibility coefficient. The equation is solved using an analytical
**  solution, using the expedient approximation that tauc is constant
**  during a given storm.
**
**  Regrowth of vegetation following a storm is computed as:
**    dV/dt = 1/Tv ( 1 - V )
**  where Tv is the timescale of vegetation regrowth. An analytical
**  solution is used to update the vegetation cover.
**
**  Finally, the critical shear stress is updated at each node using:
**    Tc = Tcb + V Tcv
**  where Tcb is critical shear stress in the absence of cover (bare)
**  and Tcv is critical shear stress under 100% cover.
**
**  Created January 2000, GT
**
\**************************************************************************/

void tVegetation::UpdateVegetation( tMesh<class tLNode> *meshPtr,
				    double dt,
                                    double interstormdur ) const
{
  ErodeVegetation( meshPtr, dt );
  GrowVegetation( meshPtr, interstormdur );

}


void tVegetation::ErodeVegetation(  tMesh<class tLNode> *meshPtr,
				    double dt ) const
{
   tMesh<tLNode>::nodeListIter_t niter( meshPtr->getNodeList() ); // Node iterator
   tLNode * cn;   // Ptr to current node
   double tauex,  // Excess shear stress
       veg;       // Fractional vegetation cover
   
   // Loop on active nodes, computing erosion during the time step of
   // duration dt.
   //   For both erosion and regrowth, we use an analytical solution for
   // veg cover given its initial value and duration of erosion/regrowth.
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
      // Erosion of vegetation during storm (if any)
      tauex = cn->getTau() - cn->getTauCrit();
      veg = cn->getVegCover().getVeg();
      if( tauex>0.0 )
      {
          veg = veg * exp( -mdKvd * tauex * dt );
          cn->getVegCover().mdVeg = veg;
	  assert( veg >= 0.0 );
	  assert( veg <= 1.0 );
          //cout << "veg after erosion: " << veg << endl;
      }

      // Update critical shear stress
      cn->setTauCrit( mdTauCritBare + veg*mdTauCritVeg );
      //cout << "tau crit: " << mdTauCritBare + veg*mdTauCritVeg << endl;
   }
   
}


void tVegetation::GrowVegetation(  tMesh<class tLNode> *meshPtr,
				    double interstormdur ) const
{
   tMesh<tLNode>::nodeListIter_t niter( meshPtr->getNodeList() ); // Node iterator
   tLNode * cn;   // Ptr to current node
   double veg;       // Fractional vegetation cover
   
   // Loop on active nodes, computing regrowth during the interstorm period.
   //   For both erosion and regrowth, we use an analytical solution for
   // veg cover given its initial value and duration of erosion/regrowth.
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
      // Regrowth following storm
      //cout << "veg before regrowth: " << veg << endl;
      veg = cn->getVegCover().getVeg();
      veg = 1.0 - (1.0 - veg) * exp( -interstormdur / mdTVeg );
      cn->getVegCover().mdVeg = veg;
      //cout << "veg after regrowth: " << veg << endl;
      assert( veg >= 0.0 );
      assert( veg <= 1.0 );

      // Update critical shear stress
      cn->setTauCrit( mdTauCritBare + veg*mdTauCritVeg );
      //cout << "tau crit: " << mdTauCritBare + veg*mdTauCritVeg << endl;
   }
}

