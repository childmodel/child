/**************************************************************************\
**
**  tVegetation.cpp: Functions for tVegetation and tVegCover classes
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
**  
**  $Id: tVegetation.cpp,v 1.5 2002-04-30 17:17:27 arnaud Exp $
\**************************************************************************/

#include <assert.h>
#include "tVegetation.h"
#include "../tMesh/tMesh.h"
#include "../tMeshList/tMeshList.h"

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

tVegetation::tVegetation( tMesh<tLNode> * meshPtr, tInputFile &infile )
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
   tMeshListIter<tLNode> niter( meshPtr->getNodeList() );
   tLNode * cn;   

   // unused
   // intlVegCover = infile.ReadItem( intlVegCover, "INTLVEGCOV" );

   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
       cn->getVegCover().mdVeg = 1.0;
       cn->setTauCrit( mdTauCritBare + mdTauCritVeg );
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

void tVegetation::UpdateVegetation( tMesh<tLNode> *meshPtr, double stormdur,
                                    double interstormdur )
{
   tMeshListIter< tLNode > niter( meshPtr->getNodeList() ); // Node iterator
   tLNode * cn;   // Ptr to current node
   double tauex,  // Excess shear stress
       veg;       // Fractional vegetation cover
   
   // Loop on active nodes, computing erosion from the previous storm
   // and regrowth during the subsequent interstorm period.
   // Key assumptions: we "cheat" by changing the critical shear stress only
   // "after" the storm, using the initial value to calculate the
   // subsequent, eroded value, rather than solving the full quadratic
   // equation for vegetation destruction.
   //   For both erosion and regrowth, we use an analytical solution for
   // veg cover given its initial value and duration of erosion/regrowth.
   for( cn=niter.FirstP(); niter.IsActive(); cn=niter.NextP() )
   {
      // Erosion of vegetation during storm (if any)
      tauex = cn->getTau() - cn->getTauCrit();
      veg = cn->getVegCover().getVeg();
      if( tauex>0.0 )
      {
          veg = veg * exp( -mdKvd * tauex * stormdur );
          cn->getVegCover().mdVeg = veg;
          //cout << "veg after erosion: " << veg << endl;
      }
      
      // Regrowth following storm
      //cout << "veg before regrowth: " << veg << endl;
      veg = 1.0 - (1.0 - veg) * exp( -interstormdur / mdTVeg );
      cn->getVegCover().mdVeg = veg;
      //cout << "veg after regrowth: " << veg << endl;

      // Update critical shear stress
      cn->setTauCrit( mdTauCritBare + veg*mdTauCritVeg );
      //cout << "tau crit: " << mdTauCritBare + veg*mdTauCritVeg << endl;
   }
   
}



