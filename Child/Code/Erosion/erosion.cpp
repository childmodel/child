/***************************************************************************\
**
**  erosion.cpp
**
**  Functions for sediment transport and bed erosion (detachment) objects.
**  Transport objects:
**    tSedTransPwrLaw
**  Detachment objects:
**    tBedErodePwrLaw
**
**    Created 1/98 gt
**
**  $Id: erosion.cpp,v 1.2 1998-01-15 19:34:31 gtucker Exp $
\***************************************************************************/

#include "../Inclusions.h"
#include "erosion.h"

/***************************************************************************\
**  Constructors:
**    Given a tInputFile as an argument, will read relevant parameters from
**  the input file.
\***************************************************************************/
tSedTransPwrLaw::tSedTransPwrLaw( tInputFile &infile )
{
   kf = infile.ReadItem( kf, "KF" );
   mf = infile.ReadItem( mf, "MF" );
   nf = infile.ReadItem( nf, "NF" );
}

//tBedErodePwrLaw::tBedErodePwrLaw()
//{ kb = 0; mb=0; nb=0; }

tBedErodePwrLaw::tBedErodePwrLaw( tInputFile &infile )
{
   kb = infile.ReadItem( kb, "KB" );
   mb = infile.ReadItem( mb, "MB" );
   nb = infile.ReadItem( nb, "NB" );
}



/***************************************************************************\
**  tBedErodePwrLaw functions
\***************************************************************************/

/***************************************************************************\
**  tBedErode::DetachCapacity
**
**  Computes the depth of erosion over a time interval dt assuming the
**  erosion rate = kb Q^mb S^nb
**
**  Input: n -- node at which to compute detachment capacity
**         dt -- time interval
**  Returns: the detachment depth
**  Assumptions: n->GetSlope() does not return a negative value; kb, mb,
**               and nb all >=0.
\***************************************************************************/
float tBedErodePwrLaw::DetachCapacity( tLNode * n, float dt )
{
   float slp = n->GetSlope();
   return( kb*pow( n->GetQ(), mb )*pow( slp, nb )*dt );
}

/***************************************************************************\
**  tBedErode::SetTimeStep
**
**  Estimates a maximum time step as a function of the Courant stability
**  criterion: dt <= dx/v, where dx is the node spacing and v is the
**  wave velocity. If nb=1, v = kb Q^mb. If the equation is nonlinear
**  with nb!=1, the wave velocity is estimated as if the equation were
**  linear with a gradient term in the coefficient, ie
**      v S = [kb Q^mb S^(nb-1)] S  (recall S = -dz/dx)
**  The estimated time step limitation is therefore given by:
**      dt = 0.2 * ( dx / (kb Q^mb S^(nb-1) ) )
**  (the 0.2 is a factor to keep us comfortably below the Courant #).
**
**  Input: n -- node for which to estimate time step
**  Returns: the estimated maximum time step size
**  Assumptions: GetSlope() returns a value >=0, edge length>0.
**
\***************************************************************************/
float tBedErodePwrLaw::SetTimeStep( tLNode * n )
{
   return( 0.2*n->GetFlowEdg()->getLength() /
           ( kb * pow( n->GetQ(), mb ) * pow( n->GetSlope(), nb-1.0 ) ) );
}
