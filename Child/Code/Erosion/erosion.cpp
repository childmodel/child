/***************************************************************************\
**
**  tSedTrans.cpp
**
**  Functions for tSedTrans objects.
**
**    Created 1/98 gt
**
**  $Id: erosion.cpp,v 1.1 1998-01-14 20:50:12 gtucker Exp $
\***************************************************************************/

#include "../Inclusions.h"
#include "tSedTrans.h"

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
