/**************************************************************************\
**
**  tStreamNet.cpp
**
**  Functions for class tStreamNet.
**
**  $Id: tStreamNet.cpp,v 1.3 1998-01-15 19:36:03 gtucker Exp $
\**************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <math.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "../GlobalFns.h"
#include "../tArray/tArray.h"
#include "../tPtrListNode/tPtrListNode.h"
#include "../tPtrList/tPtrList.h"
#include "../tPtrListIter/tPtrListIter.h"
#include "../tNode/tNode.h"
#include "../tEdge/tEdge.h"
#include "../tTriangle/tTriangle.h"
#include "../tListNode/tListNode.h"
#include "../tList/tList.h"
      for( i=0, cn=ni.FirstP(); i<nActNodes; i++, cn=ni.NextP() )
          cn->EroDep( dz[i] );
      dtg -= dt;  // decrease remaining time
   } while( dtg>0 );
   
}

   
