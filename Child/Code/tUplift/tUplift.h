/************************************************************************\
**
**  tUplift.h:  Header file for tUplift objects.
**
**  tUplift objects manage data and functions related to all types of
**  deformation and baselevel (including sea level) change. Thus, the
**  name "tUplift" is a bit misleading in that the relative motions that
**  the tUplift class models need not be solely vertical, but could for
**  example include strike-slip deformation. Moreover, the class
**  model _relative_ motion, which of course need not be of tectonic
**  origin but could instead reflect sealevel or other baselevel changes.
**  Since we have no good term that encompasses all these possibilities,
**  "tUplift" seems the most convenient catch-all.
**
**  The constructor reads in a code for the type of boundary
**  condition desired, and then reads in parameters appropriate to that
**  type. The "uplift" is then applied via a general DoUplift()
**  function, which in turn calls the appropriate function (e.g.,
**  UpliftUniform()) to implement the desired behavior.
**
**  $Id: tUplift.h,v 1.2 1998-02-24 01:42:01 stlancas Exp $
\************************************************************************/
#ifndef TUPLIFT_H
#define TUPLIFT_H
#include <iostream.h>
#include "../tInputFile/tInputFile.h"
#include "../tGrid/tGrid.h"


class tUplift
{
public:
    tUplift( tInputFile &infile );
    void DoUplift( tGrid<tLNode> *gp, double delt );
    double GetDuration();
    void UpliftUniform( tGrid<tLNode> *gp, double delt );
   double GetRate() const;

private:
    int typeCode;
    double duration;
    double rate;

};

#endif








