/*************************************************************************\
**
**  tOutput.h: Header file for output objects
**
**  NB: inherit from this basic class to do output for tLNode objects.
**
**  $Id: tOutput.h,v 1.1 1998-01-21 01:25:34 gtucker Exp $
\*************************************************************************/

#ifndef TOUTPUT_H
#define TOUTPUT_H

#include <iostream.h>
#include "../GridElements/gridElements.h"
#include "../tInputFile/tInputFile.h"
#include "../tListInputData/tListInputData.h"
#include "../tGrid/tGrid.h"


#define kMaxNameSize 80

class tOutput
{
public:
    tOutput( tGrid<tNode> * gridPtr, tInputFile &infile );
    void WriteOutput( float time );
    virtual void WriteNodeData( float time );

private:
    tGrid<tNode> * g;
    char baseName[kMaxNameSize];
    ofstream nodeofs;
    ofstream edgofs;
    ofstream triofs;
    ofstream zofs;
};

#endif
