/*************************************************************************\
**
**  tOutput.h: Header file for output objects
**
**  NB: inherit from this basic class to do output for tLNode objects.
**
**  $Id: tOutput.h,v 1.3 1998-01-29 19:53:01 stlancas Exp $
\*************************************************************************/

#ifndef TOUTPUT_H
#define TOUTPUT_H

#include <iostream.h>
#include <string.h>
#include <assert.h>
#include "../errors/errors.h"
#include "../tGridList/tGridList.h"
#include "../GridElements/gridElements.h"
#include "../tInputFile/tInputFile.h"
#include "../tGrid/tGrid.h"

#define kMaxNameSize 80

template< class tSubNode >
class tOutput
{
public:
    tOutput( tGrid<tSubNode> * gridPtr, tInputFile &infile );
    void WriteOutput( double time );
    virtual void WriteNodeData( double time );

protected:
    tGrid<tSubNode> * g;
    char baseName[kMaxNameSize];
    ofstream nodeofs;
    ofstream edgofs;
    ofstream triofs;
    ofstream zofs;
};

template< class tSubNode >
class tLOutput : public tOutput<tSubNode>
{
public:
    tLOutput( tGrid<tSubNode> * gridPtr, tInputFile &infile );
    void WriteNodeData( double time );
private:
    ofstream drareaofs;  // Drainage areas
    ofstream netofs;     // Downstream neighbor IDs
    
};

#endif
