/*************************************************************************\
**
**  tOutput.h: Header file for output objects
**
**  NB: inherit from this basic class to do output for tLNode objects.
**
**  $Id: tOutput.h,v 1.2 1998-01-27 23:41:39 gtucker Exp $
\*************************************************************************/

#ifndef TOUTPUT_H
#define TOUTPUT_H

#include <iostream.h>
#include "../tInputFile/tInputFile.h"
#include "../tGrid/tGrid.h"

#define kMaxNameSize 80

template< class tSubNode >
class tOutput
{
public:
    tOutput( tGrid<tSubNode> * gridPtr, tInputFile &infile );
    void WriteOutput( float time );
    virtual void WriteNodeData( float time );

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
    void WriteNodeData( float time );
private:
    ofstream drareaofs;  // Drainage areas
    ofstream netofs;     // Downstream neighbor IDs
    
};

#endif
