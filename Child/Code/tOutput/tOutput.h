/*************************************************************************\
**
**  tOutput.h: Header file for output objects
**
**  NB: inherit from this basic class to do output for tLNode objects.
**
**  $Id: tOutput.h,v 1.4 1998-02-02 17:53:07 gtucker Exp $
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
    void CreateAndOpenFile( ofstream * theOFStream, char * extension );

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
    ofstream slpofs;     // Slopes in the direction of flow
    
};

#endif
