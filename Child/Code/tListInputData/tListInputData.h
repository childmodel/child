/**************************************************************************\
**
**  tListInputData.h: header file for classes tListInputData and
**  tListIFStreams.
**
**  These objects are used to read in lists of Delaunay-triangulated
**  mesh elements from three user-provided input files, which contain
**  the nodes, directed edges, and triangles in the mesh, respectively.
**
**  Modifications:
**   - changed .tri file format from points-edges-triangles to
**     points-triangles-edges, compatible with earlier format (gt 1/98)
**
**  $Id: tListInputData.h,v 1.4 1998-02-01 00:57:00 stlancas Exp $
\**************************************************************************/

#ifndef TLISTINPUTDATA_H
#define TLISTINPUTDATA_H

#include <assert.h>
#include <iostream.h>
#include <string.h>
#include "../tArray/tArray.h"
#include "../tInputFile/tInputFile.h"

/** class tListIFStreams *****************************************************/
class tListIFStreams
{
  public:
      /*tListIFStreams();*/
   tListIFStreams( const char * );/*take file name from tInputFile::ReadItem*/
   tListIFStreams( tInputFile & );
   ~tListIFStreams();
   int getNNodes() const;
   int getNEdges() const;
   int getNTri() const;
   ifstream &getNodeInFile();
   ifstream &getEdgeInFile();
   ifstream &getTriInFile();
  private:
   int nnodes, nedges, ntri;
   ifstream nodeinfile;
   ifstream edgeinfile;
   ifstream triinfile;
};

/** class tListInputData ************************************************************/
template< class tSubNode >
class tListInputData
{
   friend class tGrid< tSubNode >;
  public:
   tListInputData();
   tListInputData( tListIFStreams & );
   ~tListInputData();
   void GetKeyEntry( tListIFStreams &  );
   void GetFileEntry( tListIFStreams &  );
  private:
   tArray< double > x;
   tArray< double > y;
   tArray< double > z;
   tArray< int > edgid;
   tArray< int > boundflag;
   tArray< int > orgid;
   tArray< int > destid;
   tArray< int > nextid;
   tArray< int > p0;
   tArray< int > p1;
   tArray< int > p2;
   tArray< int > e0;
   tArray< int > e1;
   tArray< int > e2;
   tArray< int > t0;
   tArray< int > t1;
   tArray< int > t2;
};

#endif
