//-*-c++-*-

/**************************************************************************/
/**
 **  @file tListInputData.h
 **  @brief header file for class tListInputData
 **
 **  This class is used to read in lists of Delaunay-triangulated
 **  mesh elements from three user-provided input files, which contain
 **  the nodes, directed edges, and triangles in the mesh, respectively.
 **
 **  Modifications:
 **   - changed .tri file format from points-edges-triangles to
 **     points-triangles-edges, compatible with earlier format (gt 1/98)
 **   - GT merged tListIFStreams and tListInputData into a single class
 **     to avoid multiple definition errors resulting from mixing
 **     template & non-template classes (1/99)
 **
 **  $Id: tListInputData.h,v 1.21 2003-08-01 17:14:55 childcvs Exp $
 */
/**************************************************************************/

#ifndef TLISTINPUTDATA_H
#define TLISTINPUTDATA_H

#include <assert.h>
#include <string.h>
#include "../Classes.h"
#include "../Definitions.h"
#include "../errors/errors.h"
#include "../tArray/tArray.h"
#include "../tInputFile/tInputFile.h"

#if !defined(HAVE_NO_NAMESPACE)
# include <fstream>
using namespace std;
#else
# include <fstream.h>
#endif

#define kTimeLineMark ' '


/**************************************************************************/
/**
 **  @class tListInputData
 **
 **  tListInputData reads an established triangulation from a set of four
 **  files, and stores the data in a series of arrays. The files are:
 **
 **    <name>.nodes  --  node (point) data
 **    <name>.edges  --  directed edge data
 **    <name>.tri    --  triangle data
 **    <name>.z      --  "z" value data (elevation or other)
 **  The files are ASCII text, and use the same format as the corresponding
 **  output files (see tOutput.h/.cpp), and in fact can be output files
 **  from a previous run. Each file contains triangulation data for one
 **  or more time-slices, with each time slice proceeded by two values:
 **  the time and the number of elements (nodes, edges, or triangles,
 **  respectively for the .nodes, .edges, and .tri files, and nodes for
 **  the .z file). The time value must be preceded by a space. Following
 **  the time and # of elements are the elements themselves. Examples of
 **  the format for each file are shown below. Note that time is only
 **  meaningful if the files are output files from a previous run; if the
 **  files represent a triangulation obtained from another source (e.g.,
 **  a DEM), the time value can simply be set to zero (and in this case
 **  the files would normally contain only one "time slice" of course).
 **  In each file the elements are listed in order by ID number, starting
 **  from zero (the ID #s are not actually list in the files).
 **
 **  NODE file format:  (note: the _ character represents a space)
 **    _time
 **    #nodes
 **    x0 y0 edg0 bnd0
 **    x1 y1 edg1 bnd1
 **    ...etc.
 **      where x0 and y0 are the x and y coords of node 0, edg0 is the
 **      ID number of one of the directed edges that radiates from node 0,
 **      and bnd0 is the boundary code (0, 1, or 2) for node 0, etc.
 **
 **  Z file format:
 **    _time
 **    #nodes
 **    z0
 **    z1
 **    ...etc.
 **      where z0 is the z value (e.g., elevation) for node 0, etc.
 **
 **  EDGE file format:
 **    _time
 **    #edges
 **    orgid0 destid0 nextid0
 **    orgid1 destid1 nextid1
 **    ...etc.
 **      where orgid0 and destid0 are the ID numbers of the origin and
 **      destination nodes, respectively, and nextid0 is the ID # of the
 **      neighboring edge in the counterclockwise direction. Note that
 **      complementary edge pairs must be stored together, so that
 **      orgid0=destid1 and destid0=orgid1, etc.
 **
 **  TRI file format:
 **    _time
 **    #triangles
 **    p00 p01 p02 t00 t01 t02 e00 e01 e02
 **    p10 p11 p12 t10 t11 t12 e10 e11 e12
 **    ...etc.
 **      where p00, p01, p02 are the ID #s of the 3 nodes (points) in
 **      triangle 0 (in counter-clockwise order); e00, e01, and e02 are
 **      the ID #s of the 3 clockwise-oriented directed edges in triangle 0
 **      (with origin nodes p00, p01, and p02, respectively); and
 **      t00, t01, t02 are the ID #s of the 3 adjacent triangles that
 **      lie opposite nodes p00, p01, and p02, respectively (a value of -1
 **      indicates that there is no opposite triangle across a given face).
 **
 **  The class includes arrays to store all of these data as well as file
 **  streams for each input file and the # of nodes, edges, and triangles.
 **  A key entry function is provided, but is not supported in this version.
 **
 **  Note that the class is templated only because of its friendship with
 **  tMesh.
 **
 */
/**************************************************************************/
template< class tSubNode >
class tListInputData
{
  friend class tMesh< tSubNode >;  // gives tMesh direct access

public:
  tListInputData( tInputFile &, tRand & ); // Read filename & time from main inp file

private:
  static void findRightTime(ifstream &, int &, double,
			    const char *, const char *, const char *);
  void GetKeyEntry();        // not currently supported
  void GetFileEntry(tRand &);       // read data from files

  int nnodes, nedges, ntri;  // # nodes, edges, & triangles
  ifstream nodeinfile;   // node input file
  ifstream edgeinfile;   // edge input file
  ifstream triinfile;    // triangle input file
  ifstream zinfile;      // "z" input file
  ifstream randominfile; // random generator input file

  tArray< double > x;      // node x coords
  tArray< double > y;      // node y coords
  tArray< double > z;      // node z values
  tArray< int > edgid;     // node edge ID #s
  tArray< int > boundflag; // node boundary codes
  tArray< int > orgid;     // directed edge origin node ID #s
  tArray< int > destid;    // directed edge destination node ID #s
  tArray< int > nextid;    // ID #s of next counter-clockwise edges
  tArray< int > p0;     // IDs of triangle node 0
  tArray< int > p1;     // IDs of triangle node 1
  tArray< int > p2;     // IDs of triangle node 2
  tArray< int > e0;     // IDs triangle clockwise-oriented edge 0
  tArray< int > e1;     // IDs triangle clockwise-oriented edge 1
  tArray< int > e2;     // IDs triangle clockwise-oriented edge 2
  tArray< int > t0;     // IDs of neighboring tri's opposite node 0
  tArray< int > t1;     // IDs of neighboring tri's opposite node 1
  tArray< int > t2;     // IDs of neighboring tri's opposite node 2

  // IO Error handling
  typedef enum {
    IOTime,
    IOSize,
    IORecord
  } IOErrorType;
  static void ReportIOError(IOErrorType t, const char *filename,
			    const char *suffix, int n=-1);
};


/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tListInputData.cpp"
#endif

#endif

