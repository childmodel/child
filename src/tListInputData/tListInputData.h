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
 **  $Id: tListInputData.h,v 1.28 2004-06-16 13:37:35 childcvs Exp $
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

class tRand;

#include <fstream>

#define kTimeLineMark ' '


/**************************************************************************/
/**
 **  @class tListInputDataBase
 **
 **  base services for tListInputData
 */
/**************************************************************************/
class tListInputDataBase
{
protected:
  static void findRightTime(std::ifstream &, int &, double,
			    const char *, const char *, const char *);
  static void openFile(std::ifstream &, const char *, const char *);
  // IO Error handling
  typedef enum {
    IOTime,
    IOSize,
    IORecord
  } IOErrorType;
  static void ReportIOError(IOErrorType t, const char *filename,
			    const char *suffix, int n=-1)
    ATTRIBUTE_NORETURN ;
};

/**************************************************************************/
/**
 **  @class tListInputDataMesh
 **
 **  tListInputData reads an established triangulation from a set of four
 **  files, and stores the data in a series of arrays. The files are:
 **
 **    &lt;name&gt;.nodes  --  node (point) data
 **    &lt;name&gt;.edges  --  directed edge data
 **    &lt;name&gt;.tri    --  triangle data
 **    &lt;name&gt;.z      --  "z" value data (elevation or other)
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
class tListInputDataMesh : private tListInputDataBase
{
  friend class tMesh< tSubNode >;  // gives tMesh direct access

  tListInputDataMesh();
public:
  tListInputDataMesh( const tInputFile & ); // Read filename & time from main inp file

private:
  void GetFileEntry();       // read data from files

  int nnodes, nedges, ntri;  // # nodes, edges, & triangles
  std::ifstream nodeinfile;   // node input file
  std::ifstream edgeinfile;   // edge input file
  std::ifstream triinfile;    // triangle input file
  std::ifstream zinfile;      // "z" input file

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

};

/**************************************************************************/
/**
 **  @class tListInputDataRand
 **
 **  tListInputData reads an established random number generator state.
 */
/**************************************************************************/
class tListInputDataRand : private tListInputDataBase
{
  tListInputDataRand();
public:
  tListInputDataRand( const tInputFile &, tRand & );
};

/**************************************************************************/
/**
 **  @class tListInputDataVegetation
 **
 **  tListInputData reads an established random number generator state.
 */
/**************************************************************************/
class tListInputDataVegetation : private tListInputDataBase
{
  tListInputDataVegetation();
public:
  tListInputDataVegetation( const tInputFile &);
  tArray< double > vegCov;
};

/**************************************************************************/
/**
 **  @class tListInputDataForest
 **
 **  tListInputData reads an established random number generator state.
 */
/**************************************************************************/
class tListInputDataForest : private tListInputDataBase
{
  tListInputDataForest();
public:
  tListInputDataForest( const tInputFile &);
  tArray< double > rootstrength;
  tArray< double > maxrootstrength; // root strength at time of death
  tArray< double > maxheightstand;
  tArray< double > biomassstand; // mass or volume per unit area
  tArray< double > biomassdown; // mass or volume per unit area
  tArray< double > standdeathtime;
};



/**************************************************************************\
 **
 **  tListInputDataMesh constructor
 **
 **  The constructor does most of the work. It takes an input file (i.e.,
 **  the "main" input file) and reads from it the base name of the files
 **  that contain the triangulation. It then opens <basename>.nodes,
 **  <basename>.z, <basename>.edges, and <basename>.tri. Assuming the
 **  files are valid, the desired time-slice is read from infile, and
 **  the start of data for that  time-slice is sought in each of the four
 **  triangulation files. The arrays are dimensioned as needed, and
 **  GetFileEntry() is called to read the data into the arrays. Note that
 **  the time in each file is identified by a space character preceding it
 **  on the same line.
 **
 **  Modifications:
 **   - Calls to seekg appear not to be working properly with the g++
 **     compiler, and a quick web search revealed all sorts of complaints
 **     about seekg under g++. So I recoded the file reading to avoid
 **     using seekg. (GT Feb 01)
 **   - Fixed bug in which no. edges and triangles were incorrectly
 **     assigned to nnodes, instead of nedges and ntri. (GT 04/02)
 **
\**************************************************************************/
template< class tSubNode >
tListInputDataMesh< tSubNode >::
tListInputDataMesh( const tInputFile &infile )
{
  double intime;                   // desired time
  char basename[80];               // base name of input files

  // Read base name for triangulation files from infile
  infile.ReadItem( basename, sizeof(basename), "INPUTDATAFILE" );

  // Open each of the four files
  openFile( nodeinfile, basename, SNODES);
  openFile( edgeinfile, basename, SEDGES);
  openFile( triinfile, basename, STRI);
  openFile( zinfile, basename, SZ);

  // Find out which time slice we want to extract
  intime = infile.ReadItem( intime, "INPUTTIME" );
  if (1) //DEBUG
    std::cout << "intime = " << intime << std::endl;
  if (1) //DEBUG
    std::cout << "Is node input file ok? " << nodeinfile.good()
	      << " Are we at eof? " << nodeinfile.eof() << std::endl;

  // Find specified input times in input data files and read # items.
  // First, nodes:
  findRightTime( nodeinfile, nnodes, intime,
		 basename, SNODES, "node");
  // Then elevations (or "z" values):
  findRightTime( zinfile, nnodes, intime,
		 basename, SZ, "elevation");
  // Now edges:
  findRightTime( edgeinfile, nedges, intime,
		 basename, SEDGES, "edge");
  // And finally, triangles:
  findRightTime( triinfile, ntri, intime,
		 basename, STRI, "triangle");

  // Dimension the arrays accordingly
  x.setSize( nnodes );
  y.setSize( nnodes );
  z.setSize( nnodes );
  edgid.setSize( nnodes );
  boundflag.setSize( nnodes );
  orgid.setSize( nedges );
  destid.setSize( nedges );
  nextid.setSize( nedges );
  p0.setSize( ntri );
  p1.setSize( ntri );
  p2.setSize( ntri );
  e0.setSize( ntri );
  e1.setSize( ntri );
  e2.setSize( ntri );
  t0.setSize( ntri );
  t1.setSize( ntri );
  t2.setSize( ntri );

  // Read in data from file
  GetFileEntry();

  // Close the files
  nodeinfile.close();
  edgeinfile.close();
  triinfile.close();
  zinfile.close();
}


/**************************************************************************\
 **
 **  tListInputDataMesh::GetFileEntry
 **
 **  Reads node, edge, and triangle data from the four triangulation input
 **  files. Assumes that each files is open and valid and that the current
 **  reading point in each corresponds the start of data for the desired
 **  time-slice.
 **
\**************************************************************************/
template< class tSubNode >
void tListInputDataMesh< tSubNode >::
GetFileEntry()
{
  char const * const  basename = "<file>";
  int i;

  for( i=0; i< nnodes; i++ ){
    nodeinfile >> x[i] >> y[i] >> edgid[i] >> boundflag[i];
    if (nodeinfile.fail())
      ReportIOError(IORecord, basename, SNODES, i);
    zinfile >> z[i];
    if (zinfile.fail())
      ReportIOError(IORecord, basename, SZ, i);
  }

  for( i=0; i<nedges; i++ ) {
    edgeinfile >> orgid[i] >> destid[i] >> nextid[i];
    if (edgeinfile.fail())
      ReportIOError(IORecord, basename, SEDGES, i);
  }
  for( i=0; i< ntri; i++ ) {
    triinfile >> p0[i] >> p1[i] >> p2[i] >> t0[i] >> t1[i] >> t2[i]
	      >> e0[i] >> e1[i] >> e2[i];
    if (triinfile.fail())
      ReportIOError(IORecord, basename, STRI, i);
  }
}

#endif

