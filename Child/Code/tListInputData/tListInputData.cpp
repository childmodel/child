/**************************************************************************/
/**
 **  @file tListInputData.cpp
 **  @brief Functions for class tListInputData.
 **
 **  Modifications:
 **   - changed .tri file format from points-edges-triangles to
 **     points-triangles-edges, compatible with earlier format (gt 1/98)
 **   - GT merged tListIFStreams and tListInputData into a single class
 **     to avoid multiple definition errors resulting from mixing
 **     template & non-template classes (1/99)
 **   - Bug fix in constructor: nnodes was being read from edge and
 **     triangle files -- thus arrays dimensioned incorrectly! (GT 04/02)
 **   - Remove dead code. Add findRightTime (AD 07/03)
 **   - Add Random number generator handling. (AD 08/03)
 **
 **  $Id: tListInputData.cpp,v 1.20 2003-09-02 08:44:33 childcvs Exp $
 */
/**************************************************************************/

#include "tListInputData.h"
#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif

template< class tSubNode >
void tListInputData< tSubNode >::
ReportIOError(IOErrorType t, const char *filename,
	      const char *suffix, int n) {
  cerr << "\nFile: '" << filename << suffix << "' "
       << "- Can't read ";
  switch (t){
  case IOTime:
    cerr << "time";
    break;
  case IOSize:
    cerr << "size";
    break;
  case IORecord:
    cerr << "record " << n;
    break;
  }
  cerr << "." << endl;
  ReportFatalError( "Input/Output Error." );
}


/**************************************************************************\
 **
 **  tListInputData constructor
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
tListInputData< tSubNode >::
tListInputData( tInputFile &infile, tRand &rand )          //tListInputData
{
  double intime;                   // desired time
  char basename[80],               // base name of input files
    inname[80];                  // full name of an input file

  // Read base name for triangulation files from infile
  infile.ReadItem( basename, sizeof(basename), "INPUTDATAFILE" );

  // Open each of the four files
  strcpy( inname, basename );
  strcat( inname, SNODES );
  nodeinfile.open(inname);    // Node input file pointer
  //assert( nodeinfile.good() );
  strcpy( inname, basename );
  strcat( inname, SEDGES );
  edgeinfile.open(inname);    // Edge input file pointer
  //assert( edgeinfile.good() );
  strcpy( inname, basename );
  strcat( inname, STRI );
  triinfile.open( inname );   // Triangle input file pointer
  //assert( triinfile.good() );
  strcpy( inname, basename );
  strcat( inname, SZ );
  zinfile.open( inname );     // Elevations input file pointer
  //assert( zinfile.good() );
  strcpy( inname, basename );
  strcat( inname, SRANDOM );
  randominfile.open( inname );// Random number generator input file pointer

  // Make sure we found them
  if( !nodeinfile.good() || !edgeinfile.good() || !triinfile.good()
      || !zinfile.good() || !randominfile.good() )
    {
      cerr << "Error: I can't find one or more of the following files:\n"
           << "\t" << basename << SNODES "\n"
           << "\t" << basename << SEDGES "\n"
           << "\t" << basename << STRI "\n"
           << "\t" << basename << SZ "\n"
           << "\t" << basename << SRANDOM "\n";
      ReportFatalError( "Unable to open triangulation input file(s)." );
    }

  // Find out which time slice we want to extract
  intime = infile.ReadItem( intime, "INPUTTIME" );
  cout << "intime = " << intime << endl;
  cout << "Is node input file ok? " << nodeinfile.good()
       << " Are we at eof? " << nodeinfile.eof() << endl;

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
 // And finally, random number generator:
  int nrandom;
  findRightTime( randominfile, nrandom, intime,
		 basename, SRANDOM, "random number generator");
  if ( rand.numberRecords() != nrandom ) {
    cerr << "Invalid number of records for the random number generator\n";
    ReportFatalError( "Input error" );
  }

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
  GetFileEntry(rand);

  // Close the files
  nodeinfile.close();
  edgeinfile.close();
  triinfile.close();
  zinfile.close();
  randominfile.close();
}


/**************************************************************************\
 **
 **  tListInputData::findRightTime
 **
 **  Find the right time in file and position it for reading
 **
\**************************************************************************/
template< class tSubNode >
void tListInputData< tSubNode >::
findRightTime( ifstream &infile, int &nn, double intime,
	       const char *basename, const char *ext, const char *typefile)
{
  char headerLine[kMaxNameLength]; // header line read from input file
  bool righttime = false;
  double time;
  while( !( infile.eof() ) && !righttime )
    {
      /*infile.getline( headerLine, kMaxNameLength );
	if( headerLine[0] == kTimeLineMark )
	{
	infile.seekg( -infile.gcount(), ios::cur );
	infile >> time;
	cout << "from file, time = " << time << endl;
	cout << "Read: " << headerLine << endl;
	if( time == intime ) righttime = 1;
	}*/
      infile >> time;
      if (infile.fail())
	ReportIOError(IOTime, basename, ext);
      if (0) //DEBUG
	cout << "Read time: " << time << endl;
      if( time < intime )
	{
	  infile >> nn;
	  if (infile.fail())
	    ReportIOError(IOSize, basename, ext);
	  if (0) //DEBUG
	    cout << "nn (" << typefile << ")= " << nn << endl;
	  int i;
	  for( i=1; i<=nn+1; i++ ) {
	    infile.getline( headerLine, kMaxNameLength );
	  }
	}
      else righttime = true;
      if (0) //DEBUG
	cout << " NOW are we at eof? " << infile.eof() << endl;
    }
  if( !( infile.eof() ) ) {
    infile >> nn;
    if (infile.fail())
      ReportIOError(IOSize, basename, ext);
  } else {
    cerr << "Couldn't find the specified input time in the " << typefile
	 << " file\n";
    ReportFatalError( "Input error" );
  }
}


/**************************************************************************\
 **
 **  tListInputData::GetFileEntry
 **
 **  Reads node, edge, and triangle data from the four triangulation input
 **  files. Assumes that each files is open and valid and that the current
 **  reading point in each corresponds the start of data for the desired
 **  time-slice.
 **
\**************************************************************************/
template< class tSubNode >
void tListInputData< tSubNode >::
GetFileEntry(tRand &rand)                  //tListInputData
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
  rand.readFromFile( randominfile );
}


/**************************************************************************\
 **
 **  tListInputData::GetKeyEntry
 **
 **  Provides alternative keyboard entry of triangulation data for
 **  testing purposes. Not currently supported.
 **
\**************************************************************************/
template< class tSubNode >
void tListInputData< tSubNode >::
GetKeyEntry()                   //tListInputData
{
  int i;
  for( i=0; i < nnodes; i++ )
    {
      cout << "x y z edgid boundary:" << endl;
      cin >> x[i] >> y[i] >> z[i] >> edgid[i] >> boundflag[i];
    }
  for( i=0; i < nedges; i++ )
    {
      cout << "orgid destid nextid" << endl;
      cin >> orgid[i] >> destid[i] >> nextid[i];
    }
  for( i=0; i< ntri; i++ )
    {
      cout << "nodeids (3), edgids (3), triangleids (3)" << endl;
      cin >> p0[i] >> p1[i] >> p2[i]
          >> e0[i] >> e1[i] >> e2[i]
          >> t0[i] >> t1[i] >> t2[i];
    }

}
