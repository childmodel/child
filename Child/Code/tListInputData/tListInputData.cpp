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
 **   - refactoring with multiple classes (AD 11/03)
 **
 **  $Id: tListInputData.cpp,v 1.23 2003-11-14 17:59:27 childcvs Exp $
 */
/**************************************************************************/

#include "tListInputData.h"
#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif

#include "../Mathutil/mathutil.h"

void tListInputDataBase::
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
 **  tListInputDataBase::openFile
 **
 **  Find the right time in file and position it for reading
 **
\**************************************************************************/
void tListInputDataBase::
openFile( ifstream &infile, const char *basename,
	  const char *ext)
{
  char inname[80];                  // full name of an input file

  // Open each of the four files
  assert( strlen(basename)+strlen(ext)<sizeof(inname) );
  strcpy( inname, basename );
  strcat( inname, ext );
  infile.open(inname);
  if( !infile.good() )
    {
      cerr << "Error: I can't find the following files:\n"
           << "\t" << basename << ext << "\n";
      ReportFatalError( "Unable to open triangulation input file(s)." );
    }
}

/**************************************************************************\
 **
 **  tListInputData::findRightTime
 **
 **  Find the right time in file and position it for reading
 **
\**************************************************************************/
void tListInputDataBase::
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
 **  tListInputRand::tListInputDataRand()
 **
 **  Read state from file
 **
\**************************************************************************/
tListInputDataRand::
tListInputDataRand( const tInputFile &infile, tRand &rand )
{
  double intime;                   // desired time
  char basename[80];

  // Read base name for triangulation files from infile
  infile.ReadItem( basename, sizeof(basename), "INPUTDATAFILE" );

  // Open each of the four files
  ifstream randominfile;
  openFile( randominfile, basename, SRANDOM);

  // Find out which time slice we want to extract
  intime = infile.ReadItem( intime, "INPUTTIME" );
  if (1) //DEBUG
    cout << "intime = " << intime << endl;

  // Find specified input times in input data files and read # items.
 // And finally, random number generator:
  int nrandom;
  findRightTime( randominfile, nrandom, intime,
		 basename, SRANDOM, "random number generator");
  if ( rand.numberRecords() != nrandom ) {
    cerr << "Invalid number of records for the random number generator\n";
    ReportFatalError( "Input error" );
  }

  // Read in data from file
  rand.readFromFile( randominfile );
}
