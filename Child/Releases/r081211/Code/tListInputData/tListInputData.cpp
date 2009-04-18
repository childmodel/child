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
 **   - Refactoring with multiple classes (AD 11/03)
 **   - Add tVegetation handling (AD 11/03)
 **
 **  $Id: tListInputData.cpp,v 1.25 2004/06/16 13:37:35 childcvs Exp $
 */
/**************************************************************************/

#include "tListInputData.h"
#include <iostream>

#include "../Mathutil/mathutil.h"

void tListInputDataBase::
ReportIOError(IOErrorType t, const char *filename,
	      const char *suffix, int n) {
  std::cerr << "\nFile: '" << filename << suffix << "' "
       << "- Can't read ";
  switch (t){
  case IOTime:
    std::cerr << "time";
    break;
  case IOSize:
    std::cerr << "size";
    break;
  case IORecord:
    std::cerr << "record " << n;
    break;
  }
  std::cerr << "." << std::endl;
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
openFile( std::ifstream &infile, const char *basename,
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
      std::cerr << "Error: I can't find the following files:\n"
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
findRightTime( std::ifstream &infile, int &nn, double intime,
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
	std::cout << "from file, time = " << time << std::endl;
	std::cout << "Read: " << headerLine << std::endl;
	if( time == intime ) righttime = 1;
	}*/
      infile >> time;
      if (infile.fail())
	ReportIOError(IOTime, basename, ext);
      if (0) //DEBUG
	std::cout << "Read time: " << time << std::endl;
      if( time < intime )
	{
	  infile >> nn;
	  if (infile.fail())
	    ReportIOError(IOSize, basename, ext);
	  if (0) //DEBUG
	    std::cout << "nn (" << typefile << ")= " << nn << std::endl;
	  int i;
	  for( i=1; i<=nn+1; i++ ) {
	    infile.getline( headerLine, kMaxNameLength );
	  }
	}
      else righttime = true;
      if (0) //DEBUG
	std::cout << " NOW are we at eof? " << infile.eof() << std::endl;
    }
  if( !( infile.eof() ) ) {
    infile >> nn;
    if (infile.fail())
      ReportIOError(IOSize, basename, ext);
  } else {
    std::cerr << "Couldn't find the specified input time in the " << typefile
	 << " file\n";
    ReportFatalError( "Input error" );
  }
}


/**************************************************************************\
 **
 **  tListInputDataRand::tListInputDataRand()
 **
 **  Read state from file
 **
\**************************************************************************/
tListInputDataRand::
tListInputDataRand( const tInputFile &inputfile, tRand &rand )
{
  double intime;                   // desired time
  char basename[80];

  inputfile.ReadItem( basename, sizeof(basename), "INPUTDATAFILE" );

  std::ifstream dataInfile;
  openFile( dataInfile, basename, SRANDOM);

  // Find out which time slice we want to extract
  intime = inputfile.ReadItem( intime, "INPUTTIME" );
  if (1) //DEBUG
    std::cout << "intime = " << intime << std::endl;

  // Find specified input times in input data files and read # items.
  int nn;
  findRightTime( dataInfile, nn, intime,
		 basename, SRANDOM, "random number generator");
  if ( rand.numberRecords() != nn ) {
    std::cerr << "Invalid number of records for the random number generator\n";
    ReportFatalError( "Input error" );
  }

  // Read in data from file
  rand.readFromFile( dataInfile );
}

/**************************************************************************\
 **
 **  tListInputDataVegetation::tListInputDataVegetation()
 **
 **  Read state from file
 **
\**************************************************************************/
tListInputDataVegetation::
tListInputDataVegetation( const tInputFile &inputfile )
{
  double intime;                   // desired time
  char basename[80];

  // Read base name for triangulation files from inputfile
  inputfile.ReadItem( basename, sizeof(basename), "INPUTDATAFILE" );

  std::ifstream dataInfile;
  openFile( dataInfile, basename, SVEG);

  // Find out which time slice we want to extract
  intime = inputfile.ReadItem( intime, "INPUTTIME" );
  if (1) //DEBUG
    std::cout << "intime = " << intime << std::endl;

  // Find specified input times in input data files and read # items.
  int nn;
  findRightTime( dataInfile, nn, intime,
		 basename, SVEG, "vegetation mask");
  // Read in data from file
  vegCov.setSize(nn);
  for( int i=0; i<nn; ++i ){
    dataInfile >> vegCov[i];
    if (dataInfile.fail())
      ReportIOError(IORecord, basename, SVEG, i);
  }
}
