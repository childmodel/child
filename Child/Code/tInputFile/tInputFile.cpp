/****************************************************************************\
**
**  tInputFile.cpp: Member functions for class tInputFile.
**
**  (see tInputFile.h for a description of this class)
**
**  Greg Tucker, November 1997
**
**  $Id: tInputFile.cpp,v 1.6 1999-03-31 21:24:04 gtucker Exp $
\****************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <string.h>
//#include "../Definitions.h"
//#include "../Classes.h"
#include "tInputFile.h"
#include "../errors/errors.h"


/****************************************************************************\
**
**  tInputFile Constructor
**
**  Looks for a file called filename, opens it if found or generates an
**  error if not. Then reads the base name for output files and creates
**  a file called <filename>.inputs which will act as a log file for
**  parameters read. (This is often useful when the original input file
**  gets lost or modified).
**
\****************************************************************************/
tInputFile::tInputFile( const char *filename )
{
   char inoutname[kMaxNameLength];

   // Open file
   infile.open( filename );
   if( !infile.good() )
   {
      cerr << "tInputFile::tInputFile: Unable to open '" << filename << "'." << endl;
      ReportFatalError( "The file may not exist or may be mis-named." );
   }

   // Create log file for inputs
   ReadItem( inoutname, "OUTFILENAME" );
   strcat( inoutname, ".inputs" );
   inoutfile.open( inoutname );
   assert( inoutfile.good() );
   
}


/****************************************************************************\
**
**  tInputFile::ReadItem
**
**  Reads one parameter from the file. The format is assumed to be a line
**  of text that begins with the code "itemCode", followed by a line containing
**  the parameter to be read. The function is overloaded according to the
**  type of data desired (datType simply governs which overloaded function
**  will be called; it is not used by the routines).
**
**  Inputs:  datType -- dummy variable indicating the data type to be read
**                      (in the case of the string version, the string read
**                      is placed here)
**           itemCode -- string that describes the parameter to be read
**  Returns:  the item read (except in the case of the string version)
**  Modifications:
**    - revised to allow arbitrary ordering of items in infile and/or
**      ReadItem calls in code; routine searches through
**      list until it either finds the right itemCode or reaches EOF.
**      12/23/97 SL
**
\****************************************************************************/
int tInputFile::ReadItem( const int &datType, const char *itemCode )
{
   //cout << "ReadItem( int )...";
   int item;
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );

   // NB: Should check for eof on reading each line

   infile.getline( headerLine, kMaxNameLength );
   while( !( infile.eof() ) &&
          ( headerLine[0]==kCommentMark ||
            strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
       infile.getline( headerLine, kMaxNameLength );
   if( !( infile.eof() ) )
   {
      infile >> item;
      infile.ignore( 1, '\n' );
      inoutfile << itemCode << endl << item << endl;
   }
   else

       //if( strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 )
   {
      cerr << "I expected to read the parameter '" << itemCode
           << "', but reached EOF first" << endl;
      ReportFatalError( "Missing parameter in input file" );
   }
   //infile.seekg( original );
   infile.seekg( 0, ios::beg );
   return item;
}

long tInputFile::ReadItem( const long &datType, const char *itemCode )
{
   //cout << "ReadItem( long )...";
   long item;
   char headerLine[kMaxNameLength];
  
   assert( infile.good() );
  
     // NB: Should check for eof on reading each line

   infile.getline( headerLine, kMaxNameLength );
   while( !( infile.eof() ) &&
          ( headerLine[0]==kCommentMark ||
            strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
       infile.getline( headerLine, kMaxNameLength );
   if( !( infile.eof() ) )
   {
      infile >> item;
      infile.ignore( 1, '\n' );
      inoutfile << itemCode << endl << item << endl;
   }
   else
   {
      cerr << "I expected to read the parameter '" << itemCode
           << "', but reached EOF first" << endl;
      ReportFatalError( "Missing parameter in input file" );
   }
   infile.seekg( 0, ios::beg );
   return item;
}

double tInputFile::ReadItem( const double &datType, const char *itemCode )
{
   //cout << "ReadItem( double )...";
   double item;
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );
  
     // NB: Should check for eof on reading each line

   infile.getline( headerLine, kMaxNameLength );
   while( !( infile.eof() ) &&
          ( headerLine[0]==kCommentMark ||
            strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
       infile.getline( headerLine, kMaxNameLength );
   if( !( infile.eof() ) )
   {
      infile >> item;
      infile.ignore( 1, '\n' );
      inoutfile << itemCode << endl << item << endl;
   }
   else
   {
      cerr << "I expected to read the parameter '" << itemCode
           << "', but reached EOF first" << endl;
      ReportFatalError( "Missing parameter in input file" );
   }
   infile.seekg( 0, ios::beg );
   return item;
}


void tInputFile::ReadItem(  char * theString, const char *itemCode )
{
   //cout << "ReadItem( char )...";
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );

   infile.getline( headerLine, kMaxNameLength );
   while( !( infile.eof() ) &&
          ( headerLine[0]==kCommentMark ||
            strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) )
       infile.getline( headerLine, kMaxNameLength );

   if( !( infile.eof() ) )
   {
      infile.getline( theString, kMaxNameLength );
      inoutfile << itemCode << endl << theString << endl;
   }
   else
   {
      cerr << "I expected to read the parameter '" << itemCode
           << "', but reached EOF first" << endl;
      ReportFatalError( "Missing parameter in input file" );
   }
   infile.seekg( 0, ios::beg );
}
