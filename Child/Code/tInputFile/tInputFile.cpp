/****************************************************************************\
**
**  tInputFile.cpp: Member functions for class tInputFile.
**
**  (see tInputFile.h for a description of this class)
**
**  Greg Tucker, November 1997
**
**  $Id: tInputFile.cpp,v 1.17 2002-11-07 17:42:25 childcvs Exp $
\****************************************************************************/

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include "../tAssert.h"
#include <string.h>
#include <stdlib.h>
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
**  Modifications:
**    - 2/02: error check for inoutfile added (GT)
**
\****************************************************************************/
tInputFile::tInputFile( const char *filename )
  :
  infile(), inoutfile(), inoutfile_opened(false)
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
   if( !inoutfile.good() )
     {
       cerr << "Unable to open '" << inoutname << "'.\n";
       cerr << "(Error generated in module tInputFile, function tInputFile::tInputFile( const char * ) )\n";
       ReportFatalError( "The specified path name may not exist.\n" );
     }
   inoutfile_opened = true;
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
// read a line in headerLine and discards the any remaining characters
static
void readLine(char *headerLine, ifstream& infile){
  infile.getline( headerLine, kMaxNameLength );
  // still some characters left on the current line ?
  if ( !infile.eof() && !infile.bad() && infile.rdstate() & ios::failbit){
    // clear failbit
    infile.clear(infile.rdstate() & ~ios::failbit);
    // discard characters.
    char c;
    while( infile.get(c) && c != '\n');
  }
}

static
void readUntilKeyword(char *headerLine, ifstream& infile,
		      const char *itemCode){
  do {
    readLine(headerLine, infile);
  } while( !infile.eof() && !infile.bad() &&
	   ( headerLine[0]==kCommentMark ||
	     strncmp( itemCode, headerLine, strlen( itemCode ) )!=0 ) );
}

static
void skipCommentsAndReadValue(char *headerLine, ifstream& infile){
  do {
    readLine(headerLine, infile);
  } while( !infile.eof() && !infile.bad() &&
	   ( headerLine[0]==kCommentMark ) );
}


int tInputFile::ReadItem( const int & /*datType*/, const char *itemCode )
{
   //cout << "ReadItem( int )...";
   int item;
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );

   // look for itemCode
   readUntilKeyword(headerLine, infile, itemCode);
   if ( infile.eof() || infile.bad() )
     goto fail;

   // skip any comment and read value
   skipCommentsAndReadValue(headerLine, infile);
   if( headerLine[0]!=kCommentMark )
     {
       item = atoi(headerLine);
       if (inoutfile_opened)
	 inoutfile << itemCode << endl << item << endl;
     }
   else
     {
       goto fail;
     }
   infile.seekg( 0, ios::beg );
   return item;
   
 fail:
   cerr << "I expected to read the parameter '" << itemCode
	<< "', but reached EOF first" << endl;
   ReportFatalError( "Missing parameter in input file" );
}

long tInputFile::ReadItem( const long & /*datType*/, const char *itemCode )
{
   //cout << "ReadItem( long )...";
   long item;
   char headerLine[kMaxNameLength];
  
   assert( infile.good() );
  
   // look for itemCode
   readUntilKeyword(headerLine, infile, itemCode);
   if ( infile.eof() || infile.bad() )
     goto fail;

   // skip any comment and read value
   skipCommentsAndReadValue(headerLine, infile);
   if( headerLine[0]!=kCommentMark )
     {
       item = atol(headerLine);
       if (inoutfile_opened)
	 inoutfile << itemCode << endl << item << endl;
     }
   else
     {
       goto fail;
     }
   infile.seekg( 0, ios::beg );
   return item;

 fail:
   cerr << "I expected to read the parameter '" << itemCode
	<< "', but reached EOF first" << endl;
   ReportFatalError( "Missing parameter in input file" );
}

double tInputFile::ReadItem( const double & /*datType*/, const char *itemCode )
{
   //cout << "ReadItem( double )...";
   double item;
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );
  
   // look for itemCode
   readUntilKeyword(headerLine, infile, itemCode);
   if ( infile.eof() || infile.bad() )
     goto fail;

   // skip any comment and read value
   skipCommentsAndReadValue(headerLine, infile);
   if( headerLine[0]!=kCommentMark )
     {
       item = atof(headerLine);
       if (inoutfile_opened)
	 inoutfile << itemCode << endl << item << endl;
     }
   else
     {
       goto fail;
     }
   infile.seekg( 0, ios::beg );
   return item;

 fail:
   cerr << "I expected to read the parameter '" << itemCode
	<< "', but reached EOF first" << endl;
   ReportFatalError( "Missing parameter in input file" );
}


void tInputFile::ReadItem(  char * theString, const char *itemCode )
{
   //cout << "ReadItem( char )...";
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );

   // look for itemCode
   readUntilKeyword(headerLine, infile, itemCode);
   if ( infile.eof() || infile.bad() )
     goto fail;

   // skip any comment and read value
   do {
     readLine(headerLine, infile);
   } while( !infile.eof() && !infile.bad() &&
	    ( headerLine[0]==kCommentMark ) );
   if( headerLine[0]!=kCommentMark )
     {
       strcpy(theString,headerLine);
       if (inoutfile_opened)
	 inoutfile << itemCode << endl << theString << endl;
     }
   else
     {
       goto fail;
     }
   infile.seekg( 0, ios::beg );
   return;

 fail:
   cerr << "I expected to read the parameter '" << itemCode
	<< "', but reached EOF first" << endl;
   ReportFatalError( "Missing parameter in input file" );
}
