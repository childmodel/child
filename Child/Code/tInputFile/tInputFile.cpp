/****************************************************************************\
**
**  tInputFile.cpp: Member functions for tInput objects.
**
**  Greg Tucker, November 1997
**
**  $Id: tInputFile.cpp,v 1.2 1998-01-27 22:48:03 gtucker Exp $
\****************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <string.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "tInputFile.h"
#include "../errors/errors.h"

/*
**  Constructors
*/

tInputFile::tInputFile( const char *filename )
{
   infile.open( filename );
   if( !infile.good() )
       cerr << "tInputFile::tInputFile: Unable to open '" << filename << "'." << endl;
}


/*
**  ReadItem
**
**  Reads one parameter from the file. The format is assumed to be a line
**  of text that begins with the code "code", followed by a line containing
**  the parameter to be read. The function is overloaded according to the
**  type of data desired (datType simply governs which overloaded function
**  will be called; it is not used by the routines).
**
**  Notes: add ability to ignore line followed by a comment mark,
**
**  IMPORTANT:
**  revised to allow arbitrary ordering of items in infile and/or ReadItem
**  calls in code; routine searches
**  through list until it either finds the right itemCode or reaches EOF.
**  This should make reading/entering parameters MUCH less complicated.
**  -12/23/97 SL
*/
//template <class T>  how do w/ templates?
int tInputFile::ReadItem( const int &datType, const char *itemCode )
{
   cout << "ReadItem( int )...";
   int item;
   char headerLine[kMaxNameLength];
   
   assert( infile.good() );

   //streampos original = infile.tellg();
   
   
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
      cout << headerLine << endl;
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
   cout << "ReadItem( long )...";
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
      cout << headerLine << endl;
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

float tInputFile::ReadItem( const float &datType, const char *itemCode )
{
   cout << "ReadItem( float )...";
   float item;
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
      cout << headerLine << endl;
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
   cout << "ReadItem( char )...";
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
      cout << headerLine << endl;
   }
   else
   {
      cerr << "I expected to read the parameter '" << itemCode
           << "', but reached EOF first" << endl;
      ReportFatalError( "Missing parameter in input file" );
   }
   infile.seekg( 0, ios::beg );
}
