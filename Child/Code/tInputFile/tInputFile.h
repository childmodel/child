//-*-c++-*- 

/****************************************************************************\
**
**  tInputFile.h: Header file for class tInputFile.
**
**  The tInputFile class manages an input file for single-valued parameters
**  (that is, physical constants and other parameters, but not large
**  arrays of data). The services of a tInputFile include opening a file
**  of a given name and reading items from that file.
**
**  An example of the input file format is shown below:
**
**    # This is a comment line
**    RHO: fluid density (kg/m3)
**    1000.0
**    GRAV: gravitational acceleration (m2/s)
**    9.81
**    ...etc.
**
**  Each parameter is preceded by a descriptive line of text that starts
**  with an all-caps "item code" (e.g. "RHO") which identifies the parameter.
**  This makes it possible to list parameters in any order in the input
**  file. Unused parameters and comments are simply ignored.
**
**  The overloaded ReadItem functions read and return parameters with a
**  specified item code. Possible data types include integer, double, long,
**  and string. A record of all items read is written to file
**  <filename>.inputs, where <filename> is specified by the parameter
**  having the item code OUTFILENAME.
**
**  The ReadItem function terminates with an error message if a given
**  item code can't be found in the input file. However, a version of
**  ReadItem could be added that simply returns a null value when a given
**  item does not exist.
**
**  Created by Greg Tucker, November 1997
**
**  $Id: tInputFile.h,v 1.9 2002-11-07 16:58:54 childcvs Exp $
\****************************************************************************/

#ifndef TINPUTFILE_H
#define TINPUTFILE_H

#if !defined(HAVE_NO_NAMESPACE)
# include <fstream>
using namespace std;
#else
# include <fstream.h>
#endif

#define kMaxNameLength 120
#define kCommentMark '#'

class tInputFile
{
public:
    tInputFile( const char * );   // constructor takes name of file to open
    int ReadItem( const int &, const char * );       // reads an int
    long ReadItem( const long &, const char * );     // reads a long
    double ReadItem( const double &, const char * ); // reads a double
    void ReadItem( char *, const char * );           // reads a string
    // similar overrides could be added for other data types

private:
    ifstream infile;     // the input file
    ofstream inoutfile;  // output file in which items are recorded
};

#endif
