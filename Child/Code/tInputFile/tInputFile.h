/****************************************************************************\
**
**  tInputFile.h: Header file for tInput objects.
**
**  Greg Tucker, November 1997
**
**  $Id: tInputFile.h,v 1.3 1999-02-01 21:45:22 gtucker Exp $
\****************************************************************************/

#ifndef TINPUTFILE_H
#define TINPUTFILE_H

#include <iostream.h>
#include <fstream.h>

#define kMaxNameLength 80
#define kCommentMark '#'

class tInputFile
{
public:
    tInputFile();
    tInputFile( const char * );
    int ReadItem( const int &, const char * );
    long ReadItem( const long &, const char * );
    double ReadItem( const double &, const char * );
    void ReadItem( char *, const char * );
    // similar overrides for strings, chars, etc

private:
    ifstream infile;
};

#endif
