/****************************************************************************\
**
**  tInputFile.h: Header file for tInput objects.
**
**  Greg Tucker, November 1997
**
**  $Id: tInputFile.h,v 1.2 1998-01-29 20:15:32 stlancas Exp $
\****************************************************************************/

#include <iostream.h>
#include <fstream.h>

#ifndef TINPUTFILE_H
#define TINPUTFILE_H
/*#define kMaxNameLength 80
#define kCommentMark '#'*/

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
