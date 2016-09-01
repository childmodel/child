//-*-c++-*- 

/****************************************************************************/
/**
**  @file tInputFile.h
**  @brief Header file for class tInputFile.
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
**  &lt;filename&gt;.inputs, where &lt;filename&gt; is specified by the parameter
**  having the item code OUTFILENAME.
**
**  The ReadItem function terminates with an error message if a given
**  item code can't be found in the input file. However, a version of
**  ReadItem could be added that simply returns a null value when a given
**  item does not exist.
**
**  Created by Greg Tucker, November 1997
**  Re-written, AD, July 2003
**  Modified, SL, September 2003: Can now generate non-fatal warnings and
**    return zero-values as specified; default will be fatal error.
**
**  $Id: tInputFile.h,v 1.26 2004/05/27 17:20:55 childcvs Exp $
*/
/****************************************************************************/

#ifndef TINPUTFILE_H
#define TINPUTFILE_H

#include "../tArray/tArray.h"

#define kMaxNameLength 128
#include "../Definitions.h"
#include <stddef.h>
class tTimeSeries;

/** @class tKeyPair
**
**  Hold a key/value pair
*/
class tKeyPair
{
public:
  tKeyPair() : key_(0), value_(0) {}
  tKeyPair(int) : key_(0), value_(0) {} // to make tArray happy
  tKeyPair(const char *, const char *);
  tKeyPair(tKeyPair const&);
  tKeyPair& operator=(tKeyPair const&);
  ~tKeyPair();
  void setKey(const char *);
  void setValue(const char *);
  void clear();
  const char *key() const {return key_;}
  const char *value() const {return value_;}
private:
  void clearKey();
  void clearValue();

  char *key_;
  char *value_;
};

/** @class tInputFile
*/
class tInputFile
{
public:
  tInputFile( const char * );   // constructor takes name of file to open
  bool Contain(const char *) const;
  void WarnObsoleteKeyword(const char *, const char *) const;

  bool ReadBool( const char *, bool reqParam = true  ) const;                  // reads a bool
  int ReadItem( const int &, const char *, bool reqParam = true ) const;       // reads an int
  int ReadInt( const char *, bool reqParam = true  ) const;		       // reads an int (alt)
  long ReadItem( const long &, const char *, bool reqParam = true ) const;     // reads a long
  long ReadLong( const char *, bool reqParam = true ) const;                   // reads a long (alt)
  double ReadItem( const double &, const char *, bool reqParam = true ) const; // reads a double
  double ReadDouble( const char *, bool reqParam = true ) const;               // reads a double (alt)
  void ReadItem( char *, size_t len, const char *, bool reqParam = true ) const;// reads a string
  void ReadItem(tTimeSeries &, const char *,
		bool reqParam = true ) const; // reads a time series

  // similar overrides could be added for other data types

  tArray< tKeyPair > & GetKeyWordTableRef();  // Returns a reference to the keyword table

private:
  tArray< tKeyPair > KeyWordTable; // hold key/value pair
  enum { notFound = -1 }; // must be strictly negative
  int findKeyWord(const char*) const; // find index of keyword
  void writeLogFile() const;

  tInputFile(tInputFile const&);
  tInputFile& operator=(tInputFile const&);
  tInputFile();
};

#endif
