//-*-c++-*-

//--------------------------------------------------------------
// tTimeSeries.cpp
//
// This file implements temporal variable parameters into CHILD.
// See the file TimeSeries.pdf for info how to use it.
//
// Copyright 2000 P. W. Bogaart / VUA
// January 2004: Cleaned up and integrated in Child (AD)
// Oct. 2010: Gave tTimeSeries a copy constructor and the various
//  sub-classes of tTimeSeriesImp virtual Initialize_Copy f'ns
//  so that different sub-classes are instantiated in the copy. (SL)
//--------------------------------------------------------------

#include "tTimeSeries.h"

#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>

#include "../Mathutil/mathutil.h"
#include "../errors/errors.h"

#include <stdio.h>
#include <iostream>

//--------------------------------------------------------------
// tToklist is a simple class which implements string-to-tokens
// conversion.
// It is used internally by the other classes in this file,
// but you are free to use it for other purposes.
//--------------------------------------------------------------

class tTokList {
 public:
  tTokList();
  ~tTokList();

  int parse(const char *s);
  int parse(const char *s, const char *ws);

  const char *operator[](const int i) const;
  const char *token(const int i) const;

  int size() const;
 private:
  int ntok;
  char *tok[64];
};


//--------------------------------------------------------------
// tDataFile is used for the parsing of data files.
// It is used internally by the other classes in this file,
// but you are free to use it for other purposes.
//--------------------------------------------------------------

#define DF_MAXCOLSIZE 1024
#define DF_MAXCOL 32

class tDataFile {
 public:
  tDataFile();
  ~tDataFile();
  int     cols() const { return ncol; }
  int     rows() const { return nrow; }
  void    load(const char *filename);
  void    print(const char *filename);
  void    print(FILE *fp);
  bool    hasColumn(const char *name) const;
  int     getColIndex(const char *name) const;
  double  *getColumn(int colindex) const;
  double  *getColumn(const char *name) const;
  const char   *getParam(const char *name) const;
 private:
  int     nrow, ncol;             // number of rows, columns in data table
  double  *data[DF_MAXCOL];        // list of ptrs to column data
  char   info[256];               // info string for datafile
  int    ncolname;                // number of column names in dtafile
  int    ncomment;                // number of comments in datafile
  int    nparam;                  // number of parameters in datafile
  char   colname[DF_MAXCOL][32];  // list of column names
  char   comment[DF_MAXCOL][128]; // list of comments
  char   parname[DF_MAXCOL][32];  // list of parameter names
  char   parvalue[DF_MAXCOL][32]; // list of parameter values
};


//--------------------------------------------------------------
// tTimeSeriesImp is the base of the class hierarchy.
// It's main purpose it to provide an interface to the
// virtual member functions 'configure' and 'calc'.
//
// The current hierarchy is:
//   tTimeSeriesImp
//      tConstantTimeSeriesImp
//      tLinearGrowthTimeSeriesImp
//      tWaveTimeSeriesImp
//         tSinwaveTimeSeriesImp
//         tBlockwaveTimeSeriesImp
//      tDataTimeSeriesImp
//         tFileTimeSeriesImp
//         tInlineTimeSeriesImp
//--------------------------------------------------------------

class tTimeSeriesImp {
 public:
  virtual void configure(const char *s) = 0;
  virtual double calc(double time) const = 0;
  virtual void TellAll() const = 0;
  virtual ~tTimeSeriesImp() {}
  virtual void Initialize_Copy( tTimeSeriesImp* ) =0;
};


//--------------------------------------------------------------
// tConstantTimeSeriesImp is used for parameter defs in the form of
//     @constant 3.14
//--------------------------------------------------------------

class  tConstantTimeSeriesImp : public tTimeSeriesImp {
 protected:
  double value;
 public:
  virtual void configure(const char *s);
  virtual double calc(double time) const;
  virtual void TellAll() const;
  virtual void Initialize_Copy( tTimeSeriesImp* );
};

inline void tConstantTimeSeriesImp::Initialize_Copy(tTimeSeriesImp* oPtr)
{
  tConstantTimeSeriesImp* ptr=static_cast<tConstantTimeSeriesImp*>(oPtr);
  value = ptr->value;
}


//--------------------------------------------------------------
// tLinearGrowthTimeSeriesImp is used for parameter defs in the form of
//     @linear 0. 3.1
//--------------------------------------------------------------

class tLinearGrowthTimeSeriesImp : public tTimeSeriesImp {
 protected:
  double initialValue, rate;
 public:
  virtual void configure(const char *s);
  virtual double calc(double time) const;
  virtual void TellAll() const;
  virtual void Initialize_Copy( tTimeSeriesImp* );
};

inline void tLinearGrowthTimeSeriesImp::Initialize_Copy(tTimeSeriesImp* oPtr)
{
  tLinearGrowthTimeSeriesImp* ptr=static_cast<tLinearGrowthTimeSeriesImp*>(oPtr);
  initialValue = ptr->initialValue;
  rate = ptr->rate;
}

//--------------------------------------------------------------
// tWaveTimeSeriesImp is used as a common ancestor for the
// specializes sinwave and blockwave classes. You are free to
// create your own variant.
//--------------------------------------------------------------

class tWaveTimeSeriesImp : public tTimeSeriesImp {
 protected:
  double mean, amp, period, lag;
 public:
  virtual void configure(const char *s);
  virtual void Initialize_Copy( tTimeSeriesImp* );
};

inline void tWaveTimeSeriesImp::Initialize_Copy(tTimeSeriesImp* oPtr)
{
  tWaveTimeSeriesImp* ptr=static_cast<tWaveTimeSeriesImp*>(oPtr);
  mean = ptr->mean;
  amp = ptr->amp;
  period = ptr->period;
  lag = ptr->lag;
}

//--------------------------------------------------------------
// tSinwaveTimeSeriesImp is used for parameter defs of the form
//    @sinwave <mean> <amplitude> <period> <lag>
//--------------------------------------------------------------

class tSinwaveTimeSeriesImp : public tWaveTimeSeriesImp {
 public:
  virtual double calc(double time) const;
  virtual void TellAll() const;
};

//--------------------------------------------------------------
// tBlockWaveTimeSeriesImp is used for parameter defs of the form
//    @blockwave <mean> <amplitude> <period> <lag>
//--------------------------------------------------------------

class tBlockwaveTimeSeriesImp : public tWaveTimeSeriesImp {
 public:
  virtual double calc(double time) const;
  virtual void TellAll() const;
};

//--------------------------------------------------------------
// tDataTimeSeriesImp is used as a common ancestor for all classes
// who handle parameter sets in the form of time-value tables
// Currently there are tFileTimeSeriesImp and tInlineTimeSeriesImp,
// but you are free to create you own
//--------------------------------------------------------------

class tDataTimeSeriesImp : public tTimeSeriesImp {
 public:
  tDataTimeSeriesImp();
  ~tDataTimeSeriesImp();
  virtual double calc(double time) const;
  virtual void Initialize_Copy( tTimeSeriesImp* );
 protected:
  double *xdata, *ydata;  // here's the data
  int    n;               // number of data pairs
  bool   interpolate;     // interpolate/forward mode switch
};

inline void tDataTimeSeriesImp::Initialize_Copy(tTimeSeriesImp* oPtr)
{
  tDataTimeSeriesImp* ptr=static_cast<tDataTimeSeriesImp*>(oPtr);
  n = ptr->n;
  interpolate = ptr->interpolate;
  xdata = new double[n];
  ydata = new double[n];
  for( int i=0; i<n; ++i )
    {
      xdata[i] = ptr->xdata[i];
      ydata[i] = ptr->ydata[i];
    }
}

//--------------------------------------------------------------
// tFileTimeSeriesImp is used for paramter defs of the form
//     @file: <filename> <time col> <val col> <mode>
//--------------------------------------------------------------

class tFileTimeSeriesImp : public tDataTimeSeriesImp {
 public:
  virtual void configure(const char *s);
  virtual void TellAll() const;
};


//--------------------------------------------------------------
// tInlineTimeSeriesImp is used for paramter defs of the form
//     @inline: <time1>:<val1> <time2><val2> ... <mode>
//--------------------------------------------------------------

class tInlineTimeSeriesImp : public tDataTimeSeriesImp {
 public:
  virtual void configure(const char *s);
  virtual void TellAll() const;
};


/*********************************************************************\
 * tTokList class                                                     *
 *                                                                    *
 * This class is used to parse a string into tokens, and              *
 * store them. The tokens can accessed via the [] operator.           *
 *                                                                    *
 * feb 11 2000, PWB                                                   *
 *                                                                    *
\********************************************************************/

tTokList::tTokList() :
  ntok(0)
{}

tTokList::~tTokList()
{
  int i;

  // delete all tokens, if any
  for (i=0; i<ntok; i++)
    delete[] tok[i];
  ntok = 0;
}

int tTokList::parse(const char *s)
{
  // parse string 's' into tokens
  // use default whitespace characters
  const char * const default_ws = " \f\n\r\t\v";
  return parse(s, default_ws);
}

int tTokList::parse(const char *s, const char *ws)
{
  // parse string 's' into tokens, use whitespace
  // set 'ws'
  int   i;
  char *buffer, *p;

  // delete any existing tokens
  for (i=0; i<ntok; i++)
    delete[] tok[i];
  ntok = 0;

  // parse string with help of strtok() lib function
  // and insert tokens in list
  buffer = new char[strlen(s)+1];
  strcpy(buffer,s);
  p = strtok(buffer, ws);
  while (p) {
    tok[ntok] = new char[strlen(p)+1];
    strcpy(tok[ntok], p);
    ntok++;
    p = strtok(0, ws);
  }
  delete[] buffer;

  return ntok;
}

const char *tTokList::operator[](const int i) const
{
  // return token at index 'i'
  assert(i>=0);
  assert(i<ntok);
  return tok[i];
}

const char *tTokList::token(const int i) const
{
  //same as operator[], but functional form
  assert(i>=0);
  assert(i<ntok);
  return tok[i];
}

int tTokList::size() const
{
  // return number of tokens in list
  return ntok;
}


/*********************************************************************\
 * tDataFile class                                                     *
 *                                                                     *
 * This class is used to read datafiles consisting of tabular data.    *
 * A number of extra items can be stored in these files, such as       *
 * an info-string, a list of column names, and name=value pairs        *
 *                                                                     *
 * Example file:                                                       *
 *                                                                     *
 * %info: This is an example of a data-file.                           *
 * %This is a comment                                                  *
 * %                                                                   *
 * %Now a name=value parameter:                                        *
 * %pi=3.1415                                                          *
 * %                                                                   *
 * %The column names are stored like this:                             *
 * %columns: time  discharge sediment_supply                           *
 * %units:   yr    m3/s      kg/s                                      *
 * %That makes 3 columns. Note the units definition. this is not       *
 * %required                                                           *
 * %                                                                   *
 * %now the raw data:                                                  *
 * 1970  2145.5  500.1                                                 *
 * 1971  2453.2  531.7                                                 *
 * 1972  1400.6  481.9                                                 *
 * .     .       .                                                     *
 * .     .       .                                                     *
 * .     .       .                                                     *
 *                                                                     *
 *                                                                     *
 *                                                                     *
\*********************************************************************/


/*********************************************************************\
 * char *strclean(const char *s)                                       *
 *                                                                     *
 * This function cleans string 's'. leading and trailing whitespace is *
 * removed. consecutive whitespace characters are being replaced by    *
 * a single space. whitespace charcacters are: space, tab, return,     *
 * newline, vertical tab, newpage.                                     *
 * The resulting, cleansed string, is returned.                        *
 * Examples:                                                           *
 *   "abcd"              -> "abcd"                                     *
 *   "   abc   def   "   -> "abc def"                                  *
 *                                                                     *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/
static
char *strclean(const char *s)
{
  static char s2[1024];     // new cleansed, string
  tTokList  toklist;        // token list
  int       i;              // token interator

  toklist.parse(s);
  strcpy(s2,"");
  for (i=0; i<toklist.size(); i++)
    {
      if (i>0) strcat(s2," ");
      strcat(s2, toklist[i]);
    }

  return s2;
}


/*********************************************************************\
 * tDataFile::tDataFile()                                              *
 *                                                                     *
 * Default constructor.                                                *
 * Does nothing except setting field to 0                              *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/

tDataFile::tDataFile() :
  nrow(0), ncol(0),
  ncolname(0),
  ncomment(0),
  nparam(0)
{
  info[0] = '\0';
  for(size_t i=0;i<sizeof(data)/sizeof(data[0]);++i)
    data[i] = NULL;
}


/*********************************************************************\
 * tDataFile::~tDataFile()                                             *
 *                                                                     *
 * Destructor.                                                         *
 * Frees dynamically allocated memory                                  *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/

tDataFile::~tDataFile()
{
  for(int i=0; i<ncol; i++) {
    delete[] data[i];
  }
}

/*********************************************************************\
 * void tDataFile::load(char *filename)                                *
 *                                                                     *
 * Loads the contents of file 'filename' into the objects's database.  *
 * See the object's description for the file format                    *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/


void tDataFile::load(const char *filename)
{
  FILE *fp;                     // pointer to the actual file
  char *p, *pp;                 // pointers, used for string-traversing
  static char line[1024];       // buffer to hold a single line
  static char s[1024];          // buffer to hold copy of current line, for destructive processes
  int c;                        // column-iterator
  int check;                    // holds number of items actually read.
  tTokList toklist;             // token list for data record parsing

  bool debug = false;           // if true, then produce output about what we're doing

  // open file, and report error if failure
  fp = fopen(filename, "rt");
  if (fp == NULL) {
    std::cout << "Time series file `" << filename << "' does not exist."
	      << std::endl;
    ReportFatalError("File not found.");
  }
  assert(fp);

  nrow=0;

  while (!feof(fp)) {     // iterate over all lines

    // read line
    fgets(line, 1024, fp);
    if (debug) printf("LINE: '%s'", line);

    // no, depending on line contents and context, perform some action:

    // case 1: line is empty; skip this line
    if (strlen(strclean(line))==0) {
      if (debug) printf("ACTION: skip empty line\n");
      continue;
    }

    // case 2: "%info" line; set info field
    if (strstr(line,"%info:")==line) {
      if (debug) printf("ACTION: parsing info field\n");
      strcpy(info, strclean(line+6));
      continue;
    }

    // case 3: %columns line; extract column names
    if (strstr(line,"%columns:")==line) {
      if (debug) printf("ACTION: parsing column def\n");
      strcpy(s,strclean(line+9));
      ncolname=0;
      p = strtok(s," \t,;");
      while (p) {
	strcpy(colname[ncolname++], p);
	p = strtok(0," \t,;");
      }
      continue;
    }

    // case 4: %name=value line, extract name,value pair
    if (strstr(line,"%")==line && strstr(line,"=")!=0) {
      if (debug) printf("ACTION: parsing name=val def\n");
      strcpy(s, line+1);  // line without '#'; 'name=value'
      p = strstr(s,"=");
      pp = p+1;           // points to 'value'
      *p = '\0';          // s points now to 'name'
      strcpy(parname[nparam], strclean(s));
      strcpy(parvalue[nparam], strclean(pp));
      nparam++;
      continue;
    }

    // case 5: %comment line
    if (strstr(line,"%")==line) {
      if (debug) printf("ACTION: parsing comment\n");
      strcpy(comment[ncomment++], strclean(line+1));
      continue;
    }

    // case 6: first data row; determine number of columns
    //         and store data
    if (nrow==0) {
      if (debug) printf("ACTION: parsing first data row\n");
      toklist.parse(line);
      ncol = toklist.size();

      // allocate memory to hold columns (size=DF_MAXCOLSIZE)
      for (c=0; c<ncol; c++) {
	data[c] = new double[DF_MAXCOLSIZE];
      }

      // store this first row of data
      for (c=0; c<ncol; c++) {
	check = sscanf(toklist[c], "%lf", &data[c][nrow]);
	assert(check==1);
      }

      nrow++;

      continue;
    }

    // case 7: not-first data row; read row and store it
    if (nrow!=0) {
      toklist.parse(line);
      if (debug) printf("ACTION: parsing data row\n");
      check=0;
      for (c=0; c<ncol; c++) {
	check += sscanf(toklist[c], "%lf", &data[c][nrow]);
	if (debug) printf("DATA %.3f", data[c][nrow]);
      }
      if (debug) printf("\n");
      assert(check==ncol);
      nrow++;
      continue;
    }

    // case 8: should not happen
    assert(0);
  }

  fclose(fp);
  return;
}

/*********************************************************************\
 * tDataFile::print(char *filename)                                   *
 *                                                                    *
 * Writes the contents of the datafile, as parsed, to a new file      *
 * named 'filename'. if 'filename' is NULL, output is  written to     *
 * stdout.                                                            *
 *                                                                    *
 * 26/01/2000, pwb                                                    *
 *                                                                    *
\*********************************************************************/

void tDataFile::print(const char *filename)
{
  FILE *fp;

  if (filename) {
    fp = fopen(filename,"wt");
    print(fp);
    fclose(fp);
  } else {
    print(stdout);
  }
}

/*********************************************************************\
 * tDataFile::print(FILE *fp)                                         *
 *                                                                    *
 * Writes the contents of the datafile, as parsed, to                 *
 * the file pointed to by `fp'.                                       *
 *                                                                    *
 * 26/01/2000, pwb                                                    *
 *                                                                    *
\*********************************************************************/


void tDataFile::print(FILE *fp)
{
  int c,r,i;

  if (strlen(info)>0)
    fprintf(fp,"%%info: %s\n", info);

  if (ncolname>0) {
    fprintf(fp,"%%columns:");
    for (c=0; c<ncolname; c++) {
      fprintf(fp," %s", colname[c]);
    }
    fprintf(fp,"\n");
  }

  if (ncomment>0) {
    for (i=0; i<ncomment; i++)
      fprintf(fp,"%%%s\n", comment[i]);
  }

  if (nparam>0) {
    for (i=0; i<nparam; i++)
      fprintf(fp,"#%s=%s\n", parname[i], parvalue[i]);
  }

  for (r=0; r<nrow; r++) {
    for (c=0; c<ncol; c++)
      fprintf(fp,"%.5f ", data[c][r]);
    fprintf(fp,"\n");
  }

  return;
}


/*********************************************************************\
 * bool tDataFile::hasColumn(char *name)                               *
 *                                                                     *
 * checks if a column 'name' exists                                    *
 * returns true if so, false otherwise                                 *
 *                                                                     *
 * 09/02/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/

bool tDataFile::hasColumn(const char *name) const
{
  int i;

  for (i=0; i<ncolname; i++)
    if (strcmp(colname[i], name)==0)
      return true;

  return false;
}

/*********************************************************************\
 * int tDataFile::getColIndex(char *name)                              *
 *                                                                     *
 * Returns the column number (zero-based) for the column named         *
 * 'name'.                                                             *
 * An index of -1 is returned if there is'n such a column.             *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/

int tDataFile::getColIndex(const char *name) const
{
  int index, i;

  index = -1;
  for (i=0; i<ncolname; i++)
    if (strcmp(colname[i], name)==0)
      index=i;

  return index;
}


/*********************************************************************\
 * double *tDataFile::getColumn(int colindex)                           *
 * double *tDataFile::getColumn(char *colname)                          *
 *                                                                     *
 * Returns a copy of the data in the column pointed to by the          *
 * parameter 'colindex' (zero-based)                                   *
 * The memory required to hold the data is newly allocated.            *
 * A NULL pointer is returned if the specified column doesn't exist    *
 *                                                                     *
 * double *tDataFile::getColumn(char *colname)                          *
 *                                                                     *
 * Idem, but uses 'colname'  indentify the column requested.           *
 * The memory required to hold the data is newly allocated.            *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/


double *tDataFile::getColumn(int colindex) const
{
  double *clone;        // ptr to copy of column data

  if (colindex<0 || colindex>=cols()) {
    return 0;
  }

  clone = new double[nrow];
  for (int row=0; row<nrow; row++) {
    clone[row] = data[colindex][row];
  }

  return clone;
}

double *tDataFile::getColumn(const char *colname_) const
{
  return getColumn( getColIndex(colname_) );
}

/*********************************************************************\
 * const char *tDataFile::getParam(const char *name) const             *
 *                                                                     *
 * Returns the value of parameter 'name', as read during the parsing
 * as a %name=value pair.                                              *
 * Returns an empty string ("") if there is no parameter 'name'        *
 *                                                                     *
 * 26/01/2000, pwb                                                     *
 *                                                                     *
\*********************************************************************/
const char *tDataFile::getParam(const char *name) const
{
  static char result[32];

  strcpy(result, "");
  for (int i=0; i<nparam; i++)
    if (strcmp(name, parname[i])==0)
      strcpy(result, parvalue[i]);
  return result;
}

/*********************************************************************\
 * nr_hunt                                                            *
 * copyright: Numerical Recipes                                       *
\*********************************************************************/
static
void nr_hunt(const double *xx, int n, double x, int& jlo)
  // given an array xx[1..n], and given a value x, returns a value jlo such that
  // x is between xx[jlo] and xx[jlo+1].
  // xx[1..n] must be monotonic, either increasing or decreasing.
  // jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input
  // is taken as the initial guess for jlo on output
{
  int jm, jhi, inc;

  const bool ascnd = (xx[n] >= xx[1]); // True if ascending order of tables

  if (jlo<=0 || jlo>n) { // Input use not usefull, go to bissection
    jlo = 0;
    jhi=n+1;
  } else {
    inc=1; // Set the hunting increment
    if ((x >= xx[jlo]) == ascnd) { // Hunt up
      if (jlo==n) return;
      jhi = jlo+1;
      while ((x >= xx[jhi]) == ascnd) { // Not done hunting
        jlo=jhi;
        inc += inc; // so double increment
        jhi=jlo+inc;
        if (jhi>n) { // Done hunting, since off end of table
          jhi=n+1;
          break;
        } // Try again
      } // Done hunting, value bracketed
    } else { // Hunt down
      if (jlo==1) {
        jlo=0;
        return;
      }
      jhi = jlo--;
      while ((x < xx[jlo]) == ascnd) { // Not done hunting
        jhi = jlo;
        inc <<= 1; // so double the increment
        if (inc >= jhi) { // Done hunting, since off end to table
          jlo=0;
          break;
        }
        else jlo = jhi-inc;
      } // and try again
    } // Done hunting, value bracketed
  } // Hunt is done, so begin the final bissection phase

  while ((jhi - jlo) != 1) {
    jm = (jhi+jlo) >> 1;
    if ((x >= xx[jm]) == ascnd)
      jlo = jm;
    else
      jhi = jm;
  }
  if (x == xx[n]) jlo=n;
  if (x == xx[1]) jlo=0;
}

/*********************************************************************\
 * tConstantTimeSeriesImp                                                *
\*********************************************************************/

void tConstantTimeSeriesImp::configure(const char *s)
{
  int check;
  // s must be a string rep. of a floating point, e.g. "3.1415"
  check = sscanf(s,"%lf", &value);
  assert(check==1);
}

double tConstantTimeSeriesImp::calc(double) const
{
  return value;
}

void tConstantTimeSeriesImp::TellAll() const
{
  std::cout << "TS @constant. value=" << value << std::endl;
}

/*********************************************************************\
 * tLinearGrowthTimeSeriesImp
\*********************************************************************/

void tLinearGrowthTimeSeriesImp::configure(const char *s)
{
  // s must be "mean amplitude period lag" format
  int check;
  check = sscanf(s,"%lf %lf", &initialValue, &rate);
  assert(check==2);
}

double tLinearGrowthTimeSeriesImp::calc(double time) const
{
  return initialValue + rate*time;
}

void tLinearGrowthTimeSeriesImp::TellAll() const
{
  std::cout
    << "TS @linear."
    << " initialValue=" << initialValue
    << " rate=" << rate
    << std::endl;
}
/*********************************************************************\
 * tWaveTimeSeriesImp                                                    *
\*********************************************************************/

void tWaveTimeSeriesImp::configure(const char *s)
{
  // s must be "mean amplitude period lag" format
  int check;
  check = sscanf(s,"%lf %lf %lf %lf", &mean, &amp, &period, &lag);
  assert(check==4);
}

/*********************************************************************\
 * tSinwaveTimeSeriesImp                                                    *
\*********************************************************************/

double tSinwaveTimeSeriesImp::calc(double time) const
{
  return mean + amp*sin((time-lag)*TWOPI/period);
}

void tSinwaveTimeSeriesImp::TellAll() const
{
  std::cout
    << "TS @sinwave."
    << " mean=" << mean
    << " amp=" << amp
    << " period=" << period
    << " lag=" << lag
    << std::endl;
}
/*********************************************************************\
 * tBlockwaveTimeSeriesImp                                                    *
\*********************************************************************/

double tBlockwaveTimeSeriesImp::calc(double time) const
{
  return sin((time-lag)*TWOPI/period) >= 0.0 ? mean+amp : mean-amp;
}

void tBlockwaveTimeSeriesImp::TellAll() const
{
  std::cout
    << "TS @blockwave."
    << " amp=" << amp
    << " mean=" << mean
    << " period=" << period
    << " lag=" << lag
    << std::endl;
}
/*********************************************************************\
 * tDataTimeSeriesImp                                                    *
\*********************************************************************/


tDataTimeSeriesImp::tDataTimeSeriesImp() :
  xdata(0), ydata(0),
  n(0),
  interpolate(false)
{
}

tDataTimeSeriesImp::~tDataTimeSeriesImp()
{
  if (n) {
    // beware for 1-based arrays
    delete[] ++xdata;
    delete[] ++ydata;
  }
  n = 0;
}

double tDataTimeSeriesImp::calc(double time) const
{
  double x,y,dx,dy;
  int index = -1;

  assert(xdata);
  assert(ydata);

  // hunt trough timeseries data
  nr_hunt(xdata, n, time, index);

  //printf("INDEX: %d %.3f\n", index, xdata[index]);

  if (index==0) {         // before timeseries range
    y = ydata[1];
  } else if (index==n) {  // after timeseries range
    y = ydata[n];
  } else {                // requested time time is between index and index+1
    if (interpolate) {
      dy = ydata[index+1] - ydata[index];
      dx = xdata[index+1] - xdata[index];
      x = time-xdata[index];
      y = (dy/dx) * x;
      y += ydata[index];
    } else { // forward
      y = ydata[index];
    }
  }
  return y;
}

/*********************************************************************\
 * tFileTimeSeriesImp                                                    *
\*********************************************************************/

void tFileTimeSeriesImp::configure(const char *s)
{
  tDataFile datafile;
  tTokList toklist;
  int xcol, ycol;
  char xcolid[64], ycolid[64];
  int ntok;

  // tokenize config string
  toklist.parse(s);
  ntok = toklist.size();
  assert(ntok>=1);

  // file name is first token
  datafile.load(toklist[0]);

  // now analyze tokens according to the following scheme:
  // 1 token : <filename>
  // 2 tokens: <filename> <mode>
  // 3 tokens: <filename> <xcol> <ycol>
  // 4 tokens: <filename> <xcol> <ycol> <mode>
  // <mode> defaults to "forward", options are "forward", "interpolate"
  // <xcol> defaults to "1"
  // <ycol> defaults to "2"

  assert(ntok>=1 && ntok<=4);

  if (ntok==1) {
    xcol = 1;
    ycol = 2;
  } else if (ntok==3 || ntok==4) { // 3 tokens!
    strcpy(xcolid, toklist[1]);
    strcpy(ycolid, toklist[2]);

    if (isdigit(xcolid[0])) {
      sscanf(xcolid,"%d", &xcol);
    } else {
      assert(datafile.hasColumn(xcolid));
      xcol = datafile.getColIndex(xcolid);
    }

    if (isdigit(ycolid[0])) {
      sscanf(ycolid,"%d", &ycol);
    } else {
      assert(datafile.hasColumn(ycolid));
      ycol = datafile.getColIndex(ycolid);
    }
  }

  if (ntok==2 || ntok==4) { // mode specified explicitely
    if (strcmp(toklist[ntok-1], "forward")==0) {
      interpolate = false;
    } else if (strcmp(toklist[ntok-1], "interpolate")==0) {
      interpolate = true;
    } else {
      std::cerr << "tFileTimeSeriesImp.Configure(): invalid mode specifier: `"
		<< toklist[ntok-1] << "'." << std::endl;
      ReportFatalError("Invalid mode.");
    }
  } else {
    interpolate = false;
  }

  // extract column data, convert 1-based to 0-based
  xdata = datafile.getColumn(xcol-1);
  ydata = datafile.getColumn(ycol-1);
  n     = datafile.rows();

  // make 1-based
  xdata--;
  ydata--;

  return;
}

void tFileTimeSeriesImp::TellAll() const
{
  std::cout
    << "TS @file."
    << " add details later."
    << std::endl;
}
/*********************************************************************\
 * tInlineTimeSeriesImp                                                   *
\*********************************************************************/


void tInlineTimeSeriesImp::configure(const char *s)
{
  tTokList toklist;
  int ntok, i;

  // tokenize config string, using custom separators
  toklist.parse(s," \t:");
  ntok = toklist.size();
  assert(ntok>=1);

  // analyze tokens according to the following scheme:
  // even number of tokens: <x1><y1><x2><y2>..<xn><yn>
  // odd  number of tokens: <x1><y1><x2><y2>..<xn><yn><mode>
  // <mode> defaults to "forward", options are "forward", "interpolate"

  // extract table size and optional mode specifier first
  if (ntok%2==0) { // ntok is even
    n = ntok/2;
    interpolate = false;
  } else {         // ntok is odd
    n = (ntok-1)/2;
    if (strcmp(toklist[ntok-1], "forward")==0) {
      interpolate = false;
    } else if (strcmp(toklist[ntok-1], "interpolate")==0) {
      interpolate = true;
    } else {
      std::cerr << "tInlineTimeSeriesImp.Configure(): invalid mode specifier: `"
		<< toklist[ntok-1] << "'." << std::endl;
      ReportFatalError("Invalid mode.");
    }
  }

  // allocate memory
  xdata = new double[n];
  ydata = new double[n];

  // parse table data from tokens
  for (i=0; i<n; i++)
    {
      sscanf(toklist[2*i  ], "%lf", &xdata[i]);
      sscanf(toklist[2*i+1], "%lf", &ydata[i]);
    }

  // make vectors 1-based
  xdata--;
  ydata--;

  return;
}

void tInlineTimeSeriesImp::TellAll() const
{
  std::cout
    << "TS @file."
    << " add details later."
    << std::endl;
}

/*********************************************************************\
 * tTimeSeries: main user interface
\*********************************************************************/

tTimeSeries::tTimeSeries() :
  ts(0)
{}

// copy constructor (SL, 10/10)
tTimeSeries::tTimeSeries(tTimeSeries const &orig) 
  : tagImp(orig.tagImp)
{
  if( tagImp == 0 )
    ts = new tConstantTimeSeriesImp();
  else if( tagImp == 1 )
    ts = new tSinwaveTimeSeriesImp();
  else if( tagImp == 2 )
    ts = new tBlockwaveTimeSeriesImp();
  else if( tagImp == 3 )
    ts = new tLinearGrowthTimeSeriesImp();
  else if( tagImp == 4 )
    ts = new tFileTimeSeriesImp();
  else if( tagImp == 5 )
    ts = new tInlineTimeSeriesImp();
  else
    assert(0);
  ts->Initialize_Copy( orig.ts );
  assert( ts );
}

tTimeSeries::~tTimeSeries()
{
  delete ts;
}

void tTimeSeries::configure(const char *s)
{
  assert(ts == NULL);
  const char *p = 0;

  if (strstr(s,"@sinwave")==s) {
    ts = new tSinwaveTimeSeriesImp();
    p = s+8;
    tagImp = 1;
  }
  else if (strstr(s,"@blockwave")==s) {
    ts = new tBlockwaveTimeSeriesImp();
    p = s+9;
    tagImp = 2;
  }
  else if (strstr(s,"@linear")==s) {
    ts = new tLinearGrowthTimeSeriesImp();
    p = s+7;
    tagImp = 3;
  }
  else if (strstr(s,"@file")==s) {
    ts = new tFileTimeSeriesImp();
    p = s+5;
    tagImp = 4;
  }
  else if (strstr(s,"@inline")==s) {
    ts = new tInlineTimeSeriesImp();
    p = s+7;
    tagImp = 5;
  }
  else if (strstr(s,"@constant")==s) {
    ts = new tConstantTimeSeriesImp();
    p = s+9;
    tagImp = 0;
  }
  else {
    // also parse "normal" lines, containing a float
    double f;
    if (sscanf(s,"%lf",&f)==1) {
      ts = new tConstantTimeSeriesImp();
      p = s;
      tagImp = 0;
    }
  }
  assert(ts!=0);
  ts->configure(p);
}

void tTimeSeries::reconfigure(const char * s)
{
  delete ts;
  ts = NULL;
  configure(s);
}

double tTimeSeries::calc(double time) const
{
  assert(ts != NULL);
  return ts->calc(time);
}

void tTimeSeries::TellAll() const
{
  assert(ts != NULL);
  ts->TellAll();
}

// eof
