//-*-c++-*- 

/*********************************************************************/
/**
**  @file mathutil.h
**  @brief Header file for special math utilities not in math.h.
**         All or most routines from Numerical Recipes in C by
**         Press et al.
**
**  $Id: mathutil.h,v 1.8 2003-08-01 17:14:54 childcvs Exp $
*/
/*********************************************************************/

#ifndef MATHUTIL_H
#define MATHUTIL_H

#define PI 3.14159265358979323846

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

// forward declaration
class tInputFile;
#if !defined(HAVE_NO_NAMESPACE)
# include <fstream>
using namespace std;
#else
# include <fstream.h>
#endif

/** @class tRand
**
**  A simple class that generates a random sequence
*/
class tRand
{
  tRand(tRand const &);
  tRand& operator=(tRand const &);
public:
  tRand(long);
  tRand(tInputFile const &);
  void init(long);
  double ran3();
  void dumpToFile( ofstream&  );
  void readFromFile( ifstream& );
  int numberRecords() const;
private:
  void initFromFile(tInputFile const &);
  // state of ran3()
  long ma[56];
  int inext, inextp;
};


#endif
