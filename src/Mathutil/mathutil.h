//-*-c++-*- 

/*********************************************************************/
/**
**  @file mathutil.h
**  @brief Header file for special math utilities not in math.h.
**         All or most routines from Numerical Recipes in C by
**         Press et al.
**
**  SL, 8/10: Added ExpDev as member function for use by other
**  objects (other than tStorm) that need it (tFire, tForest).
**
**  $Id: mathutil.h,v 1.11 2004-06-16 13:37:27 childcvs Exp $
*/
/*********************************************************************/

#ifndef MATHUTIL_H
#define MATHUTIL_H

#define PI 3.14159265358979323846
#define TWOPI (2*PI)

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
#include <iosfwd>
#include <math.h>

/** @class tRand
**
**  A simple class that generates a random sequence
*/
class tRand
{
public:
  tRand(tRand const &);
  tRand& operator=(tRand const &);
  tRand();
  tRand(long);
  tRand(tInputFile const &);
  void init(long);
  double ran3();
  double ExpDev();
  void dumpToFile( std::ofstream&  );
  void readFromFile( std::ifstream& );
  int numberRecords() const;
private:
  void initFromFile(tInputFile const &);
  // state of ran3()
  long ma[56];
  int inext, inextp;
};


#endif
