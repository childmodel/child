/**************************************************************************/
/**
**  @file mathutil.cpp
**  @brief Special math routines not in math libraries. Most or all
**         from Numerical Recipes in C by Press et al.
**
**  $Id: mathutil.cpp,v 1.5 2003-08-01 17:14:54 childcvs Exp $
*/
/**************************************************************************/

#include "mathutil.h"

#include "mt19937ar-cok.cpp"

#include "../tInputFile/tInputFile.h"


/*********************************************************\
**  ran3
**
**  Random number generator from Numerical Recipes.
**  Returns a uniform random number between 0.0 and 1.0.
**  Set idum to any negative value to initialize or
**  reinitialize the sequence.
**
**  Parameters: idum - random seed
**
\*********************************************************/

tRand::tRand(long seed)
{
  init(seed);
}

tRand::tRand( tInputFile const &infile )
{
  initFromFile( infile );
}

void tRand::initFromFile( tInputFile const &infile )
{
  int seed;
  seed = infile.ReadItem( seed, "SEED" );
  init(seed);
}

void tRand::dumpToFile( ofstream& outFile ){
  for(size_t i=1; i<sizeof(ma)/sizeof(ma[0]); ++i)
    outFile << ma[i] << '\n';
  outFile << inext << '\n' << inextp << '\n';
}

void tRand::readFromFile( ifstream& inFile ){
  for(size_t i=1; i<sizeof(ma)/sizeof(ma[0]); ++i)
    inFile >> ma[i];
  inFile >> inext;
  inFile >> inextp;
}

int tRand::numberRecords() const {
  return sizeof(ma)/sizeof(ma[0])-1+2;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

void tRand::init(long seed)
{
  int i,ii,k;

  long mj=MSEED-(seed < 0 ? -seed : seed);
  mj %= MBIG;
  ma[55]=mj;
  long mk=1;
  for (i=1;i<=54;i++) {
    ii=(21*i) % 55;
    ma[ii]=mk;
    mk=mj-mk;
    if (mk < MZ) mk += MBIG;
    mj=ma[ii];
  }
  for (k=1;k<=4;k++)
    for (i=1;i<=55;i++) {
      ma[i] -= ma[1+(i+30) % 55];
      if (ma[i] < MZ) ma[i] += MBIG;
    }
  inext=0;
  inextp=31;
}


double tRand::ran3()
{
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  long mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
