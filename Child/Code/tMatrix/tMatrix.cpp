/*************************************************************************/
/**
**  @file tMatrix.cpp
**  @brief Definition of class tMatrix.
**
**  Class tMatrix implements matrices (2D arrays) as 1D arrays of tArray
**  objects. The () operator is overloaded so that you can access array
**  elements in the "old fashioned" way, like "mymatrix( 42, 24 );".
**  One of the constructors takes two integer arguments representing the
**  size of the matrix.
**
**  $Id: tMatrix.cpp,v 1.5 2003-06-06 12:49:32 childcvs Exp $
*/
/*************************************************************************/

#include "tMatrix.h"

// Constructor: sets the size of the matrix to nr by nc and sets all values
// to zero.
template < class T >
tMatrix<T>::tMatrix( int nr, int nc ) :
  ptr(0),
  nrows(nr),
  ncols(nc)
{
   ptr = new tArray<T> [nr];
   for( int i=0; i<nr; i++ )
       ptr[i].setSize(nc);
}

// Destructor: deletes the array of tArray objects (whose destructors take
// care of the rest)
template < class T >
tMatrix<T>::~tMatrix()
{
   delete [] ptr;
}


// Overloaded () operator for referencing individual matrix entries
template < class T >
T &tMatrix<T>::operator()( int row, int col )
{
   return (ptr[row])[col];
}

template < class T >
const T &tMatrix<T>::operator()( int row, int col ) const
{
   return (ptr[row])[col];
}

