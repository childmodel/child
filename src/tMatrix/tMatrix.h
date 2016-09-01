//-*-c++-*-

/*************************************************************************/
/**
**  @file tMatrix.h
**  @brief Header file for class tMatrix.
**
**  Class tMatrix implements matrices (2D arrays) as 1D arrays of tArray
**  objects. The () operator is overloaded so that you can access array
**  elements in the "old fashioned" way, like "mymatrix( 42, 24 );".
**  One of the constructors takes two integer arguments representing the
**  size of the matrix.
**
**  $Id: tMatrix.h,v 1.14 2004-04-16 18:34:06 childcvs Exp $
*/
/*************************************************************************/
#ifndef TMATRIX_H
#define TMATRIX_H

#include "../tArray/tArray.h"


template < class T >
class tMatrix
{
  tMatrix(const tMatrix&);
  tMatrix& operator=(const tMatrix&);
  tMatrix();
public:
  inline tMatrix( int nr, int nc );
  inline T &operator()( int row, int col );
  inline const T &operator()( int row, int col ) const;
  int getNumRows() const {return nrows;}
  int getNumCols() const {return ncols;}

private:
  tArray<T> ptr;
  const int nrows;
  const int ncols;
};

// template < class T >
// inline tMatrix<T>::tMatrix(const tMatrix& orig)
//   : ptr( orig.ptr ),
//     nrows(orig.nrows),
//     ncols(orig.ncols)
// {}

// template < class T >
// inline tMatrix& tMatrix<T>::operator=(const tMatrix& right)

// Constructor: sets the size of the matrix to nr by nc and sets all values
// to zero.
template < class T >
inline tMatrix<T>::tMatrix( int nr, int nc ) :
  ptr(nr*nc),
  nrows(nr),
  ncols(nc)
{}

// Overloaded () operator for referencing individual matrix entries
template < class T >
inline T &tMatrix<T>::operator()( int row, int col )
{
  return ptr[col+ncols*row];
}

template < class T >
inline const T &tMatrix<T>::operator()( int row, int col ) const
{
   return ptr[col+ncols*row];
}

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tMatrix.cpp"
#endif

#endif
