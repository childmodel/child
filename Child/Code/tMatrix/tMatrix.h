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
**  $Id: tMatrix.h,v 1.10 2003-01-17 17:30:34 childcvs Exp $
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
public:
  tMatrix();
  tMatrix( int nr, int nc );
  ~tMatrix();
  T &operator()( int row, int col );
  int getNumRows() const {return nrows;}
  int getNumCols() const {return ncols;}
    
private:
  tArray<T> *ptr;
  int nrows;
  int ncols;  
};

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tMatrix.cpp"
#endif

#endif
