/*************************************************************************\
**
**  tMatrix.h:  Header file for class tMatrix.
**
**  Class tMatrix implements matrices (2D arrays) as 1D arrays of tArray
**  objects. The () operator is overloaded so that you can access array
**  elements in the "old fashioned" way, like "mymatrix( 42, 24 );".
**  One of the constructors takes two integer arguments representing the
**  size of the matrix.
**
**  $Id: tMatrix.h,v 1.3 1999-01-05 22:10:04 stlancas Exp $
\*************************************************************************/
#ifndef TMATRIX_H
#define TMATRIX_H

#include "../tArray/tArray.h"


template < class T >
class tMatrix
{
public:
    tMatrix();
    tMatrix( int nr, int nc );
    ~tMatrix();
    T &operator()( int row, int col );
   int getNumRows() {return nrows;}
   int getNumCols() {return ncols;}
    
private:
    tArray<T> *ptr;
    int nrows;
    int ncols;  
};

#endif
