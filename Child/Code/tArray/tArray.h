//-*-c++-*- 

/***************************************************************************/
/**
**  @file tArray.h
**  @brief Header file for tArray objects.
**
**  A tArray is an object implementation of a one-dimensional array.
**  A template is used so that the array may be of any specified type.
**  Unlike regular arrays, tArray objects provide bounds checking,
**  memberwise equality/inequality comparison, and memberwise copy
**  operations. The size of the array is determined either by an 
**  argument passed to the constructor or by assignment of one array
**  to another. 
**
**  $Id: tArray.h,v 1.21 2003-08-06 16:13:51 childcvs Exp $
*/
/***************************************************************************/

#ifndef TARRAY_H
#define TARRAY_H

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include "../errors/errors.h"

/***************************************************************************/
/**
**  @class tArray
**
**  The tArray template class implements 1D arrays (ie, vectors) of any
**  data type.
**
*/
/***************************************************************************/
template< class T >
class tArray
{
    //friend ostream &operator<<( ostream &, const tArray< T > & );
    //friend istream &operator>>( istream &, tArray< T > & );
    //friend ofstream &operator<<( ofstream &, const tArray< T > & );
    //friend ifstream &operator>>( ifstream &, tArray< T > & );*/
    void fatalReport( int ) const ATTRIBUTE_NORETURN; // bail out
public:
    tArray();                      // default constructor
    tArray( int );                 // constructor that initializes array size
    tArray( const tArray< T > & ); // copy constructor
    inline ~tArray();              // destructor
    const tArray< T > &operator=( const tArray< T > & ); // memberwise assignmt
    int operator==( const tArray< T > & ) const;    // memberwise comparison
    int operator!=( const tArray< T > & ) const;    // memberwise comparison
    inline T &operator[]( int );   // overloaded array index operator
    inline const T &operator[]( int ) const;
    int getSize() const {       // returns the number of elements in the array
      return npts;
    }
    void setSize( int );       // reinitializes and sets array size
    inline T *getArrayPtr();   // returns the actual array; needed for passing
                               // to fortran.
    inline const T *getArrayPtr() const; // returns the actual array
private:
    int npts;   // size of array
    T * avalue; // the array itself
};

template< class T >
ostream &operator<<( ostream &output, const tArray< T > &a );

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tArray.cpp"
#endif

#endif
