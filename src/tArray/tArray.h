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
 **  $Id: tArray.h,v 1.29 2004-06-16 13:37:30 childcvs Exp $
 */
/***************************************************************************/

#ifndef TARRAY_H
#define TARRAY_H

#include <iosfwd>
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
  void fatalReport( size_t ) const ATTRIBUTE_NORETURN; // bail out
public:
  inline tArray();               // default constructor
  tArray( size_t );              // constructor that initializes array size
  tArray( size_t, const T& );
  tArray( const tArray< T > & ); // copy constructor

  inline tArray( const T&, const T& );  // Array of size 2 with 2 elements
  inline tArray( const T&, const T&, const T& );  // Array of size 3 with 3 elements

  inline ~tArray();              // destructor
  const tArray< T > &operator=( const tArray< T > & ); // memberwise assignmt
  bool operator==( const tArray< T > & ) const;    // memberwise comparison
  bool operator!=( const tArray< T > & ) const;    // memberwise comparison
  inline T &operator[]( size_t );   // overloaded array index operator
  inline const T &operator[]( size_t ) const;
  inline T & at( size_t );
  inline const T & at( size_t ) const;
  size_t getSize() const {       // returns the number of elements in the array
    return npts;
  }
  void setSize( size_t );       // reinitializes and sets array size
  inline T *getArrayPtr();   // returns the actual array; needed for passing
  // to fortran.
  inline const T *getArrayPtr() const; // returns the actual array
private:
  T * avalue; // the array itself
  size_t npts;   // size of array
};

template< class T >
std::ostream &operator<<( std::ostream &output, const tArray< T > &a );


/**************************************************************************\
 **
 **  Constructors & destructors:
 **
 **  (1) default constructor - sets size to zero and pointer to null
 **  (2) creates an array of specified size and initializes values to zero
 **  (3) copy constructor - makes copy; assumes original array not empty
 **
 **  Destructor deletes the array elements.
 **
\**************************************************************************/
//default constructor
template< class T >
inline tArray< T >::
tArray() :
  avalue(0), npts(1)
{
  avalue = new T [1];
  avalue[0] = 0;
}

//specialized constructors
template< class T >
inline tArray< T >::
tArray( const T& e1, const T& e2 ) :
  avalue(0), npts(2)
{
  avalue = new T [2];
  avalue[0] = e1;
  avalue[1] = e2;
}

template< class T >
inline tArray< T >::
tArray( const T& e1, const T& e2, const T& e3 ) :
  avalue(0), npts(3)
{
  avalue = new T [3];
  avalue[0] = e1;
  avalue[1] = e2;
  avalue[2] = e3;
}

template< class T >
inline tArray< T >::
~tArray()
{
  delete [] avalue;
}

/**************************************************************************\
 **  Overloaded operators:
 **
 **    index - uses an assertion to check array bounds (assumed to be within
 **           bounds at runtime) and returns value
 **
\**************************************************************************/
//overloaded subscript operator:
template< class T >
inline T &tArray< T >::operator[]( size_t subscript )
{
  if ( unlikely(subscript >= npts) )
    fatalReport( subscript );
  return avalue[subscript];
}

template< class T >
inline const T &tArray< T >::operator[]( size_t subscript ) const
{
  if ( unlikely(subscript >= npts) )
    fatalReport( subscript );
  return avalue[subscript];
}

// likewise with no check
template< class T >
inline T &tArray< T >::at( size_t subscript )
{
  return avalue[subscript];
}

template< class T >
inline const T &tArray< T >::at( size_t subscript ) const
{
  return avalue[subscript];
}

/**************************************************************************\
 **  getArrayPtr: returns a pointer to the head of the array (needed for
 **               passing arrays to Fortran routines)
\**************************************************************************/
template< class T >
inline T *tArray< T >::
getArrayPtr() {return avalue;}

template< class T >
inline const T *tArray< T >::
getArrayPtr() const {return avalue;}

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tArray.cpp"
#endif

#endif
