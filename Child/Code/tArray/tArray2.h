//-*-c++-*-

/***************************************************************************/
/**
 **  @file tArray2.h
 **  @brief Header file for tArray2 objects.
 **
 **  A tArray2 is an object implementation of a 2D array.
 **  A template is used so that the array may be of any specified type.
 **  It is much faster than tArray as no memory allocation is performed.
 **
 **  Arnaud Desitter - April 2004
 **
 **  $Id: tArray2.h,v 1.1 2004-04-14 12:57:37 childcvs Exp $
 */
/***************************************************************************/

#ifndef TARRAY2_H
#define TARRAY2_H

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include "../errors/errors.h"

/***************************************************************************/
/**
 **  @class tArray2
 **
 **  The tArray2 template class implements 1D arrays (ie, vectors) of any
 **  data type.
 **
 */
/***************************************************************************/
template< class T >
class tArray2
{
  void fatalReport( int ) const ATTRIBUTE_NORETURN; // bail out
public:
  inline tArray2();                      // default constructor
  inline tArray2( const tArray2< T > & ); // copy constructor

  inline tArray2( const T&, const T& );  // Constructor with 2 elements

  inline tArray2< T > &operator=( const tArray2< T > & ); // memberwise assignmt
  bool operator==( const tArray2< T > & ) const;    // memberwise comparison
  inline bool operator!=( const tArray2< T > & ) const;    // memberwise comparison
private:
  // Do not use them (too slow) !
  inline T &operator[]( unsigned int );   // overloaded array index operator
  inline const T &operator[]( unsigned int ) const;
public:
  inline T & at( unsigned int );
  inline const T & at( unsigned int ) const;
  inline T *getArrayPtr();   // returns the actual array; needed for passing
  // to fortran.
  inline const T *getArrayPtr() const; // returns the actual array
private:
  T avalue[2]; // the array itself
};

template< class T >
ostream &operator<<( ostream &output, const tArray2< T > &a );


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
inline tArray2< T >::
tArray2()
{
  avalue[0] = 0;
  avalue[1] = 0;
}

//specialized constructors
template< class T >
inline tArray2< T >::
tArray2( const T& e1, const T& e2 )
{
  avalue[0] = e1;
  avalue[1] = e2;
}

//copy constructor
template< class T >
inline tArray2< T >::tArray2( const tArray2< T > &original )
{
  avalue[0] = original.avalue[0];
  avalue[1] = original.avalue[1];
}

//overloaded assignment operator:
template< class T >
inline tArray2< T > &tArray2< T >::operator=( const tArray2< T > &right )
{
  if( &right != this ){
    avalue[0] = right.avalue[0];
    avalue[1] = right.avalue[1];
  }
  return *this;
}

/**************************************************************************\
 **  Overloaded operators:
 **
\**************************************************************************/
//overloaded equality operator:
template< class T >
bool tArray2< T >::operator==( const tArray2< T > &right ) const
{
  if ( &right == this ) return true;
  if( avalue[0] != right.avalue[0] ) return false;
  if( avalue[1] != right.avalue[1] ) return false;
  return true;
}

//overloaded inequality operator:
template< class T >
inline bool tArray2< T >::operator!=( const tArray2< T > &right ) const
{
  return !operator==(right);
}

//overloaded subscript operator:
template< class T >
inline T &tArray2< T >::operator[]( unsigned int subscript )
{
  if ( unlikely(subscript >= 2) )
    fatalReport( subscript );
  return avalue[subscript];
}

template< class T >
inline const T &tArray2< T >::operator[]( unsigned int subscript ) const
{
  if ( unlikely(subscript >= 2) )
    fatalReport( subscript );
  return avalue[subscript];
}

// likewise with no check
template< class T >
inline T &tArray2< T >::at( unsigned int subscript )
{
  return avalue[subscript];
}

template< class T >
inline const T &tArray2< T >::at( unsigned int subscript ) const
{
  return avalue[subscript];
}

template< class T >
void tArray2< T >::fatalReport( int subscript ) const
{
  cout << "subscript is " << subscript << " > 2" <<endl;
  ReportFatalError( "Subscript out of range." );
}


//overloaded left shift operator
template< class T >
ostream &operator<<( ostream &output, const tArray2< T > &a )
{
  output << a.avalue[0] << ' ' << a.avalue[1] << ' ';
  return output;
}

/**************************************************************************\
 **  getArrayPtr: returns a pointer to the head of the array (needed for
 **               passing arrays to Fortran routines)
\**************************************************************************/
template< class T >
inline T *tArray2< T >::
getArrayPtr() {return avalue;}

template< class T >
inline const T *tArray2< T >::
getArrayPtr() const {return avalue;}

#endif
