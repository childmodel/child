/**************************************************************************/
/**
 **  @file tArray.cpp
 **  @brief Functions for template class tArray< T >
 **
 **  $Id: tArray.cpp,v 1.22 2003-08-12 13:30:29 childcvs Exp $
 */
/**************************************************************************/

#if !defined(HAVE_NO_NAMESPACE)
# include <fstream>
using namespace std;
#else
# include <fstream.h>
#endif
#include <assert.h>
#include "tArray.h"
#include "../compiler.h"

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
  npts(1), avalue(0)
{
  avalue = new T [1];
  avalue[0] = 0;
}

//constructor that initializes array size
template< class T >
tArray< T >::
tArray( int number ) :
  npts(number), avalue(0)
{
  assert( number > 0 );
  avalue = new T [npts];
  for( int i=0; i<npts; i++ )
    avalue[i] = 0;
}

//copy constructor
template< class T >
tArray< T >::
tArray( const tArray< T > &original ) :
  npts(original.npts), avalue(0)
{
  assert( npts > 0 );

  avalue = new T[npts];
  for( int i = 0; i < npts; i++ )
    avalue[i] = original.avalue[i];
}

//specialized constructors
template< class T >
inline tArray< T >::
tArray( const T& e1, const T& e2 ) :
  npts(2), avalue(0)
{
  avalue = new T [2];
  avalue[0] = e1;
  avalue[1] = e2;
}

template< class T >
inline tArray< T >::
tArray( const T& e1, const T& e2, const T& e3 ) :
  npts(3), avalue(0)
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
 **    assignment, equality, inequality - memberwise operation
 **    index - uses an assertion to check array bounds (assumed to be within
 **           bounds at runtime) and returns value
 **    left shift - sends array values to output stream, separated by
 **                 spaces and with a carriage return after every 10 items
 **
 **  Modifications:
 **   - Assignment: now allows assignment of empty arrays - GT 7/98
 **   - Do not reallocate if the size is the same - AD 8/2003
 **
\**************************************************************************/

//overloaded assignment operator:
template< class T >
const tArray< T > &tArray< T >::operator=( const tArray< T > &right )
{
  if( &right != this )
    {
      if (npts != right.npts) { // delete and reallocate
	delete [] avalue; avalue = 0;
	npts = right.npts;
	if( npts>0 )
	{
	  assert( right.avalue != 0 );
	  avalue = new T [npts];
	}
      }
      if( npts>0 )
	{
	  for( int i = 0; i < npts; i++ )
	    avalue[i] = right.avalue[i];
	}
    }
  return *this;
}

//overloaded equality operator:
template< class T >
int tArray< T >::operator==( const tArray< T > &right ) const
{
  if( npts != right.npts ) return 0;
  int i;
  for( i = 0; i < npts; i++ )
    if( avalue[i] != right.avalue[i] )
      return 0;
  return 1;
}

//overloaded inequality operator:
template< class T >
int tArray< T >::operator!=( const tArray< T > &right ) const
{
  return !operator==(right);
}

// error handler (do not make it inline)
template< class T >
void tArray< T >::fatalReport( int subscript ) const
{
  cout<<"subscript is "<<subscript<<" npts is "<<npts<<endl;
  ReportFatalError( "Subscript out of range." );
}

//overloaded subscript operator:
template< class T >
inline T &tArray< T >::operator[]( int subscript )
{
  if ( unlikely(0 > subscript || subscript >= npts) )
    fatalReport( subscript );
  return avalue[subscript];
}

template< class T >
inline const T &tArray< T >::operator[]( int subscript ) const
{
  if ( unlikely(0 > subscript || subscript >= npts) )
    fatalReport( subscript );
  return avalue[subscript];
}

//overloaded left shift operator
template< class T >
ostream &operator<<( ostream &output, const tArray< T > &a )
{
  int i;

  for( i = 0; i < a.npts; i++ )
    {
      output << a.avalue[i] << " ";
      if( (i + 1) %10 == 0 ) output<<endl;
    }
  if( i % 10 != 0 ) output<<endl;
  return output;
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

/**************************************************************************\
 **  setSize: reinitializes and resizes the array
\**************************************************************************/
template< class T >
void tArray<T>::setSize( int size )
{
  assert( size>=0 );
  delete [] avalue;
  npts = size;
  avalue = new T [npts];
  for( int i=0; i<npts; i++ ) avalue[i] = 0;
}
