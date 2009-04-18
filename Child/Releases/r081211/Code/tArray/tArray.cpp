/**************************************************************************/
/**
 **  @file tArray.cpp
 **  @brief Functions for template class tArray< T >
 **
 **  $Id: tArray.cpp,v 1.28 2004/06/16 13:37:29 childcvs Exp $
 */
/**************************************************************************/

#include <iostream>
#include <fstream>
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

//constructor that initializes array size
template< class T >
tArray< T >::
tArray( size_t number ) :
  avalue(0), npts(number)
{
  assert( number > 0 );
  avalue = new T [npts];
  for( size_t i=0; i<npts; i++ )
    avalue[i] = 0;
}

template< class T >
tArray< T >::
tArray( size_t number, const T &init) :
  avalue(0), npts(number)
{
  assert( number > 0 );
  avalue = new T [npts];
  for( size_t i=0; i<npts; i++ )
    avalue[i] = init;
}

//copy constructor
template< class T >
tArray< T >::
tArray( const tArray< T > &original ) :
  avalue(0), npts(original.npts)
{
  assert( npts > 0 );

  avalue = new T[npts];
  for( size_t i = 0; i < npts; i++ )
    avalue[i] = original.avalue[i];
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
	  for( size_t i = 0; i < npts; i++ )
	    avalue[i] = right.avalue[i];
	}
    }
  return *this;
}

//overloaded equality operator:
template< class T >
bool tArray< T >::operator==( const tArray< T > &right ) const
{
  if( npts != right.npts ) return false;
  for( size_t i = 0; i < npts; i++ )
    if( avalue[i] != right.avalue[i] )
      return false;
  return true;
}

//overloaded inequality operator:
template< class T >
bool tArray< T >::operator!=( const tArray< T > &right ) const
{
  return !operator==(right);
}

// error handler (do not make it inline)
template< class T >
void tArray< T >::fatalReport( size_t subscript ) const
{
  std::cout<<"subscript is "<<subscript<<" npts is "<<npts<<std::endl;
  ReportFatalError( "Subscript out of range." );
}


//overloaded left shift operator
template< class T >
std::ostream &operator<<( std::ostream &output, const tArray< T > &a )
{
  size_t i;

  for( i = 0; i < a.npts; i++ )
    {
      output << a.avalue[i] << " ";
      if( (i + 1) %10 == 0 ) output<<std::endl;
    }
  if( i % 10 != 0 ) output<<std::endl;
  return output;
}

/**************************************************************************\
 **  setSize: reinitializes and resizes the array
\**************************************************************************/
template< class T >
void tArray<T>::setSize( size_t size )
{
  delete [] avalue;
  npts = size;
  avalue = new T [npts];
  for( size_t i=0; i<npts; i++ ) avalue[i] = 0;
}
