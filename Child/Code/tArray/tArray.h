/***************************************************************************\
**
**  tArray.h
**
**  Header file for tArray objects.
**
**  A tArray is an object implementation of a one-dimensional array.
**  A template is used so that the array may be of any specified type.
**  Unlike regular arrays, tArray objects provide bounds checking,
**  memberwise equality/inequality comparison, and memberwise copy
**  operations. The size of the array is determined either by an 
**  argument passed to the constructor or by assignment of one array
**  to another. 
**
**  $Id: tArray.h,v 1.8 1999-01-12 21:21:38 gtucker Exp $
\***************************************************************************/

#ifndef TARRAY_H
#define TARRAY_H

#include <iostream.h>


/***************************************************************************\
**
**  class tArray
**
**  The tArray template class implements 1D arrays (ie, vectors) of any
**  data type.
**
\***************************************************************************/
template< class T >
class tArray
{
    friend ostream &operator<<( ostream &, const tArray< T > & );
    //friend istream &operator>>( istream &, tArray< T > & );
    //friend ofstream &operator<<( ofstream &, const tArray< T > & );
    //friend ifstream &operator>>( ifstream &, tArray< T > & );*/
public:
    tArray();                      // default constructor
    tArray( int );                 // constructor that initializes array size
    tArray( const tArray< T > & ); // copy constructor
    ~tArray();                     // destructor
    const tArray< T > &operator=( const tArray< T > & ); // memberwise assignmt
    int operator==( const tArray< T > & ) const;    // memberwise comparison
    int operator!=( const tArray< T > & ) const;    // memberwise comparison
    T &operator[]( int );      // overloaded array index operator
    int getSize() const;       // returns the number of elements in the array
    void setSize( int );       // reinitializes and sets array size
    T *getArrayPtr();          // returns the actual array; needed for passing
                               // to fortran.
private:
    int npts;   // size of array
    T * avalue; // the array itself
};

#endif
