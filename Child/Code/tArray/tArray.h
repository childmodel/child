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
**  $Id: tArray.h,v 1.1 1998-01-14 20:12:51 gtucker Exp $
\***************************************************************************/

#ifndef TARRAY_H
#define TARRAY_H

/** class tArray  **********************************************************/
template< class T >
class tArray
{
     /*friend ostream &operator<<( ostream &, const tArray< T > & );
   friend istream &operator>>( istream &, tArray< T > & );
   friend ofstream &operator<<( ofstream &, const tArray< T > & );*/
     /*friend ifstream &operator>>( ifstream &, tArray< T > & );*/
  public:
   tArray();
   tArray( int );
   tArray( const tArray< T > & );
   ~tArray();
   const tArray< T > &operator=( const tArray< T > & ); // memberwise assignmt
   int operator==( const tArray< T > & ) const;    // memberwise comparison
   int operator!=( const tArray< T > & ) const;    // memberwise comparison
   T &operator[]( int );               
   int getSize();             // returns the number of elements in the array
   
  private:
   int npts;
   T * avalue;
};

#endif
