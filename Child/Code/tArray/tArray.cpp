/**************************************************************************\
**
**  tArray.cpp
**
**  Functions for class tArray< T >
**
**  $Id: tArray.cpp,v 1.3 1998-02-12 23:24:08 stlancas Exp $
\**************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "tArray.h"

/**************************************************************************\
**
**  Constructors
**
**  (1) default constructor sets size to zero and pointer to null.
**  (2) creates an array of specified size and initializes values to zero.
**  (3) copy constructor
**
\**************************************************************************/
template< class T >                                               //tArray
tArray< T >::
tArray()
{
   npts = 1;
   avalue = new T [1];
   assert( avalue != 0 );
   avalue[0] = 0;
     //cout<<"tArray(): one member array"<<npts<<endl<<flush;
}

template< class T >                                               //tArray
tArray< T >::
tArray( int number )
{
   assert( number > 0 );
   int i;
   npts = number;
   avalue = new T [npts];
   assert( avalue != 0 );
   for( i=0; i<npts; i++ )
       avalue[i] = 0;
     //cout<<"tArray(npts): no. in array "<<npts<<endl<<flush;
}

//copy constructor
template< class T >                                               //tArray
tArray< T >::
tArray( const tArray< T > &original )
{
   int i;
   if( &original != 0 )
   {
      cout << flush;
      /*for( i = 0; i < original.npts; i++ )
      {
         cout << original.avalue[i] << " ";
      }
      cout << endl << flush;*/
      
      npts = original.npts;
      assert( npts > 0 );
      avalue = new T[npts];
      assert( avalue != 0 );
      for( int i = 0; i < npts; i++ )
          avalue[i] = original.avalue[i];
   }
     //cout<<"tArray(original): no. in array "<<npts<<endl<<flush;
}

/**************************************************************************\
**  Destructor
\**************************************************************************/
template< class T >                                               //tArray
tArray< T >::
~tArray()
{
     //cout<<"~tArray: no. in array "<<npts<<endl<<flush;
   delete [] avalue;
}

/**************************************************************************\
**  Overloaded operators:
**    assignment, equality, inequality: memberwise operation
**    index: uses an assertion to check array bounds (assumed to be within
**           bounds at runtime) 
\**************************************************************************/
//overloaded assignment operator:
template< class T >                                               //tArray
const tArray< T > &tArray< T >::operator=( const tArray< T > &right )
{
   assert( &right != 0 );
   if( &right != this )
   {
      delete [] avalue;
      npts = right.npts;
      avalue = new T [npts];
      assert( avalue != 0 && right.avalue != 0 && npts != 0 );
        //cout << "tArray op=: npts " << npts << "; ";
      for( int i = 0; i < npts; i++ )
      {
           //cout << right.avalue[i] << " ";
         avalue[i] = right.avalue[i];
      }
        //cout << endl;
   }
   return *this;
}

//overloaded equality operator:
template< class T >                                               //tArray
int tArray< T >::operator==( const tArray< T > &right ) const
{
   if( npts != right.npts ) return 0;
   for( int i = 0; i < npts; i++ )
       if( avalue[i] != right.avalue[i] )
           return 0;
   return 1;
}

//overloaded inequality operator:
template< class T >                                               //tArray
int tArray< T >::operator!=( const tArray< T > &right ) const
{
   if( npts != right.npts ) return 0;
   for( int i = 0; i < npts; i++ )
       if( avalue[i] != right.avalue[i] )
           return 1;
   return 0;
}

//overloaded subscript operator:
template< class T >                                               //tArray
T &tArray< T >::operator[]( int subscript )
{
   assert( 0 <= subscript && subscript < npts );
   return avalue[subscript];
}

//overloaded input operator for class tArray< T >:
/*template< class T >                                               //tArray
istream &operator>>( istream &input, tArray< T > &a )
{operator<<( ostream &output, const tArray< T > &a )
{
   for( int i = 0; i < a.npts; i++ )
   {
      output << a.avalue[i] << " ";
      if( (i + 1) %10 == 0 ) output<<endl;
   }
   if( i % 10 != 0 ) output<<endl;
   return output;
}*/
//overloaded input file operator for class tArray< T >:
/*template< class T >                                               //tArray
ifstream &operator>>( ifstream &input, tArray< T > &a )
{
   for( int i = 0; i < a.npts; i++ )
       input >> a.avalue[i];
   return input;
}*/
//overloaded output file operator for class tArray< T >:
/*template< class T >                                               //tArray
ofstream &operator<<( ofstream &output, const tArray< T > &a )
{
   for( int i = 0; i < a.npts; i++ )
   {
      output << a.avalue[i] << " ";
      if( (i + 1) %10 == 0 ) output<<endl;
   }
   if( i % 10 != 0 ) output<<endl;
   return output;
}*/

/**************************************************************************\
**  getSize: returns the number of elements in the array 
\**************************************************************************/
template< class T >                                               //tArray
int tArray< T >::
getSize() {return npts;}

template< class T >
T *tArray< T >::
getArrayPtr() {return avalue;}
