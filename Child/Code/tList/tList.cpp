/**************************************************************************\
**
**  Functions for class tList< NodeType >
**
**  $Id: tList.cpp,v 1.1 1998-01-14 20:34:17 gtucker Exp $
\**************************************************************************/
#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "../tListNode/tListNode.h"
#include "tList.h"

//default constructor
template< class NodeType >                         //tList
tList< NodeType >::tList()
{
   first = last = 0;
   nNodes = 0;
     //cout << "list instantiated" << first << endl;
}
template< class NodeType >                         //tList
tList< NodeType >::
tList( const tList< NodeType > *original )
{
   assert( original != 0 );
     //nNodes = original->nNodes;
   tListNode< NodeType > * current = original->first;
   for( int i=0; i<nNodes; i++ )
   {
      insertAtBack( current->getDataRef() );
      current = current->next;
   }
   assert( nNodes == original->nNodes );
     //cout << "list copy instantiated" << first << endl;
}


template< class NodeType >                         //tList
tList< NodeType >::
~tList()
{
   if( !isEmpty() )
   {
        //cout<<"Destroying nodes ... "<<endl;
      tListNode<NodeType > * current = first, * temp;
      while( current != 0 )
      {
         temp = current;
         //cout<<temp->data<<endl;
         current = current->next;
         delete temp;
      }
   }
     //cout<<"All nodes destroyed"<<endl<<endl;
}

//overloaded assignment operator
template< class NodeType >                         //tList
const tList< NodeType > &tList< NodeType >::
operator=( const tList< NodeType > &right )
{
   if( this != &right )
   {
        //create an equivalent object rather than pointers to the original
      Flush();
      tListNode< NodeType > *cn = right.first;
      if( cn != 0 )
      {
         insertAtBack( cn->data );
         for( cn = cn->next; cn != last->next; cn = cn->next )
             insertAtBack( cn->data );
      }
        /*first = right.first;
      last = right.last;
      nNodes = right.nNodes;*/
   }
     //cout << "list assigned" << first << endl;
   return *this;
}

//overloaded equality operator:
template< class NodeType >                         //tList
int tList< NodeType >::
operator==( const tList< NodeType > &right ) const
{
   if( nNodes != right.nNodes ) return 0;
   if( first != right.first ) return 0;
   if( last != right.last ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                         //tList
int tList< NodeType >::
operator!=( const tList< NodeType > &right ) const
{
   if( nNodes != right.nNodes ) return 1;
   if( first != right.first ) return 1;
   if( last != right.last ) return 1;
   return 0;
}

//return pointer to newly allocated node
template< class NodeType >                         //tList
tListNode< NodeType > * tList< NodeType >::
getNewNode( const NodeType &value )
{
   tListNode< NodeType > * ptr =
       new tListNode< NodeType >( value );
   assert( ptr != 0 );
     //cout << "new node created" << endl;
   nNodes++;
   return ptr;
}

//insert at front
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtFront( const NodeType &value )
{
   tListNode< NodeType > *newPtr = getNewNode( value );
   if( isEmpty() ) first = last = newPtr;
   else
   {
      newPtr->next = first;
      if( last->next == first ) last->next = newPtr;
      first = newPtr;
   }
     //nNodes++;
}

//insert at back
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtBack( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node to list" << endl;
   assert( this != 0 );
   if( isEmpty() )
   {
      first = last = newPtr;
   }
   else
   {
      newPtr->next = last->next;
      last->next = newPtr;
      last = newPtr;
   }
     //nNodes++;
     //cout << "added node to back of list" << endl;
}

//insert at next spot in list
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtNext( const NodeType &value, tListNode< NodeType > * prev )
{
   if( prev != 0 )
   {
      if( prev == last )
      {
         insertAtBack( value );
         return;
      }
      tListNode< NodeType > * newPtr = getNewNode( value );
      newPtr->next = prev->next;
      prev->next = newPtr;
        //nNodes++;
   }
}

//insert at previous spot in list
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtPrev( const NodeType &value, tListNode< NodeType > * node )
{
   tListNode< NodeType > * prev;
   if( node != 0 )
   {
      if( node == first )
      {
         insertAtFront( value );
         return;
      }
      tListNode< NodeType > * newPtr = getNewNode( value );
      for( prev = first; prev->next != node; prev = prev->next ); 
      newPtr->next = prev->next;
      prev->next = newPtr;
        //nNodes++;
   }
}

//delete from front
template< class NodeType >                         //tList
int tList< NodeType >::
removeFromFront( NodeType &value )
{
   if( isEmpty() ) return 0;
   else
   {
      tListNode< NodeType > * temp = first;
      if( first == last ) first = last = 0;
      else
      {
         if( last->next == first ) last->next = first->next;
         first = first->next;
      }
      value = temp->data;
      delete temp;
      nNodes--;
      return 1;
   }
}

//delete from back
template< class NodeType >                         //tList
int tList< NodeType >::
removeFromBack( NodeType &value )
{
   if( isEmpty() ) return 0;
   else
   {
      tListNode< NodeType > * temp = last;
      if( first == last ) first = last = 0;
      else
      {
         tListNode< NodeType > * current = first;
         while( current->next != last ) current = current->next;
         current->next = last->next;
         last = current;
      }
      value = temp->data;
      delete temp;
      nNodes--;
      return 1;
   }
}

//delete next node
template< class NodeType >                         //tList
int tList< NodeType >::
removeNext( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr->next == 0 ) return 0;
   if( ptr == 0 ) return 0;
   if( ptr->next == last ) return removeFromBack( value );
   else if( ptr->next == first ) return removeFromFront( value );
   //if( ptr == last ) return 0;
   tListNode< NodeType > * temp = ptr->next;
   ptr->next = ptr->next->next;
   value = temp->data;
   delete temp;
   nNodes--;
   return 1;
}

//delete previous node
template< class NodeType >                         //tList
int tList< NodeType >::
removePrev( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr == first && last->next == 0 ) return 0;
   if( ptr == first ) return removeFromBack( value );
   tListNode< NodeType > * temp, *prev;
   for( temp = first; temp->next->next != ptr; temp = temp->next );
   prev = temp;
   temp = temp->next;
   if( temp == first ) return removeFromFront( value );
   prev->next = prev->next->next;
   value = temp->data;
   delete temp;
   nNodes--;
   return 1;
}

//remove all nodes
template< class NodeType >                         //tList
void tList< NodeType >::
Flush()
{
   if( !isEmpty() )
   {
        //cout<<"Destroying nodes ... "<<endl;
      tListNode<NodeType > * current = first, * temp;
      while( current != 0 )
      {
         temp = current;
           //cout<<temp->data<<endl;
         current = current->next;
         delete temp;
      }
   }
   assert( isEmpty() );
   nNodes = 0;
     //cout<<"All nodes destroyed"<<endl<<endl;
}

//empty?
template< class NodeType >                         //tList
int tList< NodeType >::
isEmpty() const
{
     //cout << "checking if list empty" << endl;
   if( first == 0 )
   {
        //cout << "list is empty" << endl;
      return 1;
   }
   else
   {
        //cout << "list is not empty" << endl;
      return 0;
   }
}

//display list contents
template< class NodeType >                         //tList
void tList< NodeType >::
print() const
{
   if( isEmpty() )
   {
      cout<<"The list is empty"<<endl<<endl;
      return;
   }
   tListNode< NodeType > * current = first;
   cout<<"The list is: ";
   while( current != 0 )
   {
      //cout<<current->data<<' ';
      current = current->next;
   }
   cout<<endl<<endl;
}

//input list contents
/*template< class NodeType >                         //tList
void tList< NodeType >::
input( int size )
{
   NodeType tempnode;
   for( int i=0; i<size; i++ )
   {
      cout<<"input node data:" << endl;
      cin >> tempnode;
      assert( &tempnode != 0 );
      tempnode.setID( i );
      insertAtBack( tempnode );
   }
   nNodes = size;
}*/

//'get' utilities
//return size
template< class NodeType >                         //tList
int tList< NodeType >::
getSize() const {return nNodes;}

template< class NodeType >                         //tList
tListNode< NodeType > * tList< NodeType >::
getFirst() const {return first;}

template< class NodeType >                         //tList
tListNode< NodeType > * tList< NodeType >::
getLast() const {return last;}

//'move' utilities
template< class NodeType >                         //tList
void tList< NodeType >::
moveToBack( tListNode< NodeType > * mvnode ) 
{
   if( mvnode != last )
   {  
      if( mvnode == first ) first = first->next;
      else
      {
         for( tListNode< NodeType > * prev = first;
              prev->next != mvnode;
              prev = prev->next );
         prev->next = mvnode->next;
      }
      mvnode->next = last->next;
      last->next = mvnode;
      last = mvnode;
      if( last->next != 0 ) last->next = first;
   }
}

template< class NodeType >                         //tList
void tList< NodeType >::
moveToFront( tListNode< NodeType > * mvnode ) 
{
   if( mvnode != first )
   {
      for( tListNode< NodeType > * prev = first;
           prev->next != mvnode;
           prev = prev->next );
      prev->next = mvnode->next;
      mvnode->next = first;
      first = mvnode;
      if( last == mvnode ) last = prev;
      if( last->next != 0 ) last->next = first;
   }
}

template< class NodeType >                         //tList
void tList< NodeType >::
setNNodes( int val ) {nNodes = ( val >= 0 ) ? val : 0;}

template< class NodeType >                         //tList
void tList< NodeType >::
makeCircular() {last->next = first;}

template< class NodeType >                         //tList
const NodeType tList< NodeType >::
getIthData( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getData();
}

template< class NodeType >                         //tList
const NodeType &tList< NodeType >::
getIthDataRef( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataRef();
}

template< class NodeType >                         //tList
const NodeType *tList< NodeType >::
getIthDataPtr( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0;i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataPtr();
}
//set
template< class NodeType >                         //tList
NodeType tList< NodeType >::
getIthDataNC( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getData();
}

template< class NodeType >                         //tList
NodeType &tList< NodeType >::
getIthDataRefNC( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataRefNC();
}

template< class NodeType >                         //tList
NodeType *tList< NodeType >::
getIthDataPtrNC( int num ) const
 {
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataPtrNC();
}

