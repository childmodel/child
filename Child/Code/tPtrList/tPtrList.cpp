/**************************************************************************\
**
**  tPtrList.cpp: Functions for classes tPtrList and tPtrListNode.
**
**  $Id: tPtrList.cpp,v 1.1 1998-01-21 00:49:35 gtucker Exp $
\**************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>

#include "tPtrList.h"

/**************************************************************************\
**
**         Utilities for class tPtrListNode
**
\**************************************************************************/
template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::tPtrListNode()
{
   Ptr = 0;
   next = 0;
}

template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::
tPtrListNode( const tPtrListNode< NodeType > & init )
{
   if( &init != 0 )
   {
      Ptr = init.Ptr;
      next = init.next;
   }
}

template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::
tPtrListNode( NodeType * NTPtr )
{
   Ptr = NTPtr;
   next = 0;
}

template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::
~tPtrListNode()
{
   Ptr = 0;  //redirect the pointer away from its target
   next = 0;
}

template< class NodeType >                  //tPtrListNode
const tPtrListNode< NodeType > &tPtrListNode< NodeType >::
operator=( const tPtrListNode< NodeType > &right )
{
   if( &right != this )
   {
      assert( &right != 0 );
      Ptr = right.Ptr;
      next = right.next;
   }
   return *this;
}

//overloaded equality operator:
template< class NodeType >                  //tPtrListNode
int tPtrListNode< NodeType >::
operator==( const tPtrListNode< NodeType > &right ) const
{
   if( next != right.next ) return 0;
   if( Ptr != right.Ptr ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                  //tPtrListNode
int tPtrListNode< NodeType >::
operator!=( const tPtrListNode< NodeType > &right ) const
{
   if( next != right.next ) return 1;
   if( Ptr != right.Ptr ) return 1;
   return 0;
}

template< class NodeType >                  //tPtrListNode
const NodeType * tPtrListNode< NodeType >::
getPtr() const {return Ptr;}

template< class NodeType >                  //tPtrListNode
NodeType * tPtrListNode< NodeType >::
getPtrNC() {return Ptr;}

template< class NodeType >                  //tPtrListNode
const tPtrListNode< NodeType > *
tPtrListNode< NodeType >::
getNext() const {return next;}

template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType > *
tPtrListNode< NodeType >::
getNextNC( ) {return next;}

/*#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "../Definitions.h"
#include "../Classes.h"
#include "../tPtrListNode/tPtrListNode.h"
#include "tPtrList.h"*/

/**************************************************************************\
**
**         Utilities for class tPtrList
**
\**************************************************************************/
template< class NodeType >                      //tPtrList
tPtrList< NodeType >::
tPtrList()
{
   first = last = 0;
   nNodes = 0;
     //cout << "tPtrList() instantiated" << endl;
}

template< class NodeType >                      //tPtrList
tPtrList< NodeType >::
tPtrList( const tPtrList< NodeType > & orig )
{
   tPtrListNode< NodeType > *node;
   NodeType *NTPtr;
   if( &orig != 0 )
   {
      node = orig.first;
      NTPtr = node->Ptr;
      insertAtBack( NTPtr );
      for( node = node->next; node != orig.first; node = node->next )
      {
         NTPtr = node->Ptr;
         insertAtBack( NTPtr );
      }
      if( orig.last->next == orig.first ) last->next = first;
   }
}


template< class NodeType >                      //tPtrList
tPtrList< NodeType >::
~tPtrList()
{
   first = last = 0;
     //cout << "    ~tPtrList()" << endl;
}

template< class NodeType >                      //tPtrList
const tPtrList< NodeType > &tPtrList< NodeType >::
operator=( const tPtrList< NodeType > &right )
{
   if( this != &right )
   {
      Flush();
      tPtrListNode< NodeType > *cn = right.first;
      if( cn != 0 )
      {
         insertAtBack( cn->Ptr );
         for( cn = cn->next; cn != right.last->next; cn = cn->next )
             insertAtBack( cn->Ptr );
         if( right.last->next == right.first ) last->next = first;
      }
   }
   return *this;
}

template< class NodeType >                      //tPtrList
tPtrListNode< NodeType > * tPtrList< NodeType >::
getNewNode( NodeType *NTPtr )
{
   tPtrListNode< NodeType > * newptr =
       new tPtrListNode< NodeType >( NTPtr );
   assert( newptr != 0 );
     //newptr->Ptr = NTPtr;
     //newptr->setPtr( NTPtr );
     //cout << "new ptr node created" << endl;
   nNodes++;
   return newptr;
}

//empty?
template< class NodeType >                      //tPtrList
int tPtrList< NodeType >::
isEmpty() const
{
     //cout << "checking if tPtrList empty" << endl;
   if( first == 0 )
   {
        //cout << "tPtrList is empty" << endl;
      return 1;
   }
   else
   {
        //cout << "tPtrList is not empty" << endl;
      return 0;
   }
}

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
insertAtFront( NodeType *NTPtr )
{
   tPtrListNode< NodeType > *newPtr = getNewNode( NTPtr );
   if( isEmpty() ) first = last = newPtr;
   else
   {
      newPtr->next = first;
      if( last->next == first ) last->next = newPtr;
      first = newPtr;
   }
   //nNodes++;
}

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
insertAtBack( NodeType *NTPtr )
{
   tPtrListNode< NodeType > * newPtr = getNewNode( NTPtr );
     //cout << "add new node to tPtrList" << endl;
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
     //cout << "added node to back of tPtrList" << endl;
}

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
insertAtNext( NodeType *NTPtr, tPtrListNode< NodeType > * prev )
{
   if( prev != 0 )
   {
      if( prev == last )
      {
         insertAtBack( NTPtr );
         return;
      }
      tPtrListNode< NodeType > * newPtr = getNewNode( NTPtr );
      newPtr->next = prev->next;
      prev->next = newPtr;
   }
}

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
insertAtPrev( NodeType *NTPtr, tPtrListNode< NodeType > * node )
{
   tPtrListNode< NodeType > * prev;
   if( node != 0 )
   {
      if( node == first )
      {
         insertAtFront( NTPtr );
         return;
      }
      tPtrListNode< NodeType > * newPtr = getNewNode( NTPtr );
      for( prev = first; prev->next != node; prev = prev->next ); 
      newPtr->next = prev->next;
      prev->next = newPtr;
      //nNodes++;
   }
}

template< class NodeType >                      //tPtrList
int tPtrList< NodeType >::
removeFromFront( NodeType *NTPtr )
{
   if( isEmpty() ) return 0;
   else
   {
      tPtrListNode< NodeType > * temp = first;
      if( first == last ) first = last = 0;
      else
      {
         if( last->next == first ) last->next = first->next;
         first = first->next;
      }
      NTPtr = temp->Ptr;
      delete temp;
      nNodes--;
      return 1;
   }
}

template< class NodeType >                      //tPtrList
int tPtrList< NodeType >::
removeFromBack( NodeType *NTPtr )
{
   if( isEmpty() ) return 0;
   else
   {
      tPtrListNode< NodeType > * temp = last;
      if( first == last ) first = last = 0;
      else
      {
         tPtrListNode< NodeType > * current = first;
         while( current->next != last ) current = current->next;
         current->next = last->next;
         last = current;
      }
      NTPtr = temp->Ptr;
      delete temp;
      nNodes--;
      return 1;
   }
}

template< class NodeType >                      //tPtrList
int tPtrList< NodeType >::
removeNext( NodeType *NTPtr, tPtrListNode< NodeType > * ptr )
{
   if( ptr->next == 0 ) return 0;
   if( ptr == 0 ) return 0;
   if( ptr->next == last ) return removeFromBack( NTPtr );
   else if( ptr->next == first ) return removeFromFront( NTPtr );
   tPtrListNode< NodeType > * temp = ptr->next;
   ptr->next = ptr->next->next;
   NTPtr = temp->Ptr;
   delete temp;
   nNodes--;
   return 1;
}

template< class NodeType >                      //tPtrList
int tPtrList< NodeType >::
removePrev( NodeType *NTPtr, tPtrListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr == first && last->next == 0 ) return 0;
   if( ptr == first ) return removeFromBack( NTPtr );
   tPtrListNode< NodeType > * temp, *prev;
   for( temp = first; temp->next->next != ptr; temp = temp->next );
   prev = temp;
   temp = temp->next;
   if( temp == first ) return removeFromFront( NTPtr );
   prev->next = prev->next->next;
   NTPtr = temp->Ptr;
   delete temp;
   nNodes--;
   return 1;
}

//empty list
template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
Flush()
{
   NodeType *data;
   while( removeFromBack( data ) );
   assert( isEmpty() );
   nNodes = 0;
}

//display list contents
template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
print() const
{
   if( isEmpty() )
   {
      cout<<"The list is empty"<<endl<<endl;
      return;
   }
   tPtrListNode< NodeType > * current = first;
   cout<<"The list is: ";
   while( current != 0 )
   {
      cout<<current->Ptr->getID() <<' ';
      current = current->next;
   }
   cout<<endl<<endl;
}

/*template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
input( int size, tList< NodeType > *list )
{
   NodeType *temptr;
   int idin;
   for( int i=0; i<size; i++ )
   {
      cout<<"input node data:" << endl;
      cin >> idin;
      if( idin >= 0 ) insertAtBack( list->getIthDataPtrNC( idin ) );
      else
      {
         temptr = 0;
         insertAtBack( temptr );
      }
   }
}*/

//return size
template< class NodeType >                      //tPtrList
int tPtrList< NodeType >::
getSize() const
{
   /*int i;
   tPtrListNode< NodeType > *temp;
   for( i = 1, temp = first;
        temp != last; i++, temp = temp->next );*/
   return nNodes; //i;
}

template< class NodeType >                      //tPtrList
tPtrListNode< NodeType > * tPtrList< NodeType >::
getFirstNC() {return first;}

template< class NodeType >                      //tPtrList
const tPtrListNode< NodeType > * tPtrList< NodeType >::
getFirst() const {return first;}

template< class NodeType >                      //tPtrList
tPtrListNode< NodeType > * tPtrList< NodeType >::
getLast() const {return last;}

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
moveToBack( tPtrListNode< NodeType > * mvnode ) 
{
   if( mvnode != last )
   {  
      if( mvnode == first ) first = first->next;
      else
      {
         for( tPtrListNode< NodeType > * prev = first;
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

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
moveToFront( tPtrListNode< NodeType > * mvnode ) 
{
   if( mvnode != first )
   {
      for( tPtrListNode< NodeType > * prev = first;
           prev->next != mvnode;
           prev = prev->next );
      prev->next = mvnode->next;
      mvnode->next = first;
      first = mvnode;
      if( last == mvnode ) last = prev;
      if( last->next != 0 ) last->next = first;
   }
}

template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
makeCircular() {assert( first != 0 ); last->next = first;}

template< class NodeType >                      //tPtrList
const NodeType *tPtrList< NodeType >::
getIthPtr( int num ) const
{
   int i;
   tPtrListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getPtr();
}

template< class NodeType >                      //tPtrList
NodeType *tPtrList< NodeType >::
getIthPtrNC( int num ) const
{
   int i;
   tPtrListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getPtrNC();
}

