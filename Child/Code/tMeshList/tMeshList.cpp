/**************************************************************************\
**
**  tGridList.cpp
**
**  Functions for derived class tGridList. The class is declared in
**  tGridList.h.
**
**  Modifications:
**   - added "MoveToActiveBack()" function, 12/97 GT
**
**  $Id: tMeshList.cpp,v 1.1 1998-01-14 20:19:49 gtucker Exp $
\**************************************************************************/

#include <iostream.h>
#include <fstream.h>
#include <assert.h>

#include "../Definitions.h"
#include "../Classes.h"
#include "../tListNode/tListNode.h"
#include "../tList/tList.h"

#include "tGridList.h"


template< class NodeType >                     //tGridtList
tGridList< NodeType >::
tGridList()
{
   nActiveNodes = 0;
   lastactive = 0;
     //cout << "                  from tGridList()" << first << endl;
}

template< class NodeType >                     //tGridtList
tGridList< NodeType >::
tGridList( const tGridList< NodeType > *original )
        : tList< NodeType >( original )
{
   nActiveNodes = original->nActiveNodes;
   if( original->lastactive != 0 ) lastactive = original->lastactive;
   else lastactive = 0;
     //cout << "                  from tGridList( original )" << first << endl;
}

template< class NodeType >                     //tGridtList
tGridList< NodeType >::
~tGridList()
{
     //cout << "                  from ~tGridList()" << first << endl;
}

//overloaded assignment operator
template< class NodeType >                     //tGridList
const tGridList< NodeType > &tGridList< NodeType >::
operator=( const tGridList< NodeType > &right )
{
   if( this != &right )
   {
      tList< NodeType >::operator=( right );
      lastactive = right.lastactive;
      nActiveNodes = right.nActiveNodes;
   }
     //cout << "tGridList assigned" << first << endl;
   return *this;
}

//overloaded equality operator:
template< class NodeType >                      //tGridList
int tGridList< NodeType >::
operator==( const tGridList< NodeType > &right ) const
{
   if( tList< NodeType >::operator!=( right ) ) return 0;
   if( nActiveNodes != right.nActiveNodes ) return 0;
   if( lastactive != right.lastactive ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                      //tGridList
int tGridList< NodeType >::
operator!=( const tGridList< NodeType > &right ) const
{
   if( tList< NodeType >::operator!=( right ) ) return 1;
   if( nActiveNodes != right.nActiveNodes ) return 1;
   if( lastactive != right.lastactive ) return 1;
   return 0;
}

template< class NodeType >                      //tGridList
int tGridList< NodeType >::
getActiveSize() const {return nActiveNodes;}

template< class NodeType >                      //tGridList
tListNode< NodeType > *
tGridList< NodeType >::
getLastActive() const {return lastactive;}

template< class NodeType >                     //tGridtList
void tGridList< NodeType >::
setNActiveNodes( int val ) {nActiveNodes = ( val >= 0 ) ? val : 0;}

template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
isActiveEmpty() const
{
     //cout << "checking if tGridList empty of active nodes" << endl;
   if( lastactive == 0 )
   {
        //cout << "tGridList is empty of active nodes" << endl;
      return 1;
   }
   else
   {
        //cout << "tGridList is not empty of active nodes" << endl;
      return 0;
   }
}

template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
isBoundEmpty() const
{
   if( lastactive == last ) return 1;
   else return 0;
}

template< class NodeType >                     //tGridtList
void tGridList< NodeType >::
insertAtBoundFront( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node at boundary front of tGridList" << endl;
   assert( newPtr>0 );
   assert( this != 0 );
   
   if( isEmpty() )  // Case list empty
   {
      first = last = newPtr;
   }
   else if( lastactive==0 ) // Case active part of list empty
   {
        //cout << "Case no active nodes" << endl;
      newPtr->next = first;
      first = newPtr;
   }
   else  // Usual case: list and active part of list NOT empty
   {
      newPtr->next = lastactive->next;
      lastactive->next = newPtr;
      if( lastactive==last ) last = newPtr; // Case: new is last (only) bdy
   }
   
     //nNodes++;  // why commented out? why not increment nnodes & nactive?
     //cout << "added node to front of boundary tGridList" << endl;
    
}


template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
removeFromBoundFront( NodeType &value )
{
   assert( &value != 0 );
   if( isEmpty() ) return 0;
   else if( last == lastactive ) return 0;
   else
   {
      tListNode< NodeType > * temp = lastactive->next;
      if( first == last ) first = last = 0;
      else lastactive->next = lastactive->next->next;
      value = temp->data;
      delete temp;
      nNodes--;
      return 1;
   }
}
   

template< class NodeType >                     //tGridtList
void tGridList< NodeType >::
insertAtActiveBack( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node at active back of tGridList" << endl;
   assert( this != 0 );
   //cout << " isActiveEmpty() = " << isActiveEmpty() << endl;
   if( isEmpty() )
   {
      first = last = lastactive = newPtr;
   }
   else if( !( isEmpty() ) && isActiveEmpty() && !( isBoundEmpty() ) )
   {
      lastactive = newPtr;
      lastactive->next = first;
      first = lastactive;
      //cout << "first = lastactive: " << first->getDataRef().getID()
      //   << " " << lastactive->getDataRef().getID() << endl;
      //newPtr->next = first;
      //first = lastactive = newPtr;
   }
   else if( !( isEmpty() ) && isBoundEmpty() )
   {
      newPtr->next = lastactive->next;
      lastactive->next = newPtr;
      lastactive = newPtr;
      last = lastactive;
   }
   else
   {
      newPtr->next = lastactive->next;
      lastactive->next = newPtr;
      lastactive = newPtr;
   }
   if( isBoundEmpty() ) last = lastactive;
     //nNodes++;
   nActiveNodes++;
     //cout << "added node to back of active tGridList" << endl;
}

template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
removeFromActiveBack( NodeType &value )
{
   if( isEmpty() ) return 0;
   else
   {
      tListNode< NodeType > * temp = lastactive;
      if( first == lastactive ) lastactive = 0;
      if( first == last ) first = last = 0;
      else
      {
         tListNode< NodeType > * current = first;
         while( current->next != lastactive ) current = current->next;
         current->next = lastactive->next;
         lastactive->next = 0;
         lastactive = current;
      }
      value = temp->data;
      delete temp;
      nNodes--;
      nActiveNodes--;
      return 1;
   }
}

//delete next node
template< class NodeType >                         //tList
int tGridList< NodeType >::
removeNext( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr->next == 0 ) return 0;
   if( ptr == 0 ) return 0;
   if( ptr->next == lastactive ) return removeFromActiveBack( value );
   if( ptr == lastactive ) return removeFromBoundFront( value );
   if( tList< NodeType >::removeNext( value, ptr ) )
   {
      if( !( value.getBoundaryFlag() ) ) nActiveNodes--;
      return 1;
   }
   return 0;
}

//delete previous node
template< class NodeType >                         //tList
int tGridList< NodeType >::
removePrev( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr == first && last->next == 0 ) return 0;
   if( lastactive->next == ptr ) return removeFromActiveBack( value );
   if( tList< NodeType >::removePrev( value, ptr ) )
   {
      if( !( value.getBoundaryFlag() ) ) nActiveNodes--;
      return 1;
   }
   return 0;
}

//'move' utilities
template< class NodeType >                         //tList
void tGridList< NodeType >::
moveToBack( tListNode< NodeType > * mvnode ) 
{
   if( mvnode != last )
   {
      if( mvnode == lastactive )
      {
         for( tListNode< NodeType > * prev = first;
              prev->next != mvnode;
              prev = prev->next );
         lastactive = prev;
      }
      tList< NodeType >::moveToBack( mvnode );
   }
}

template< class NodeType >                         //tList
void tGridList< NodeType >::
moveToFront( tListNode< NodeType > * mvnode ) 
{
   if( mvnode != first )
   {
      if( mvnode == lastactive )
      {
         for( tListNode< NodeType > * prev = first;
              prev->next != mvnode;
              prev = prev->next );
         lastactive = prev;
      }
      tList< NodeType >::moveToFront( mvnode );
   }
}

template< class NodeType >                         //tList
void tGridList< NodeType >::
moveToActiveBack( tListNode< NodeType > * mvnode ) 
{
   if( mvnode != lastactive )
   {
      // Detach mvnode from its position on the list
      if( mvnode == first ) first = first->next;
      else
      {
         tListNode< NodeType > * prev = first;
         while( prev->next != mvnode ) prev = prev->next;
         prev->next = mvnode->next;
      }

      // Insert it at the end of the active part of the list
      mvnode->next = lastactive->next;
      lastactive->next = mvnode;
      if( lastactive == last )
      {
         last = mvnode;
         // If it's a circular list, make sure to preserve circularity
         if( last->next != 0 ) last->next = first;
      }
      lastactive = mvnode;
   }
}

template< class NodeType >                         //tList
void tGridList< NodeType >::
insertAtFront( const NodeType &value )
{
   tList< NodeType >::insertAtFront( value );
   if( isActiveEmpty() ) lastactive = first;
   nActiveNodes++;
}

template< class NodeType >                         //tList
int tGridList< NodeType >::
removeFromFront( NodeType &value )
{
   if( !( isActiveEmpty() ) )
   {
      nActiveNodes--;
      if( lastactive == first ) lastactive = 0;
   }
   return tList< NodeType >::removeFromFront( value );
}

template< class NodeType >                         //tList
void tGridList< NodeType >::
Flush()
{
   tList< NodeType >::Flush();
   lastactive = 0;
   nActiveNodes = 0;
}

