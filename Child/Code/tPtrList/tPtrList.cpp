/**************************************************************************\
**
**  tPtrList.cpp: Functions for classes tPtrList, tPtrListNode, and
**                tPtrListIter.
**
**  $Id: tPtrList.cpp,v 1.12 1999-01-05 22:36:23 stlancas Exp $
\**************************************************************************/

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
   if( !isEmpty() )
   {
      //cout<<"Destroying nodes ... "<<endl;
      tPtrListNode<NodeType > * current = first, * temp;
      first = 0;
      while( last != 0 )
      {
         temp = current;
         //cout<<temp->data<<endl;
         if( current != last ) current = current->next;
         else
         {
            current = 0;
            last = 0;
         }
         delete temp;
      }
   }
     //cout<<"All nodes destroyed"<<endl<<endl;
   //first = last = 0;
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
   //cout << "checking if tPtrList empty" << endl << flush;
   if( first == 0 )
   {
      //cout << "tPtrList is empty" << endl << flush;
      return 1;
   }
   else
   {
      //cout << "tPtrList is not empty" << endl << flush;
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
   assert( this!=0 );
   if( !isEmpty() )
   {
        //cout<<"Destroying nodes ... "<<endl;
      tPtrListNode<NodeType > * current = first, * temp;
      first = 0;
      while( last != 0 )
      {
         temp = current;
         //cout<<temp->data<<endl;
         if( current != last ) current = current->next;
         else
         {
            current = 0;
            last = 0;
         }
         delete temp;
      }
   }
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
   tPtrListNode< NodeType > *prev;
   if( mvnode != last )
   {  
      if( mvnode == first ) first = first->next;
      else
      {
         for( prev = first;
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
   tPtrListNode< NodeType > *prev;
   if( mvnode != first )
   {
      for( prev = first;
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


/**************************************************************************\
**
**         Utilities for class tPtrListIter
**
\**************************************************************************/
template< class NodeType >     //tPtrListIter
tPtrListIter< NodeType >::
tPtrListIter()
{
   ptrlistPtr = 0;
   curptrnode = 0;
   counter = 0;
     //cout << "tPtrListIter()" << endl;
}

template< class NodeType >    //tPtrListIter
tPtrListIter< NodeType >::
tPtrListIter( tPtrList< NodeType > &ptrlist )
{
   assert( &ptrlist != 0 );
   ptrlistPtr = &ptrlist;
   curptrnode = ptrlist.first;
   counter = 0;
     //cout << "tPtrListIter( ptrlist )" << endl;
}

template< class NodeType >    //tPtrListIter
tPtrListIter< NodeType >::
~tPtrListIter()
{
   ptrlistPtr = 0;
   curptrnode = 0;
     //cout << "    ~tPtrListIter()" << endl;
}

template< class NodeType >     //tPtrListIter
void tPtrListIter< NodeType >::
Reset( tPtrList< NodeType > &ptrlist )
{
   assert( &ptrlist != 0 );
   ptrlistPtr = &ptrlist;
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->first;
   counter = 0;
   //assert( curptrnode != 0 );
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
First()
{
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->first;
   counter = 0;
   if( curptrnode != 0 ) return 1;
   else return 0;
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
Last()
{
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->last;
   counter = -1;
   if( curptrnode != 0 ) return 1;
   else return 0;
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
Get( int num )
{
     //cout << "Get: num " << num << "; ";
   assert( ptrlistPtr != 0 );
   //if( num < 0 ) return 0;
   tPtrListNode< NodeType > *tempnodeptr;// = ptrlistPtr->first;
     //counter = 0;
   for( tempnodeptr = ptrlistPtr->first, counter = 0;
        counter <= ptrlistPtr->getSize() && tempnodeptr != 0;
        tempnodeptr = tempnodeptr->next, counter++ )
         //while( tempnodeptr != 0 )
   {
        //cout << "Get: tempnodeptr->id " << tempnodeptr->getPtr()->getID()
        //   << "; ";
        //tempnodeptr = tempnodeptr->next;
        //if( tempnodeptr != 0 ) counter++;
      if( tempnodeptr->Ptr->getID() == num ) break;
   }
     //cout << endl;
   if( tempnodeptr == 0 ) return 0;
   if( tempnodeptr->Ptr->getID() != num ) return 0;
   curptrnode = tempnodeptr;
   return 1;
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
Next()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 )
   {
      curptrnode = ptrlistPtr->first;
      counter = 0;
      if( curptrnode != 0 ) return 1;
      else return 0;
   }
   curptrnode = curptrnode->next;
   counter++;
   if( curptrnode != 0 ) return 1;
   else return 0;
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
Prev()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 )
   {
      curptrnode = ptrlistPtr->last;
      counter = -1;
      if( curptrnode != 0 ) return 1;
      else return 0;
   }
   if( curptrnode == ptrlistPtr->first )
   {
      if( ptrlistPtr->last->next == 0 ) return 0;
      else
      {
         assert( ptrlistPtr->last->next == curptrnode );
         curptrnode = ptrlistPtr->last;
         counter = -1;
         return 1;
      }
   }
   tPtrListNode< NodeType > *tempnode;
   for( tempnode = ptrlistPtr->first;
        tempnode->next->Ptr->getID() != curptrnode->Ptr->getID();
        tempnode = tempnode->next );
   curptrnode = tempnode;
   assert( curptrnode != 0 );
   counter--;
   return 1;
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
Where()
{
   if( curptrnode == 0 ) return -1;
   return curptrnode->getPtr()->getID();
}

template< class NodeType >     //tPtrListIter
NodeType *tPtrListIter< NodeType >::
DatPtr()
{
   if( curptrnode == 0 ) return 0;
   return curptrnode->Ptr;
}


template< class NodeType >     //tPtrListIter
tPtrListNode< NodeType > *tPtrListIter< NodeType >::
NodePtr()
{
   return curptrnode;
}

template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
NextIsNotFirst()
{
   assert( curptrnode != 0 );
   assert( ptrlistPtr != 0 );
   if( curptrnode->next == ptrlistPtr->first ) return 0;
   return 1;
}

template< class NodeType >        //tListIter
NodeType * tPtrListIter< NodeType >::
FirstP()
{
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->first;
   counter = 0;
   if( curptrnode != 0 ) return curptrnode->Ptr;
   else return 0;
}
   
template< class NodeType >        //tListIter
NodeType * tPtrListIter< NodeType >::
LastP()
{
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->last;
   counter = 0;
   if( curptrnode != 0 ) return curptrnode->Ptr;
   else return 0;
}
   
template< class NodeType >       //tPtrListIter
NodeType *tPtrListIter< NodeType >::
PrevP()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 )
   {
      curptrnode = ptrlistPtr->last;
      counter = -1;
      if( curptrnode != 0 ) return curptrnode->Ptr;
      else return 0;
   }
   if( curptrnode == ptrlistPtr->first )
   {
      if( ptrlistPtr->last->next == 0 ) return 0;
      else
      {
         assert( curptrnode == ptrlistPtr->last->next );
         curptrnode = ptrlistPtr->last;
         counter = -1;
         return curptrnode->Ptr;
      }
   }
   tPtrListNode< NodeType > *tempnode;
   for( tempnode = ptrlistPtr->first;
        tempnode->next->Ptr != curptrnode->Ptr;
        tempnode = tempnode->next );
   curptrnode = tempnode;
   assert( curptrnode != 0 );
   counter--;
   return curptrnode->Ptr;
}

template< class NodeType >        //tListIter
NodeType * tPtrListIter< NodeType >::
NextP()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 )
   {
      curptrnode = ptrlistPtr->first;
      counter = 0;
      if( curptrnode != 0 ) return curptrnode->Ptr;
      else return 0;
   }
   curptrnode = curptrnode->next;
   counter++;
   if( curptrnode != 0 ) return curptrnode->Ptr;
   else return 0;
}

template< class NodeType >       //tListIter
NodeType * tPtrListIter< NodeType >::
GetP( int num )
{
   assert( ptrlistPtr != 0 );
   //if( num < 0 ) return 0;
   //cout << "Get: num " << num << "; ";
   int i;
   tPtrListNode< NodeType > *tempnodeptr = ptrlistPtr->first;
   counter = 0;
   while( tempnodeptr->Ptr->getID() != num && tempnodeptr != 0 )
   {
      //cout << "Get: tempnodeptr->id " << tempnodeptr->getDataPtr()->getID()
      //     << "; ";
      tempnodeptr = tempnodeptr->next;
      assert( tempnodeptr != 0 );
      counter++;
   }
   //cout << endl;
   if( tempnodeptr == 0 ) return 0;
   if( tempnodeptr->Ptr->getID() != num ) return 0;
   curptrnode = tempnodeptr;
   return tempnodeptr->Ptr;
}

template< class NodeType >        //tListIter
NodeType * tPtrListIter< NodeType >::
ReportNextP()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 ) return 0;
   if( curptrnode->next != 0 ) return curptrnode->next->Ptr;
   else return 0;
}

template< class NodeType >        //tListIter
NodeType * tPtrListIter< NodeType >::
ReportPrevP()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 ) return 0;
   if( curptrnode == ptrlistPtr->first )
   {
      if( ptrlistPtr->last->next == 0 ) return 0;
      else
      {
         assert( ptrlistPtr->last->next == curptrnode );
         return ptrlistPtr->last->Ptr;
      }
   }
   tPtrListNode< NodeType > *tempnode;
   for( tempnode = ptrlistPtr->first;
        tempnode->next != curptrnode;
        tempnode = tempnode->next );
   assert( tempnode != 0 );
   return tempnode->Ptr;
}

template< class NodeType >       //tListIter
int tPtrListIter< NodeType >::
AtEnd()
{
   if( ptrlistPtr->last == 0 ) return 1;
   if( ptrlistPtr->last->next == 0 ) return ( curptrnode==0 );
   else return ( curptrnode == ptrlistPtr->first && counter != 0 );
}

