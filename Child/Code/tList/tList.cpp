/**************************************************************************\
**
**  tList.cpp:  Functions for class tList and related classes tListNode
**              and tListIter.
**
**  $Id: tList.cpp,v 1.9 1998-06-10 19:18:32 nmgaspar Exp $
\**************************************************************************/

#include "tList.h"

/**************************************************************************\
**
**         Utilities for class tListNode< NodeType >
**
**         tListNodes contain a data item and a pointer to the next
**         tListNode in the tList (or tGridList). The data item may
**         be of any type, specified in the template brackets.
**
**         Some of the functions for retrieving the data are made
**         obsolete by tListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/
//default constructor
template< class NodeType >                     //tListNode
tListNode< NodeType >::
tListNode()
{
     //data = 0;
   next = 0;
}

//copy constructor with data reference
template< class NodeType >                     //tListNode
tListNode< NodeType >::
tListNode( const tListNode< NodeType > &original )
{
   if( &original != 0 )
   {
      data = original.data;
      next = original.next;
   }
}

//value (by reference) constructor 
template< class NodeType >                     //tListNode
tListNode< NodeType >::
tListNode( const NodeType &info )
{
   data = info;
   next = 0;
}

//overloaded assignment operator
template< class NodeType >                     //tListNode
const tListNode< NodeType > &tListNode< NodeType >::
operator=( const tListNode< NodeType > &right )
{
   if( &right != this )
   {
        //delete data;
        //data = new NodeType;
      assert( &data != 0 );
      data = right.data;
      next = right.next;
   }
   return *this;
}

//overloaded equality operator:
template< class NodeType >                     //tListNode
int tListNode< NodeType >::
operator==( const tListNode< NodeType > &right ) const
{
   if( next != right.next ) return 0;
   if( &data != &(right.data) ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                     //tListNode
int tListNode< NodeType >::
operator!=( const tListNode< NodeType > &right ) const
{
   if( next != right.next ) return 1;
   if( &data != &(right.data) ) return 1;
   return 0;
}

//set data by returning non-const
template< class NodeType >                     //tListNode
NodeType tListNode< NodeType >::
getDataNC() {return data;}

template< class NodeType >                     //tListNode
NodeType &tListNode< NodeType >::
getDataRefNC() {return data;}

template< class NodeType >                     //tListNode
NodeType *tListNode< NodeType >::
getDataPtrNC() {return &data;}

template< class NodeType >                     //tListNode
tListNode< NodeType > * tListNode< NodeType >::
getNextNC() const {return next;}

//return data by value
template< class NodeType >                     //tListNode
NodeType tListNode< NodeType >::
getData() const {return data;}

//return data by reference
template< class NodeType >                     //tListNode
const NodeType &tListNode< NodeType >::
getDataRef() const {return data;}

//return data by pointer
template< class NodeType >                     //tListNode
const NodeType *tListNode< NodeType >::
getDataPtr() const {return &data;}

//return next pointer
template< class NodeType >                     //tListNode
const tListNode< NodeType > * tListNode< NodeType >::
getNext() const {return next;}


/**************************************************************************\
**
**  Functions for class tList< NodeType >
**
**   tList is a linked list of tListNodes. Functions to add and remove
**   items from list, etc.
**
**   Again, some functions are made obsolete by tListIter.
**
**   Created: fall, 97, SL
**
\**************************************************************************/

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
   int i;
   assert( original != 0 );
   nNodes = 0;
     //nNodes = original->nNodes;
   tListNode< NodeType > * current = original->first;
   for( i=0; i<original->nNodes; i++ )
   {
      insertAtBack( current->data );
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
      last = right.last;*/
      assert( nNodes == right.nNodes );
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
   //Xcout << "new node created" << endl;
   nNodes++;
   return ptr;
}

//insert at front
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtFront( const NodeType &value )
{
   //cout << "ADD NEW NODE TO LIST AT FRONT" << endl;
   
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
   //  cout << "add new node to list in back" << endl;
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
   tListNode<NodeType > * current = first,* temp;
   //while( !isEmpty() ) removeFromBack( temp );
   
   if( !isEmpty() )
   {
      //cout<<"Destroying nodes ... "<<endl;
      //while( removeFromBack( temp ) );
      while( current != 0 )
      {
         temp = current;
         //cout<<temp->data<<endl;
         current = current->next;
         delete temp;
      }
      first = last = 0;
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
      //cout<< current->data <<' ';
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
   tListNode< NodeType > * prev;
   if( mvnode != last )
   {  
      if( mvnode == first ) first = first->next;
      else
      {
         for( prev = first; prev->next != mvnode; prev = prev->next );
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
   tListNode< NodeType > * prev;
   if( mvnode != first )
   {
      for( prev = first; prev->next != mvnode; prev = prev->next );
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


/**************************************************************************\
**
**  tList::getListNode
**
**  Finds and returns the list node containing the data pointed to by
**  desiredDatPtr, or zero if not found.
**
**  Parameters:  desiredDatPtr -- pointer to the data item sought after
**  Returns:  pointer to the tListNode containing desiredDatPtr, or zero
**            if not found
**  Notes: might be safer to implement with a const return type
**  Created: 4/29/98 GT
**
\**************************************************************************/
template< class NodeType >
tListNode< NodeType > * tList< NodeType >::
getListNode( NodeType * desiredDatPtr )
{
   tListNode< NodeType > * listnode = first;

   if( listnode==0 ) return 0;
   while( &(listnode->data) != desiredDatPtr )
   {
      listnode = listnode->next;
      if( listnode==0 ) return 0;
   }
   return listnode;
   
}


/**************************************************************************\
**
**  CODE FOR tListIter OBJECTS.
**
**  A tListIter is an iterator for the linked list tList objects (and their
**  descendants). Its services include fetching data from the current entry
**  on the list, advancing to the next or previous item on the list, etc.
**
**  See also tGridList.
**
**  Created: fall, 97, SL.
**  Modifications:  added an "AtEnd" function that signals whether the
**  the iterator has fallen off the end of the list, 11/17/97 gt.
*/

/**************************************************************************\
**
**         Utilities for class tListIter
**
\**************************************************************************/
template< class NodeType >        //tListIter
tListIter< NodeType >::
tListIter()
{
   listPtr = 0;
   curnode = 0;
   counter = 0;
   assert( &this != 0 );
     //cout << "tListIter()" << endl;
}

template< class NodeType >        //tListIter
tListIter< NodeType >::
tListIter( tList< NodeType > &list )
{
   assert( &list != 0 );
   listPtr = &list;
   curnode = list.first;
   counter = 0;
   //assert( curnode != 0 );
     //cout << "tListIter( list )" << endl;
   
}

template< class NodeType >        //tListIter
tListIter< NodeType >::
tListIter( tList< NodeType > *ptr )
{
   assert( ptr != 0 );
   listPtr = ptr;
   curnode = ptr->first;
   //assert( curnode != 0 );
     //cout << "tListIter( ptr )" << endl;
   
}

template< class NodeType >        //tListIter
tListIter< NodeType >::
~tListIter()
{
   listPtr = 0;
   curnode = 0;
     //cout << "~tListIter()" << endl;
}

template< class NodeType >        //tListIter
void tListIter< NodeType >::
Reset( tList< NodeType > &list )
{
   assert( &list != 0 );
   listPtr = &list;
   curnode = list.first;
   counter = 0;
}


template< class NodeType >        //tListIter
int tListIter< NodeType >::
First()
{
   assert( listPtr != 0 );
   curnode = listPtr->first;
   counter = 0;
   if( curnode != 0 ) return 1;
   else if( curnode == 0 && listPtr->isEmpty() ) return 1;
   return 0;
}

template< class NodeType >       //tListIter
int tListIter< NodeType >::
Last()
{
   assert( listPtr != 0 );
   curnode = listPtr->last;
   counter = -1;
   if( curnode != 0 ) return 1;
     //else if( curnode == 0 && listPtr->isEmpty() ) return 1;
   return 0;
}

template< class NodeType >     //tListIter
int tListIter< NodeType >::
Get( int num )
{
   assert( listPtr != 0 );
   if( num < 0 )
   {
      cout << "tListIter::Get(num): num < 0" << endl;
      return 0;
   }
   tListNode< NodeType > *tempnodeptr;
   for( tempnodeptr = listPtr->first, counter = 0;
        counter <= listPtr->nNodes && tempnodeptr != 0;
        tempnodeptr = tempnodeptr->next, counter++ )
   {
      if( tempnodeptr->data.getID() == num ) break;
   }
   if( tempnodeptr == 0 )
   {
      cout << "tListIter::Get(num): tempnodeptr == 0" << endl;
      return 0;
   }
   if( tempnodeptr->data.getID() != num )
   {
      cout << "tListIter::Get(num): tempnodeptr->data.getID() != num" << endl;
      return 0;
   }
   curnode = tempnodeptr;
   return 1;
}

//template< class NodeType >       //tListIter
/*int tListIter< NodeType >::
Get( int num )
{
   assert( listPtr != 0 );
   if( num < 0 ) return 0;
   int i;
   tListNode< NodeType > *tempnodeptr = listPtr->first;
   i = 0;
   while( tempnodeptr->getDataPtr()->getID() != num && tempnodeptr != 0 )
   {
      tempnodeptr = tempnodeptr->next;
      assert( tempnodeptr != 0 );
      i++;
   }
   if( tempnodeptr == 0 ) return 0;
   if( tempnodeptr->getDataPtr()->getID() != num ) return 0;
   curnode = tempnodeptr;
   return 1;
}*/
   
template< class NodeType >        //tListIter
int tListIter< NodeType >::
Next()
{
   assert( listPtr != 0 );
   if( curnode == 0 )
   {
      curnode = listPtr->first;
      counter = 0;
      if( curnode != 0 ) return 1;
      else return 0;
   }
   curnode = curnode->next;
   counter++;
   if( curnode != 0 ) return 1;
   else return 0;
}

template< class NodeType >       //tListIter
int tListIter< NodeType >::
Prev()
{
   assert( listPtr != 0 );
   if( curnode == 0 )
   {
      curnode = listPtr->last;
      counter = -1;
      if( curnode != 0 ) return 1;
      else return 0;
   }
   if( curnode == listPtr->first )
   {
      if( listPtr->last->next == 0 ) return 0;
      else
      {
         assert( curnode == listPtr->last->next );
         curnode = listPtr->last;
         counter = -1;
         return 1;
      }
   }
   tListNode< NodeType > *tempnode;
   //int id = curnode->data.getID();
   for( tempnode = listPtr->first;
        tempnode->next != curnode; //tempnode->next->data.getID() != id;
        tempnode = tempnode->next );
   curnode = tempnode;
   counter--;
   assert( curnode != 0 );
   return 1;
}

template< class NodeType >        //tListIter
NodeType * tListIter< NodeType >::
FirstP()
{
   assert( listPtr != 0 );
   curnode = listPtr->first;
   counter = 0;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}
   
template< class NodeType >        //tListIter
NodeType * tListIter< NodeType >::
LastP()
{
   assert( listPtr != 0 );
   curnode = listPtr->last;
   counter = 0;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}
   
template< class NodeType >        //tListIter
NodeType * tListIter< NodeType >::
NextP()
{
   assert( listPtr != 0 );
   if( curnode == 0 )
   {
      curnode = listPtr->first;
      counter = 0;
      if( curnode != 0 ) return curnode->getDataPtrNC();
      else return 0;
   }
   curnode = curnode->next;
   counter++;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}

template< class NodeType >       //tListIter
NodeType *tListIter< NodeType >::
PrevP()
{
   assert( listPtr != 0 );
   if( curnode == 0 )
   {
      curnode = listPtr->last;
      counter = -1;
      if( curnode != 0 ) return curnode->getDataPtrNC();
      else return 0;
   }
   if( curnode == listPtr->first )
   {
      if( listPtr->last->next == 0 ) return 0;
      else
      {
         assert( curnode == listPtr->last->next );
         curnode = listPtr->last;
         counter = -1;
         return curnode->getDataPtrNC();
      }
   }
   tListNode< NodeType > *tempnode;
   //int id = curnode->data.getID();
   for( tempnode = listPtr->first;
        tempnode->next != curnode;  //tempnode->next->data.getID() != id;
        tempnode = tempnode->next );
   curnode = tempnode;
   assert( curnode != 0 );
   counter--;
   return curnode->getDataPtrNC();
}

template< class NodeType >       //tListIter
NodeType * tListIter< NodeType >::
GetP( int num )
{
   assert( listPtr != 0 );
   if( num < 0 ) return 0;
   //cout << "Get: num " << num << "; ";
   int i;
   tListNode< NodeType > *tempnodeptr = listPtr->first;
   counter = 0;
   while( tempnodeptr->getDataPtr()->getID() != num && tempnodeptr != 0 )
   {
      //cout << "Get: tempnodeptr->id " << tempnodeptr->getDataPtr()->getID()
      //     << "; ";
      tempnodeptr = tempnodeptr->next;
      assert( tempnodeptr != 0 );
      counter++;
   }
   //cout << endl;
   if( tempnodeptr == 0 ) return 0;
   if( tempnodeptr->getDataPtr()->getID() != num ) return 0;
   curnode = tempnodeptr;
   return tempnodeptr->getDataPtrNC();
}

template< class NodeType >       //tListIter
int tListIter< NodeType >::
Where()
{
   if( curnode == 0 ) return -1;
   return curnode->getDataPtr()->getID();
}

template< class NodeType >       //tListIter
int tListIter< NodeType >::
AtEnd()
{
   if( listPtr->isEmpty() ) return 1;
   if( listPtr->last->next == 0 ) return ( curnode==0 );
   else return ( curnode == listPtr->first && counter != 0 );
   //return curnode==0;
}

template< class NodeType >       //tListIter
NodeType &tListIter< NodeType >::
DatRef()
{
     //if( curnode == 0 ) return 0;
   return curnode->getDataRefNC();
}

template< class NodeType >       //tListIter
NodeType *tListIter< NodeType >::
DatPtr()
{
   if( curnode == 0 ) return 0;
   return curnode->getDataPtrNC();
}

template< class NodeType >       //tListIter
tListNode< NodeType > *tListIter< NodeType >::
NodePtr()
{
   //assert( curnode != 0 );
   return curnode;
}



