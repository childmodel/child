/**************************************************************************\
**
**  tPtrList.cpp: Functions for classes tPtrList, tPtrListNode, and
**                tPtrListIter. (see tPtrList.h for a description of
**                what these classes do)
**
**  $Id: tPtrList.cpp,v 1.17 2000-01-13 19:21:07 gtucker Exp $
\**************************************************************************/

#include "tPtrList.h"


/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tPtrListNode< NodeType >
**
**         tPtrListNodes contain a data pointer and a pointer to the next
**         tPtrListNode in the tPtrList (or descendent). The data pointer
**         may be of any type, specified in the template brackets.
**
**         Some of the functions for retrieving the data are duplicated
**         by tPtrListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/

/**************************************************************************\
**
**  tPtrListNode constructors & destructor:
**
**  Default constructor: sets ptr & next to null
**  Copy constructor #1: makes a copy of a given tPtrListNode
**  Copy constructor #2: init's data ptr with NTPtr
**  Destructor: sets both pointers to zero
**
\**************************************************************************/

//default constructor
template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::tPtrListNode()
{
   Ptr = 0;
   next = 0;
}

//copy constructor with tPtrListNode
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

// copy constructor with data ptr
template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::
tPtrListNode( NodeType * NTPtr )
{
   Ptr = NTPtr;
   next = 0;
}

//destructor
template< class NodeType >                  //tPtrListNode
tPtrListNode< NodeType >::
~tPtrListNode()
{
   Ptr = 0;  //redirect the pointer away from its target
   next = 0;
}


/**************************************************************************\
**
**  tPtrListNode overloaded operators:
**
**  Assignment: makes a copy (including next ptr)
**  Equality: compares both data and next ptr (note: must point to the
**            same data location; identical data in different locations
**            are not considered equal)
**  Inequality: compares both data and next ptr (as above)
**
\**************************************************************************/

//overloaded assignment operator
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


/**************************************************************************\
**
**  tPtrListNode "get" functions:
**  (note: to "set" an item, use non-const "get")
**
**  getPtr: returns copy of data pointer as const
**  getPtrNC: returns a non-const (modifiable) copy of data ptr
**  getNext: returns const ptr to next item on list
**  getNextNC: returns non-const ptr to next item on list
**
\**************************************************************************/

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



/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tPtrList< NodeType >
**
**         tPtrList is a linked list of tPtrListNodes. Functions to add &
**         remove items to/from list, etc.
**
**         Again, some functions are duplicated by tPtrListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/

/**************************************************************************\
**
**  tPtrList constructors & destructor:
**
**  Default constructor: initializes all values to 0 (empty list)
**  Copy constructor: makes a complete copy of another tPtrList
**  Destructor: deletes all nodes on list. NOTE: does not destroy the
**              data items themselves!
**
\**************************************************************************/

//default constructor
template< class NodeType >                      //tPtrList
tPtrList< NodeType >::
tPtrList()
{
   first = last = 0;
   nNodes = 0;
     //cout << "tPtrList() instantiated" << endl;
}

//copy constructor
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

//destructor
template< class NodeType >                      //tPtrList
tPtrList< NodeType >::
~tPtrList()
{
   if( !isEmpty() )
   {
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
}


/**************************************************************************\
**
**  tPtrList overloaded operators:
**
**  Assignment: clears the list and makes a copy of the right-hand list
**
\**************************************************************************/

//overloaded assignment operator
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


/**************************************************************************\
**
**  tPtrList::getNewNode
**
**  Creates a new tPtrListNode and returns a pointer to it. Used by list
**  insertion routines (see below); not publically accessible.
**
\**************************************************************************/
template< class NodeType >                      //tPtrList
tPtrListNode< NodeType > * tPtrList< NodeType >::
getNewNode( NodeType *NTPtr )
{
   tPtrListNode< NodeType > * newptr =
       new tPtrListNode< NodeType >( NTPtr );
   assert( newptr != 0 );
   nNodes++;
   return newptr;
}


/**************************************************************************\
**
**  tPtrList: list insertion routines
**
**  A collection of routines to add items to the list. 
**
**    insertAtFront: new item with given ptr at top of list
**    insertAtBack: new item with given ptr at bottom of list
**    insertAtNext: place new node on the list after _prev_
**    insertAtPrev: place new node on the list before _node_
**
\**************************************************************************/

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


/**************************************************************************\
**
**  tPtrList::removeFromFront
**
**  Removes the first item on the list and points NTPtr to the new first
**  item. Returns 0 if the list is already empty, 1 otherwise. Note that
**  if the list is empty, NTPtr is unchanged.
**
**  ALERT: There is a potential bug here: if the list is circular but
**  contains only one item (which points to itself), NTPtr will contain
**  a dangling pointer! TODO (gt)
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrList::removeFromBack
**
**  Removes the last item on the list and points NTPtr to the new last
**  item. Returns 0 if the list is already empty, 1 otherwise. Note that
**  if the list is empty, NTPtr is unchanged.
**
**  ALERT: There is a potential bug here: if the list is circular but
**  contains only one item (which points to itself), NTPtr will contain
**  a dangling pointer! TODO (gt)
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrList::removeNext
**
**  Removes the item after _ptr_ on the list, returning the pointer in
**  NTPtr.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrList::removePrev
**
**  Removes the item before _ptr_ on the list, returning the pointer in
**  NTPtr.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrList::moveToFront and moveToBack
**
**  Moves the list item pointed to by mvnode to the front or back of
**  the list, respectively.
**
\**************************************************************************/

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


/**************************************************************************\
**
**  tPtrList::Flush
**
**  Deletes all nodes on list. NOTE: destroys only the pointers, not the
**  data.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrList::isEmpty
**
**  Returns TRUE if first points to null; FALSE otherwise.
**
\**************************************************************************/
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
   do
   {
      cout<<current->Ptr->getID() <<' ';
      current = current->next;
   }while( current != first && current != last );
   cout<<current->Ptr->getID() <<endl;
}


/**************************************************************************\
**
**  tPtrList "get" functions:
**
**  getSize: returns # of items on list
**  getFirstNC: returns non-const ptr to first tPtrListNode
**  getFirst: returns const ptr to first tPtrListNode
**  getLast: returns const ptr to last tPtrListNode
**  getIthPtr: returns a const copy of the Ith data pointer
**  getIthPtrNC: returns a non-const copy of the Ith data pointer
**  (see also getListNode)
**
\**************************************************************************/

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
**  tPtrList::makeCircular
**
**  Converts the list into a circular list by having the last item point
**  to the first.
**
\**************************************************************************/
template< class NodeType >                      //tPtrList
void tPtrList< NodeType >::
makeCircular() {assert( first != 0 ); last->next = first;}


/**************************************************************************\
**
**  tPtrList::DataCopy
**
**  Creates and returns a complete copy of the list, including copies of
**  the items that are pointed to.
**
**  Created 12/7/99, GT
**  Assumes: NodeType has a copy constructor defined
**
\**************************************************************************/
template< class NodeType >                      //tPtrList
tPtrList< NodeType > *tPtrList<NodeType>::
DataCopy()
{
   tPtrListIter iter( this );
   NodeType * curr, * newitem;
   tPtrList * newlist = new tPtrList();
   
   for( curr=iter.FirstP(); !(iter.AtEnd()); curr=iter.NextP() )
   {
      newitem = new NodeType( *curr );
      newlist.insertAtBack( newitem );
   }
   return newlist;
}



/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tPtrListIter
**
**  A tPtrListIter is an iterator for the linked list tPtrList objects (and
**  descendants). Its services include fetching data from the current entry
**  on the list, advancing to the next or previous item on the list, etc.
**
**  Created: fall, 97, SL.
**
\**************************************************************************/

/**************************************************************************\
**
**  tPtrListIter constructors & destructor:
**
**  Default constructor: initializes all values to 0
**  Constructor: attaches the iterator to _ptrlist_ and moves to the first
**               item
**  Destructor: resets values to zero (needed?)
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

/**************************************************************************\
**
**  tPtrListIter::Reset
**
**  Points iterator at the 1st node on _ptrlist_ (provides a way of telling
**  an iterator which list to work on).
**
\**************************************************************************/
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

/**************************************************************************\
**
**  tPtrListIter::First and Last
**
**  Move to the first or last item on the current list. Return TRUE if
**  pointing to a valid list item, FALSE otherwise.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrListIter::Get
**
**  Move to list item with ID number _num_. Note: assumes that list items
**  have a member function getID()! Returns 1 if found, 0 if not.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrListIter::Next and Prev
**
**  Move to the next or previous item on the current list. Return TRUE if
**  pointing to a valid list item, FALSE otherwise. If we're not 
**  initially pointing to any item, then move to the first or last item,
**  respectively. Both assume we're working on a valid list.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrListIter::Where
**
**  Returns the ID number of the current data item, or -1 if there is
**  no current data item. Assumes data item has a getID() mbr function!
**
\**************************************************************************/
template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
Where()
{
   if( curptrnode == 0 ) return -1;
   return curptrnode->getPtr()->getID();
}


/**************************************************************************\
**
**  tPtrListIter::DatPtr
**
**  Returns copy of current data pointer.
**
\**************************************************************************/
template< class NodeType >     //tPtrListIter
NodeType *tPtrListIter< NodeType >::
DatPtr()
{
   if( curptrnode == 0 ) return 0;
   return curptrnode->Ptr;
}


/**************************************************************************\
**
**  tPtrListIter::NodePtr
**
**  Returns pointer to current list node.
**
\**************************************************************************/
template< class NodeType >     //tPtrListIter
tPtrListNode< NodeType > *tPtrListIter< NodeType >::
NodePtr()
{
   return curptrnode;
}


/**************************************************************************\
**
**  tPtrListIter::NextIsNotFirst
**
**  Tests whether we're at the end of a circular list by checking whether
**  the next item is the first item (which might be true if the list is
**  circular). Returns 0 if the next item is the first on the list.
**
\**************************************************************************/
template< class NodeType >     //tPtrListIter
int tPtrListIter< NodeType >::
NextIsNotFirst()
{
   assert( curptrnode != 0 );
   assert( ptrlistPtr != 0 );
   if( curptrnode->next == ptrlistPtr->first ) return 0;
   return 1;
}


/**************************************************************************\
**
**  tPtrListIter::FirstP and tPtrListIter::LastP
**
**  Move to the first or last item on the list and return a copy of the
**  data pointer, or 0 if first/last item is empty.
**
\**************************************************************************/
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

  
/**************************************************************************\
**
**  tPtrListIter::NextP and tPtrListIter::PrevP
**
**  Same as Next and Prev, except that the functions return a copy of the
**  data pointer rather than a pointer to the list item.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrListIter::GetP
**
**  Similar to Get, but returns a copy of the current data pointer rather
**  than a pointer to the list item (or 0 if undefined).
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrListIter::ReportNextP and ReportPrevP
**
**  Returns a copy of the next or previous data pointer without actually
**  moving to the next or previous item.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tPtrListIter::AtEnd
**
**  Returns TRUE if:
**   - the list is empty (last==0)
**   - the list is non-circular and the current item is null
**   - the list is circular, the current item is the first, and the 
**     counter is nonzero (meaning we've gone all the way through the
**     list and come back to the start)
**
\**************************************************************************/
template< class NodeType >       //tListIter
int tPtrListIter< NodeType >::
AtEnd()
{
   if( ptrlistPtr->last == 0 ) return 1;
   if( ptrlistPtr->last->next == 0 ) return ( curptrnode==0 );
   else return ( curptrnode == ptrlistPtr->first && counter != 0 );
}

