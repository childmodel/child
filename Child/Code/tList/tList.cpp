/**************************************************************************\
**
**  tList.cpp:  Functions for class tList and related classes tListNode
**              and tListIter. (See tList.h for a description of what
**              tLists do).
**
**  Changes:
**    - GT added currentItem, FirstP(), and NextP(), plus modifications
**      to prevent currentItem from getting corrupted (1/22/99)
**
**  $Id: tList.cpp,v 1.15 1999-01-26 20:35:34 gtucker Exp $
\**************************************************************************/

#include "tList.h"

/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tListNode< NodeType >
**
**         tListNodes contain a data item and a pointer to the next
**         tListNode in the tList (or tGridList). The data item may
**         be of any type, specified in the template brackets.
**
**         Some of the functions for retrieving the data are duplicated
**         by tListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/

/**************************************************************************\
**
**  tListNode constructors:
**
**  Default constructor: sets next to null
**  Copy constructor #1: makes a copy of a given tListNode
**  Copy constructor #2: fills in data item w/ copy of given NodeType
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


/**************************************************************************\
**
**  tListNode overloaded operators:
**
**  Assignment: makes a copy (including next ptr)
**  Equality: compares both data contents and next ptr
**  Inequality: compares both data contents and next ptr
**
\**************************************************************************/

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


/**************************************************************************\
**
**  tListNode "get" functions:
**  (note: to "set" an item, use non-const "get")
**
**  getDataNC: returns a non-const (modifiable) copy of data
**  getDataRefNC: returns a non-const (modifiable) reference to data
**  getDataPtrNC: returns a non-const (modifiable) pointer to data
**  getNextNC: returns non-const ptr to next item on list
**  getData: returns const copy of data
**  getDataRef: returns const reference to data
**  getDataPtr: returns const ptr to data
**  getNext: returns const ptr to next item on list
**
\**************************************************************************/

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
**         FUNCTIONS FOR CLASS tList< NodeType >
**
**         tList is a linked list of tListNodes. Functions to add and
**         remove items to/from list, etc.
**
**         Again, some functions are duplicated by tListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/

/**************************************************************************\
**
**  tList constructors & destructor:
**
**  Default constructor: initializes all values to 0 (empty list)
**  Copy constructor: makes a complete copy of another tList
**  Destructor: deletes all nodes on list
**
\**************************************************************************/

//default constructor
template< class NodeType >                         //tList
tList< NodeType >::tList()
{
   first = last = currentItem = 0;
   nNodes = 0;
     //cout << "list instantiated" << first << endl;
}

//copy constructor
template< class NodeType >                         //tList
tList< NodeType >::
tList( const tList< NodeType > *original )
{
   int i;

   assert( original != 0 );
   nNodes = 0;
     //nNodes = original->nNodes;
   tListNode<NodeType> * current = original->first;
   for( i=0; i<original->nNodes; i++ )
   {
      insertAtBack( current->data );
      current = current->next;
   }
   assert( nNodes == original->nNodes );
   cout << "list copy instantiated" << first << endl;
   current = first;
   
}

//destructor
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

/**************************************************************************\
**
**  tList overloaded operators:
**
**  Assignment: clears the list and makes a copy of the right-hand list
**  Equality: returns TRUE if first & last pointers are identical and
**            nNodes is the same. Note that two lists with identical
**            contents are still not considered equal! (TODO--makes sense?)
**  Inequality: opposite of equality
**
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >                         //tList
const tList< NodeType > &tList< NodeType >::
operator=( const tList< NodeType > &right )
{
   if( this != &right )
   {
      // create an equivalent object rather than pointers to the original
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


/**************************************************************************\
**
**  tList::getNewNode
**
**  Creates a new tListNode and returns a pointer to it. Used by list
**  insertion routines (see below); not publically accessible.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tList: list insertion and removal routines
**
**  A collection of routines to add nodes to the list, remove them from
**  the list, and change the position of nodes on the list. In the case
**  of insertion routines, a new copy of the given node is created first.
**  The list removal routines take a reference to a NodeType variable;
**  on return, this variable contains a copy of the data item removed.
**  List removal routines return TRUE if there was an item to remove;
**  FALSE otherwise.
**
**    insertAtFront: make a new node containing _value_ and put at top
**    insertAtBack: make a new node containing _value_ and put at bottom
**    insertAtNext: make new node containing _value_ and place it on the
**                  list after _prev_
**    insertAtPrev: make new node containing _value_ and place it on the
**                  list before _node_
**    removeFromFront: remove 1st item on list and place a copy in _value_
**    removeFromBack: remove last item on list and place a copy in _value_
**    removeNext: remove the node following node _ptr_ and place a copy
**                in _value_
**    removePrev: remove the node before node _ptr_ and place a copy
**                in _value_
**
\**************************************************************************/

//insert at front
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtFront( const NodeType &value )
{
   //cout << "ADD NEW NODE TO LIST AT FRONT" << endl;
   
   tListNode< NodeType > *newPtr = getNewNode( value );
   if( isEmpty() ) first = last = currentItem = newPtr;
   else
   {
      newPtr->next = first;
      if( last->next == first ) last->next = newPtr;
      first = newPtr;
   }
     //nNodes++; now handled by getNewNode
}

//insert at back
template< class NodeType >                         //tList
void tList< NodeType >::
insertAtBack( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
   //cout << "add new node to list in back" << endl;
   assert( this != 0 );
   if( isEmpty() )
   {
      first = last = currentItem = newPtr;
   }
   else
   {
      newPtr->next = last->next;
      last->next = newPtr;
      last = newPtr;
   }

   //cout << "data is " << newPtr->data << endl;
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
      if( first == last ) first = last = currentItem = 0;
      else
      {
         if( last->next == first ) last->next = first->next;
         if( currentItem==first ) currentItem = first->next;
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
      if( first == last ) first = last = currentItem = 0;
      else
      {
         tListNode< NodeType > * nexttolast = first;
         while( nexttolast->next != last ) nexttolast = nexttolast->next;
         nexttolast->next = last->next;
         if( currentItem==last ) currentItem = nexttolast;
         last = nexttolast;
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
   if( currentItem == temp ) currentItem = ptr;
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
   if( currentItem == temp ) currentItem = prev;
   delete temp;
   nNodes--;
   return 1;
}


/**************************************************************************\
**
**  tList::Flush
**
**  Deletes all nodes on list.
**
\**************************************************************************/
template< class NodeType >                         //tList
void tList< NodeType >::
Flush()
{
   tListNode<NodeType > * current = first, * temp;

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
      first = last = currentItem = 0;
   }
   assert( isEmpty() );
   nNodes = 0;
     //cout<<"All nodes destroyed"<<endl<<endl;
}


/**************************************************************************\
**
**  tList::isEmpty
**
**  Returns TRUE if first points to null; FALSE otherwise.
**
\**************************************************************************/
//empty?
template< class NodeType >                         //tList
int tList< NodeType >::
isEmpty() const
{
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

//display list contents -- for debugging only
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


/**************************************************************************\
**
**  tList "get" functions:
**
**  getSize: returns # of items on list
**  getFirst: returns const ptr to first tListNode
**  getLast: returns const ptr to last tListNode
**  FirstP: returns ptr to first data item and sets currentItem to first
**  NextP: returns ptr to data item following currentItem and advances
**         currentItem to the next one on the list
**  getIthData: returns const copy of list data item number _num_
**  getIthDataRef: returns const reference to list data item number _num_
**  getIthDataRef: returns const ptr to list data item number _num_
**  getIthDataNC: returns non-const copy of list data item number _num_
**  getIthDataRefNC: returns non-const ref to list data item number _num_
**  getIthDataRefNC: returns non-const ptr to list data item number _num_
**  (see also getListNode, below)
**
\**************************************************************************/

template< class NodeType >                         //tList
int tList< NodeType >::
getSize() const {return nNodes;}

template< class NodeType >                         //tList
tListNode< NodeType > * tList< NodeType >::
getFirst() const {return first;}

template< class NodeType >                         //tList
tListNode< NodeType > * tList< NodeType >::
getLast() const {return last;}

// Added by gt 1/99
template< class NodeType >                         //tList
NodeType * tList< NodeType >::
FirstP() 
{
   assert( first!=0 );
   currentItem = first;
   return &first->data;
}

// Added by gt 1/99
template< class NodeType >                         //tList
NodeType * tList< NodeType >::
NextP() 
{
   assert( currentItem!=0 );
   currentItem = currentItem->next;
   if( currentItem!=0 )
       return &currentItem->data;
   else return 0;
}

template< class NodeType >                         //tList
const NodeType tList< NodeType >::
getIthData( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
//    if(num>= nNodes)
//        {
//           cout<<"using an index which is too large"<<endl;
//           cout<<"you have "<<nNodes<<endl;
//           cout<<"you wanted list member number "<<num<<endl;
//        }
//    if(num<0)
//    {
//       cout<<"using a negative index"<<endl;
//       cout<<"you have "<<nNodes<<endl;
//       cout<<"you wanted list member number "<<num<<endl;
//    }
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
**  tList "move" functions:
**
**  moveToBack: looks for _mvnode_ on the list and moves it to the back
**              of the list. Assumes that _mvnode_ IS an item on the list.
**  moveToFront: same thing, but moves it to the front
**  
\**************************************************************************/

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



/* GT commented this out 1/99 -- doesn't seem to be used anywhere, and
   is a bit dangerous bcs the actual # of items on list is untouched
template< class NodeType >                         //tList
void tList< NodeType >::
setNNodes( int val ) {nNodes = ( val >= 0 ) ? val : 0;}*/

/**************************************************************************\
**
**  tList::makeCircular
**
**  Converts the list into a circular list by having the last item point
**  to the first.
**
\**************************************************************************/
template< class NodeType >                         //tList
void tList< NodeType >::
makeCircular() {last->next = first;}




/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tListIter
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
**
\**************************************************************************/

/**************************************************************************\
**
**  tListIter constructors & destructor:
**
**  Default constructor: initializes all values to 0
**  Constructor (reference version): sets listPtr to point to _list_ and
**                                   points curnode to 1st node on list
**  Constructor (pointer version): same thing, but takes a pointer to a
**                                 tList as an argument
**
\**************************************************************************/

template< class NodeType >        //tListIter
tListIter< NodeType >::
tListIter()
{
   listPtr = 0;
   curnode = 0;
   counter = 0;
   assert( this != 0 );
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


/**************************************************************************\
**
**  tListIter::Reset
**
**  Points iterator at the 1st node on _list_ (provides a way of telling
**  an iterator which list to work on).
**
\**************************************************************************/
template< class NodeType >        //tListIter
void tListIter< NodeType >::
Reset( tList< NodeType > &list )
{
   assert( &list != 0 );
   listPtr = &list;
   curnode = list.first;
   counter = 0;
}


/**************************************************************************\
**
**  tListIter::First and tListIter::Last
**
**  Move to the first or last item on the current list. Return TRUE if
**  pointing to a valid list item, FALSE otherwise.
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tListIter::Get
**
**  Move to list item with ID number _num_. Note: assumes that list items
**  have a member function getID()! Returns 1 if found, 0 if not.
**
\**************************************************************************/
template< class NodeType >     //tListIter
int tListIter< NodeType >::
Get( int num )
{
   assert( listPtr != 0 );
   if( num < 0 )
   {
      cout << "tListIter::Get(num): num < 0" << endl;
      //return 0;
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

   
/**************************************************************************\
**
**  tListIter::Next and tListIter::Prev
**
**  Move to the next or previous item on the current list. Return TRUE if
**  pointing to a valid list item, FALSE otherwise. If we're not 
**  initially pointing to any item, then move to the first or last item,
**  respectively. Both assume we're working on a valid list.
**
\**************************************************************************/
template< class NodeType >        //tListIter
int tListIter< NodeType >::
Next()
{
   assert( listPtr != 0 );

   // if current position undefined, move to first node...
   if( curnode == 0 )
   {
      curnode = listPtr->first;
      counter = 0;
      if( curnode != 0 ) return 1;
      else return 0;
   }

   // ...otherwise just move to the next one
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

   // if current position undefined, move to the last one
   if( curnode == 0 )
   {
      curnode = listPtr->last;
      counter = -1; // why -1 and not nNodes? TODO
      if( curnode != 0 ) return 1;
      else return 0;
   }

   // if we're at the first node, the previous one is only defined if we're
   // a circular list, in which case last points to first -- so move to last
   if( curnode == listPtr->first )
   {
      if( listPtr->last->next == 0 ) return 0;
      else
      {
         assert( curnode == listPtr->last->next );
         curnode = listPtr->last;
         counter = -1; // why -1?
         return 1;
      }
   }

   // general case: search through the list until we reach the one before
   // curnode, then set curnode to that one
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


/**************************************************************************\
**
**  tListIter::FirstP and tListIter::LastP
**
**  Move to the first or last item on the list and return a pointer to the
**  data, or 0 if first/last item is empty.
**
\**************************************************************************/
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
   

/**************************************************************************\
**
**  tListIter::NextP and tListIter::PrevP
**
**  Same as Next and Prev, except that the functions return a pointer to
**  the current data item (or 0 if none exists).
**
\**************************************************************************/

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


/**************************************************************************\
**
**  tListIter::GetP
**
**  Similar to Get, but returns a pointer to the current data item (or
**  0 if undefined).
**
\**************************************************************************/
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


/**************************************************************************\
**
**  tListIter::Where
**
**  Returns the ID number of the current data item, or -1 if there is
**  no current data item. Assumes data item has a getID() mbr function!
**
\**************************************************************************/
template< class NodeType >       //tListIter
int tListIter< NodeType >::
Where()
{
   if( curnode == 0 ) return -1;
   return curnode->getDataPtr()->getID();
}


/**************************************************************************\
**
**  tListIter::AtEnd
**
**  Returns TRUE if:
**   - the list is empty
**   - the list is non-circular and the current item is null
**   - the list is circular, the current item is the first, and the 
**     counter is nonzero (meaning we've gone all the way through the
**     list and come back to the start)
**
\**************************************************************************/
template< class NodeType >       //tListIter
int tListIter< NodeType >::
AtEnd()
{
   if( listPtr->isEmpty() ) return 1;
   if( listPtr->last->next == 0 ) return ( curnode==0 );
   else return ( curnode == listPtr->first && counter != 0 );
   //return curnode==0;
}


/**************************************************************************\
**
**  tListIter::DatRef, DatPtr, and NodePtr
**
**  Returns a non-constant reference or pointer to the current data item
**  or the current list node.
**
\**************************************************************************/
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



