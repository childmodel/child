/**************************************************************************\
**
**  tPtrList.h: Header file for tPtrList, tPtrListNode, and tPtrListIter
**              objects.
**
**  A tPtrList is an object that implements a general linked list of 
**  pointers to NodeType objects, where NodeType can be any type (double, 
**  int, other objects, etc). Lists of pointers are handled separately
**  from lists of non-pointer data types (which are handled by tList)
**  because of the requirements of pointers (specifically, in a normal
**  list you want to retrieve the actual item in the list, whereas in a 
**  pointer list you want to retrieve the item to which the list entry
**  points).
**
**  Pointer lists can be either linear or circular. The tPtrList class 
**  provides a variety of methods for adding, moving, and retrieving list 
**  elements. For moving back and forth in a tList and retrieving items, 
**  it's often most useful to use a tPtrListIter object (q.v.).
**
**  tPtrListNode objects are the nodes on the list; each contains a pointer
**  to the given data type (double, int, class, etc) and a pointer to the
**  next node in the list.
**
**  A tPtrListIter is an iterator for the linked list items (and their
**  descendants). Its services include fetching data from the current entry
**  on the list, advancing to the next or previous item on the list, etc.
**
**  See also tList, tArray, tMatrix
**
**  $Id: tPtrList.h,v 1.9 2000-01-13 23:56:06 gtucker Exp $
\**************************************************************************/

#ifndef TPTRLIST_H
#define TPTRLIST_H

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
//#include "../Classes.h" // TODO: include only needed stuff

template < class NodeType > class tPtrList;
template < class NodeType > class tPtrListIter;
template < class NodeType > class tMeshList;
template < class NodeType > class tMeshListIter;


/**************************************************************************\
** class tPtrListNode *****************************************************
**
** Class tPtrListNode represents the items (or "nodes") on the list. Each
** tPtrListNode object has two parts: the data pointer (of type NodeType *)
** and a pointer to the next item on the list. Capabilities include copy
** construction (from either another tPtrListNode or a NodeType *),
** returning a pointer or reference to either the data or the tPtrListNode
** itself, and assignment and equality/inequality operations.
**
\**************************************************************************/
template< class NodeType >
class tPtrListNode
{
    friend class tPtrList< NodeType >;
    friend class tPtrListIter< NodeType >;
public:
    tPtrListNode();                                    // default constructor
    tPtrListNode( const tPtrListNode< NodeType > & );  // copy constr #1
    tPtrListNode( NodeType * );                        // copy constr #2
    ~tPtrListNode();                                   // destructor
    const tPtrListNode< NodeType >
    &operator=( const tPtrListNode< NodeType > & );          // assignment
    int operator==( const tPtrListNode< NodeType > & ) const;  // equality
    int operator!=( const tPtrListNode< NodeType > & ) const;  // inequality
    NodeType *getPtrNC();                              // return data ptr
    const NodeType *getPtr() const;                    // return const data ptr
    tPtrListNode< NodeType > *getNextNC();             // return next item
    const tPtrListNode< NodeType > * getNext() const;  // return next as const
private:
    NodeType * Ptr;                   // ptr to data
    tPtrListNode< NodeType > * next;  // ptr to next list item
};



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
inline tPtrListNode< NodeType >::tPtrListNode()
{
   Ptr = 0;
   next = 0;
}

//copy constructor with tPtrListNode
template< class NodeType >                  //tPtrListNode
inline tPtrListNode< NodeType >::
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
inline tPtrListNode< NodeType >::
tPtrListNode( NodeType * NTPtr )
{
   Ptr = NTPtr;
   next = 0;
}

//destructor
template< class NodeType >                  //tPtrListNode
inline tPtrListNode< NodeType >::
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
inline const tPtrListNode< NodeType > &tPtrListNode< NodeType >::
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
inline int tPtrListNode< NodeType >::
operator==( const tPtrListNode< NodeType > &right ) const
{
   if( next != right.next ) return 0;
   if( Ptr != right.Ptr ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                  //tPtrListNode
inline int tPtrListNode< NodeType >::
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
inline const NodeType * tPtrListNode< NodeType >::
getPtr() const {return Ptr;}

template< class NodeType >                  //tPtrListNode
inline NodeType * tPtrListNode< NodeType >::
getPtrNC() {return Ptr;}

template< class NodeType >                  //tPtrListNode
inline const tPtrListNode< NodeType > *
tPtrListNode< NodeType >::
getNext() const {return next;}

template< class NodeType >                  //tPtrListNode
inline tPtrListNode< NodeType > * tPtrListNode< NodeType >::
getNextNC( ) {return next;}



/**************************************************************************\
** class tPtrList *********************************************************
**
** Class tPtrList implements a linked list of pointers. The class includes
** pointers to the first and last list nodes (see tPtrListNode) and the 
** number of items on the list.
**
\**************************************************************************/
/** class tPtrList **********************************************************/
template< class NodeType >
class tPtrList
{
    friend class tPtrListIter< NodeType >;
public:
    tPtrList();                               // default constructor
    tPtrList( const tPtrList< NodeType > & ); // copy constructor
    tPtrList( const tPtrList< NodeType > * ); // copy constructor
    ~tPtrList();                              // destructor
    const tPtrList< NodeType > 
    &operator=( const tPtrList< NodeType > & );  // assignment
    void insertAtFront( NodeType * ); // puts ptr at list front
    void insertAtBack( NodeType * );  // puts ptr at list back
    void insertAtNext( NodeType *, tPtrListNode< NodeType > * );
    void insertAtPrev( NodeType *, tPtrListNode< NodeType > * );
    int removeFromFront( NodeType * );  // removes 1st item, puts in ptr
    NodeType * removeFromFront();       // removes & returns 1st item
    int removeFromBack( NodeType * );   // removes last item, puts in ptr
    int removeNext( NodeType *, tPtrListNode< NodeType > * );
    int removePrev( NodeType *, tPtrListNode< NodeType > * );
    void Flush();         // clears & reinitializes list
    int isEmpty() const;  // returns 1 if list empty, 0 otherwise
    void print() const;   // prints list contents -- DEBUG ONLY
    /*void input( int, tList< NodeType > * );*/
    int getSize() const;  // returns size of list
    tPtrListNode< NodeType > * getFirstNC();  // returns ptr to 1st list node
    const tPtrListNode< NodeType > * getFirst() const;  // return const "
    void moveToBack( tPtrListNode< NodeType > *  );   // moves item to back
    void moveToFront( tPtrListNode< NodeType > *  );  // moves item to front
    tPtrListNode< NodeType > * getLast() const;
    void makeCircular();
    const NodeType *getIthPtr( int ) const;
    NodeType *getIthPtrNC( int ) const;
    tPtrList<NodeType> *DataCopy();   // copies AND CONTENTS POINTED TO (gt)
    
private:
    int nNodes;
    tPtrListNode< NodeType > * first;
    tPtrListNode< NodeType > * last;
    tPtrListNode< NodeType > * getNewNode( NodeType * );
};



/**************************************************************************\
**
**  tPtrList constructors & destructor:
**
**  Default constructor: initializes all values to 0 (empty list)
**  Copy constructors: each makes a complete copy of another tPtrList.
**                     (one takes a reference, the other a ptr)
**  Destructor: deletes all nodes on list. NOTE: does not destroy the
**              data items themselves!
**
**  Modifications:
**   - 2nd copy constructor added 1/2000, GT
\**************************************************************************/

//default constructor
template< class NodeType >                      //tPtrList
inline tPtrList< NodeType >::
tPtrList()
{
   first = last = 0;
   nNodes = 0;
     //cout << "tPtrList() instantiated" << endl;
}

//copy constructor
template< class NodeType >                      //tPtrList
inline tPtrList< NodeType >::
tPtrList( const tPtrList< NodeType > & orig )
{
   tPtrListNode< NodeType > *curNode = orig.first;
   first = last = 0;

   if( curNode != 0 )
   {
      insertAtBack( curNode->Ptr );
      for( curNode=curNode->next; curNode!=orig.first; curNode=curNode->next )
         insertAtBack( NTPtr );
      if( orig.last->next == orig.first ) last->next = first;
   }
}

//second copy constructor
template< class NodeType >
inline tPtrList< NodeType >::
tPtrList( const tPtrList< NodeType > * origptr )
{
   tPtrListNode< NodeType > *curNode = origptr->first;
   first = last = 0;

   if( curNode != 0 )
   {
      insertAtBack( curNode->Ptr );
      for( curNode=curNode->next; curNode!=orig.first; curNode=curNode->next )
         insertAtBack( NTPtr );
      if( origptr->last->next == origptr->first ) last->next = first;
   }
}
   

//destructor
template< class NodeType >                      //tPtrList
inline tPtrList< NodeType >::
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
inline const tPtrList< NodeType > &tPtrList< NodeType >::
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
inline tPtrListNode< NodeType > * tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
insertAtBack( NodeType *NTPtr )
{
   tPtrListNode< NodeType > * newPtr = getNewNode( NTPtr );
     //cout << "add new node to tPtrList" << endl;
   assert( this != 0 );
   if( isEmpty() )
      first = last = newPtr;
   else
   {
      newPtr->next = last->next;
      last->next = newPtr;
      last = newPtr;
   }

}

template< class NodeType >                      //tPtrList
inline void tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
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
**  Version 1: Removes the first item on the list and points NTPtr to the
**  new 1st item. Returns 0 if the list is already empty, 1 otherwise. Note
**  that if the list is empty, NTPtr is unchanged.
**
**  Version 2: Removes the first item on the list and returns it.
**
**  ALERT: There is a potential bug here: if the list is circular but
**  contains only one item (which points to itself), NTPtr will contain
**  a dangling pointer! TODO (gt)
**
**  Modifications:
**   - version 2 added 1/2000, GT
**
\**************************************************************************/
template< class NodeType >                      //tPtrList
inline int tPtrList< NodeType >::
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
inline NodeType* tPtrList< NodeType >::
removeFromFront( )
{
   NodeType *NTPtr;
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
      return NTPtr;
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
inline int tPtrList< NodeType >::
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
inline int tPtrList< NodeType >::
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
inline int tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
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
inline int tPtrList< NodeType >::
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
inline int tPtrList< NodeType >::
getSize() const
{
   /*int i;
   tPtrListNode< NodeType > *temp;
   for( i = 1, temp = first;
        temp != last; i++, temp = temp->next );*/
   return nNodes; //i;
}

template< class NodeType >                      //tPtrList
inline tPtrListNode< NodeType > * tPtrList< NodeType >::
getFirstNC() {return first;}

template< class NodeType >                      //tPtrList
inline const tPtrListNode< NodeType > * tPtrList< NodeType >::
getFirst() const {return first;}

template< class NodeType >                      //tPtrList
inline tPtrListNode< NodeType > * tPtrList< NodeType >::
getLast() const {return last;}

template< class NodeType >                      //tPtrList
inline const NodeType *tPtrList< NodeType >::
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
inline NodeType *tPtrList< NodeType >::
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
inline void tPtrList< NodeType >::
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
inline tPtrList< NodeType > *tPtrList<NodeType>::
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

/*X
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
*/

/**************************************************************************\
** class tPtrListIter ******************************************************
**
** Helper class for tPtrList. tPtrListIters are "iterators" that walk up &
** down a tPtrList, fetching items. Their chief advantage is that you can
** have multiple iterators on any given list at once, and thus multiple
** access points. Use of iterator classes is discussed by Deitel and
** Deitel, _C++ How to Program_, first edition, Prentice Hall, 1994.
**
** Note that in the current implementation, list items are fetched by
** ID number, which presupposes that the list items have a member function
** getID. This restricts the generality of tPtrList, and should be moved
** to tGridList. (TODO)
**
** Modifications:
**  - added 2nd copy constructor, GT, 1/2000
**
\**************************************************************************/
//TO DO: make Get, Where, GetP, refer to place in list rather than use getID()
template< class NodeType >
class tPtrListIter
{
  public:
   tPtrListIter();
   tPtrListIter( tPtrList< NodeType > & );
   tPtrListIter( tPtrList< NodeType > * );
   ~tPtrListIter();
   int First();
   int Last();
   int Get( int );
   int Get( NodeType * );
   int Where();
   NodeType *DatPtr();
   tPtrListNode< NodeType > *NodePtr();
   int Next();
   int Prev();
   int NextIsNotFirst();
   void Reset( tPtrList< NodeType > & );
   NodeType *NextP();
   NodeType *GetP( int );
   NodeType *FirstP();
   NodeType *LastP();
   NodeType *PrevP();
   NodeType *ReportNextP();
   NodeType *ReportPrevP();
   int AtEnd();
  private:
   tPtrList< NodeType > * ptrlistPtr;
   tPtrListNode< NodeType > * curptrnode;
   int counter;
};


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
inline tPtrListIter< NodeType >::
tPtrListIter()
{
   ptrlistPtr = 0;
   curptrnode = 0;
   counter = 0;
     //cout << "tPtrListIter()" << endl;
}

template< class NodeType >    //tPtrListIter
inline tPtrListIter< NodeType >::
tPtrListIter( tPtrList< NodeType > &ptrlist )
{
   assert( &ptrlist != 0 );
   ptrlistPtr = &ptrlist;
   curptrnode = ptrlist.first;
   counter = 0;
     //cout << "tPtrListIter( ptrlist )" << endl;
}

template< class NodeType >    //tPtrListIter
inline tPtrListIter< NodeType >::
tPtrListIter( tPtrList< NodeType > * ptrlist )
{
   assert( ptrlist != 0 );
   ptrlistPtr = ptrlist;
   curptrnode = ptrlist->first;
   counter = 0;
}

template< class NodeType >    //tPtrListIter
inline tPtrListIter< NodeType >::
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
inline void tPtrListIter< NodeType >::
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
inline int tPtrListIter< NodeType >::
First()
{
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->first;
   counter = 0;
   if( curptrnode != 0 ) return 1;
   else return 0;
}

template< class NodeType >     //tPtrListIter
inline int tPtrListIter< NodeType >::
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
**  (1) Move to list item with ID number _num_. Note: assumes that list items
**  have a member function getID()! Returns 1 if found, 0 if not.
**
**  (2) Move to list item containing pointer to desiredItemPtr; return
**  zero if not found. (added by GT, 1/2000)
**
\**************************************************************************/
template< class NodeType >     //tPtrListIter
inline int tPtrListIter< NodeType >::
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

template< class NodeType >
inline int tPtrListIter< NodeType >::
Get( NodeType *desiredItemPtr )
{
   tPtrListNode<NodeType> *cln = 0;
   int listSize = ptrListPtr->getSize();

   // Go down the list looking for a node that points to desiredItemPtr
   // If found, move to it and return 1; otherwise return 0
   for( cln=ptrListPtr->first, ctr=0; ctr<listSize; cln=cln->next, ctr++ )
       if( cln->Ptr == desiredItemPtr )
       {
          curptrnode = cln;
          return 1;
       }
   return 0;
   
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
inline int tPtrListIter< NodeType >::
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
inline int tPtrListIter< NodeType >::
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
inline int tPtrListIter< NodeType >::
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
inline NodeType *tPtrListIter< NodeType >::
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
inline tPtrListNode< NodeType > *tPtrListIter< NodeType >::
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
inline int tPtrListIter< NodeType >::
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
inline NodeType * tPtrListIter< NodeType >::
FirstP()
{
   assert( ptrlistPtr != 0 );
   curptrnode = ptrlistPtr->first;
   counter = 0;
   if( curptrnode != 0 ) return curptrnode->Ptr;
   else return 0;
}
   
template< class NodeType >        //tListIter
inline NodeType * tPtrListIter< NodeType >::
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
inline NodeType *tPtrListIter< NodeType >::
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
inline NodeType * tPtrListIter< NodeType >::
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
inline NodeType * tPtrListIter< NodeType >::
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
inline NodeType * tPtrListIter< NodeType >::
ReportNextP()
{
   assert( ptrlistPtr != 0 );
   if( curptrnode == 0 ) return 0;
   if( curptrnode->next != 0 ) return curptrnode->next->Ptr;
   else return 0;
}

template< class NodeType >        //tListIter
inline NodeType * tPtrListIter< NodeType >::
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
inline int tPtrListIter< NodeType >::
AtEnd()
{
   if( ptrlistPtr->last == 0 ) return 1;
   if( ptrlistPtr->last->next == 0 ) return ( curptrnode==0 );
   else return ( curptrnode == ptrlistPtr->first && counter != 0 );
}

#endif


