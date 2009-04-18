//-*-c++-*-

/**************************************************************************/
/**
 **  @file tPtrList.h
 **  @brief Header file for tPtrList, tPtrListNode, and tPtrListIter
 **         objects.
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
 **  Modifications:
 **    - 3/31/00 bug fix to tPtrList copy constructors (GT)
 **    - 5/10/00 typo fix in DataCopy (GT)
 **    - 1/99: (SL)added member "prev" to tPtrListNode to make tPtrList a
 **             doubly linked list; other than addition of
 **             tPtrListNode::getPrev(), getPrevNC(), interface is unchanged
 **      9/02: (AD)merge in main Child version
 **
 **  $Id: tPtrList.h,v 1.49 2004/06/16 13:37:41 childcvs Exp $
 */
/**************************************************************************/

#ifndef TPTRLIST_H
#define TPTRLIST_H

#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include "../compiler.h"

// do not support these ill-defined functions
#undef SUPPORT_DEPRECATED

template < class NodeType > class tPtrList;
template < class NodeType > class tPtrListIter;

/*
  class PtrListDebug
  {
  public:
  static void TellAll();
  static int dbgPRV;
  static int dbgPRP;
  static int dbgIAP;
  static int dbgRP;
  static int dbgMTB;
  static int dbgMTF;
  };


  int PtrListDebug::dbgIAP = 0;
  int PtrListDebug::dbgRP = 0;
  int PtrListDebug::dbgMTB = 0;
  int PtrListDebug::dbgMTF = 0;
  int PtrListDebug::dbgPRV = 0;
  int PtrListDebug::dbgPRP = 0;

  void PtrListDebug::TellAll()
  {
  //debug
  std::cout << "IAP " << dbgIAP << std::endl;
  std::cout << "RP " << dbgRP << std::endl;
  std::cout << "MTB " << dbgMTB << std::endl;
  std::cout << "MTF " << dbgMTF << std::endl;
  std::cout << "PRV " << dbgPRV << std::endl;
  std::cout << "PRP " << dbgPRP << std::endl;

  }*/




/**************************************************************************/
/**
 ** @class tPtrListNode
 **
 ** Class tPtrListNode represents the items (or "nodes") on the list. Each
 ** tPtrListNode object has two parts: the data pointer (of type NodeType *)
 ** and a pointer to the next item on the list. Capabilities include copy
 ** construction (from either another tPtrListNode or a NodeType *),
 ** returning a pointer or reference to either the data or the tPtrListNode
 ** itself, and assignment and equality/inequality operations.
 **
 */
/**************************************************************************/
template< class NodeType >
class tPtrListNode
{
  friend class tPtrList< NodeType >;
  friend class tPtrListIter< NodeType >;
public:
  inline tPtrListNode();                                    // default constructor
  inline tPtrListNode( const tPtrListNode< NodeType > & );  // copy constr #1
  inline tPtrListNode( NodeType * );                        // copy constr #2
  inline ~tPtrListNode();                            // destructor
  const tPtrListNode< NodeType >
  &operator=( const tPtrListNode< NodeType > & );          // assignment
  inline bool operator==( const tPtrListNode< NodeType > & ) const;  // equality
  inline bool operator!=( const tPtrListNode< NodeType > & ) const;  // inequality
  inline NodeType *getPtrNC();                       // return data ptr
  inline const NodeType *getPtr() const;             // return const data ptr
  inline tPtrListNode< NodeType > *getNextNC();             // return next item
  inline const tPtrListNode< NodeType > * getNext() const;  // return next as const
  inline tPtrListNode< NodeType > *getPrevNC();
  inline const tPtrListNode< NodeType > * getPrev() const;
private:
  NodeType * Ptr;                   // ptr to data
  tPtrListNode< NodeType > * next;  // ptr to next list item
  tPtrListNode< NodeType > * prev;
};


/**************************************************************************/
/**
 ** @class tPtrList
 **
 ** Class tPtrList implements a linked list of pointers. The class includes
 ** pointers to the first and last list nodes (see tPtrListNode) and the
 ** number of items on the list.
 **
 */
/**************************************************************************/
template< class NodeType >
class tPtrList
{
  friend class tPtrListIter< NodeType >;
public:
  inline tPtrList();                        // default constructor
  tPtrList( const tPtrList< NodeType > & ); // copy constructor
  tPtrList( const tPtrList< NodeType > * ); // copy constructor
  ~tPtrList();                              // destructor
  const tPtrList< NodeType >
  &operator=( const tPtrList< NodeType > & );  // assignment
  inline void insertAtFront( NodeType * ); // puts ptr at list front
  inline void insertAtBack( NodeType * );  // puts ptr at list back
  inline void insertAtNext( NodeType *, tPtrListNode< NodeType > * );
  inline void insertAtPrev( NodeType *, tPtrListNode< NodeType > * );
#if defined(SUPPORT_DEPRECATED)
  int removeFromFront( NodeType * );  // removes 1st item, puts in ptr
#endif
  inline NodeType * removeFromFront();       // removes & returns 1st item
#if defined(SUPPORT_DEPRECATED)
  int removeFromBack( NodeType * );   // removes last item, puts in ptr
#endif
  inline NodeType * removeFromBack();        // removes & returns last item
#if defined(SUPPORT_DEPRECATED)
  int removeNext( NodeType *, tPtrListNode< NodeType > * );
#endif
  inline NodeType* removeNext(tPtrListNode< NodeType > * );
#if defined(SUPPORT_DEPRECATED)
  int removePrev( NodeType *, tPtrListNode< NodeType > * );
#endif
  inline NodeType* removePrev( tPtrListNode< NodeType > * );
  void Flush();         // clears & reinitializes list
  inline bool isEmpty() const;  // returns true if list empty, false otherwise
  void print() const;   // prints list contents -- DEBUG ONLY
  inline int getSize() const;  // returns size of list
  inline tPtrListNode< NodeType > * getFirstNC();  // returns ptr to 1st list node
  inline const tPtrListNode< NodeType > * getFirst() const;  // return const "
  inline void moveToBack( tPtrListNode< NodeType > *  );   // moves item to back
  inline void moveToFront( tPtrListNode< NodeType > *  );  // moves item to front
  inline tPtrListNode< NodeType > const * getLast() const;
  inline void makeCircular();
  inline const NodeType *getIthPtr( int ) const;
  inline NodeType *getIthPtrNC( int );
  const tPtrListNode< NodeType >* getIth( int ) const;
  tPtrListNode< NodeType >* getIthNC( int );

private:
  int nNodes;
  tPtrListNode< NodeType > * first;
  tPtrListNode< NodeType > * last;
  tPtrListNode< NodeType > * getNewNode( NodeType * );
};


/**************************************************************************/
/**
 ** @class tPtrListIter
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
 ** to tMeshList. (TODO)
 **
 ** Modifications:
 **  - added 2nd copy constructor, GT, 1/2000
 **
 */
/**************************************************************************/
//TO DO: make Get, Where, GetP, refer to place in list rather than use getID()
template< class NodeType >
class tPtrListIter
{
  tPtrListIter& operator=( const tPtrListIter& );
  int NextIfNoCurrent();  // set 1st as current undefined
  int PrevIfNoCurrent();  // set last as current undefined
public:
  inline tPtrListIter();
  inline tPtrListIter( const tPtrListIter< NodeType > & );
  inline tPtrListIter( tPtrList< NodeType > & );
  inline tPtrListIter( tPtrList< NodeType > * );
  inline ~tPtrListIter();
  inline int First();
  inline int Last();
  int Get( int );
  int Get( const NodeType * );
  inline int Where() const;
  inline NodeType *DatPtr() const;
  inline tPtrListNode< NodeType > *NodePtr();
  inline int Next();
  inline int Prev();
  inline int NextIsNotFirst() const;
  inline void Reset( tPtrList< NodeType > & );
  inline NodeType *NextP();
  inline NodeType *GetP( int ); //use only if NodeType has member getID()!!
  inline NodeType* GetP( NodeType* );
  inline NodeType *FirstP();
  inline NodeType *LastP();
  inline NodeType *PrevP();
  inline NodeType *ReportNextP() const;
  inline NodeType *ReportPrevP() const;
  inline bool AtEnd() const;
private:
  tPtrList< NodeType > * ptrlistPtr;
  tPtrListNode< NodeType > * curptrnode;
  int counter;
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
template< class NodeType >
inline tPtrListNode< NodeType >::tPtrListNode() :
  Ptr(0),
  next(0),
  prev(0)
{}

//copy constructor with tPtrListNode
template< class NodeType >
inline tPtrListNode< NodeType >::
tPtrListNode( const tPtrListNode< NodeType > & init ) :
  Ptr(init.Ptr),
  next(init.next),
  prev(init.prev)
{}

// copy constructor with data ptr
template< class NodeType >
inline tPtrListNode< NodeType >::
tPtrListNode( NodeType * NTPtr ) :
  Ptr(NTPtr),
  next(0),
  prev(0)
{}

//destructor
template< class NodeType >
inline tPtrListNode< NodeType >::
~tPtrListNode()
{
  Ptr = 0;  //redirect the pointer away from its target
  next = prev = 0;
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
template< class NodeType >
inline const tPtrListNode< NodeType > &tPtrListNode< NodeType >::
operator=( const tPtrListNode< NodeType > &right )
{
  if( &right != this )
    {
      Ptr = right.Ptr;
      next = right.next;
      prev = right.prev;
    }
  return *this;
}

//overloaded equality operator:
template< class NodeType >
inline bool tPtrListNode< NodeType >::
operator==( const tPtrListNode< NodeType > &right ) const
{
  if( next != right.next ) return false;
  if( prev != right.prev ) return false;
  if( Ptr != right.Ptr ) return false;
  return true;
}

//overloaded inequality operator:
template< class NodeType >
inline bool tPtrListNode< NodeType >::
operator!=( const tPtrListNode< NodeType > &right ) const
{
  return ! operator==(right);
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

template< class NodeType >
inline const NodeType * tPtrListNode< NodeType >::
getPtr() const {return Ptr;}

template< class NodeType >
inline NodeType * tPtrListNode< NodeType >::
getPtrNC() {return Ptr;}

template< class NodeType >
inline const tPtrListNode< NodeType > *
tPtrListNode< NodeType >::
getNext() const {return next;}

template< class NodeType >
inline tPtrListNode< NodeType > * tPtrListNode< NodeType >::
getNextNC( ) {return next;}

template< class NodeType >
inline const tPtrListNode< NodeType >*tPtrListNode< NodeType >::
getPrev() const {return prev;}

template< class NodeType >
inline tPtrListNode< NodeType >*tPtrListNode< NodeType >::
getPrevNC( ) {return prev;}


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
 **   - bug in both copy constructors: use of undeclared variable NTPtr.
 **     Fixed 3/00, GT
\**************************************************************************/

//default constructor
template< class NodeType >
inline tPtrList< NodeType >::
tPtrList() :
  nNodes(0), first(0), last(0)
{
  if (0) //DEBUG
    std::cout << "tPtrList() instantiated" << std::endl;
}

//copy constructor
template< class NodeType >
tPtrList< NodeType >::
tPtrList( const tPtrList< NodeType > & orig ) :
  nNodes(0), first(0), last(0)
{
  tPtrListNode< NodeType > *curNode = orig.first;

  if( curNode != 0 )
    {
      insertAtBack( curNode->Ptr );
      for( curNode=curNode->next; curNode!=orig.last->next; curNode=curNode->next )
	insertAtBack( curNode->Ptr );
      if( orig.last->next == orig.first ) last->next = first;
      if( orig.first->prev == orig.last ) first->prev = last;
    }
}

//second copy constructor
template< class NodeType >
tPtrList< NodeType >::
tPtrList( const tPtrList< NodeType > * origptr ) :
  nNodes(0), first(0), last(0)
{
  tPtrListNode< NodeType > *curNode = origptr->first;

  if( curNode != 0 )
    {
      insertAtBack( curNode->Ptr );
      for( curNode=curNode->next; curNode!=origptr->last->next; curNode=curNode->next )
	insertAtBack( curNode->Ptr );
      if( origptr->last->next == origptr->first ) last->next = first;
      if( origptr->first->prev == origptr->last ) first->prev = last;
    }
}


//destructor
template< class NodeType >
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
	  if( current != last ) current = current->next;
	  else
	    {
	      current = 0;
	      last = 0;
	    }
	  delete temp;
	}
    }
  first = 0;
  last = 0;
}



/**************************************************************************\
 **
 **  tPtrList overloaded operators:
 **
 **  Assignment: clears the list and makes a copy of the right-hand list
 **
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >
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
	  if( right.first->prev == right.last ) first->prev = last;
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
template< class NodeType >
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

template< class NodeType >
inline void tPtrList< NodeType >::
insertAtFront( NodeType *NTPtr )
{
  tPtrListNode< NodeType > *newPtr = getNewNode( NTPtr );
  if( isEmpty() )
    first = last = newPtr;
  else
    {
      newPtr->next = first;
      newPtr->prev = first->prev;
      first->prev = newPtr;
      if( last->next == first )
	last->next = newPtr;
      first = newPtr;
    }
}

template< class NodeType >
inline void tPtrList< NodeType >::
insertAtBack( NodeType *NTPtr )
{
  tPtrListNode< NodeType > * newPtr = getNewNode( NTPtr );
  assert( this != 0 );
  if( isEmpty() )
    first = last = newPtr;
  else
    {
      newPtr->next = last->next;
      newPtr->prev = last;
      last->next = newPtr;
      if( first->prev == last )
	first->prev = newPtr;
      last = newPtr;
    }
}

template< class NodeType >
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
      newPtr->prev = prev;
      prev->next->prev = newPtr;
      prev->next = newPtr;
    }
}

template< class NodeType >
inline void tPtrList< NodeType >::
insertAtPrev( NodeType *NTPtr, tPtrListNode< NodeType > * node )
{
  if( node != 0 )
    {
      if( node == first )
	{
	  insertAtFront( NTPtr );
	  return;
	}
      tPtrListNode< NodeType > * newPtr = getNewNode( NTPtr );
      newPtr->next = node;
      newPtr->prev = node->prev;
      node->prev->next = newPtr;
      node->prev = newPtr;
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
#if defined(SUPPORT_DEPRECATED)
template< class NodeType >
inline int tPtrList< NodeType >::
removeFromFront( NodeType * /*NTPtr*/ )
{
  NodeType * temp = removeFromFront();
  //NTPtr = temp;
  return temp != 0 ? 1:0;
}
#endif

template< class NodeType >
inline NodeType* tPtrList< NodeType >::
removeFromFront()
{
  NodeType *NTPtr = 0;
  if( isEmpty() ) return 0;
  else
    {
      tPtrListNode< NodeType > * temp = first;
      if( first == last ) first = last = 0;
      else
	{
	  first->next->prev = first->prev;
	  if( last->next == first ) last->next = first->next;
	  first = first->next;
	}
      NTPtr = temp->Ptr;
      delete temp;
      --nNodes;
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
#if defined(SUPPORT_DEPRECATED)
template< class NodeType >
inline int tPtrList< NodeType >::
removeFromBack( NodeType * /*NTPtr*/ )
{
  NodeType * temp = removeFromBack();
  //NTPtr = temp;
  return temp != 0 ? 1:0;
}
#endif

template< class NodeType >
inline NodeType* tPtrList< NodeType >::
removeFromBack()
{
  NodeType *NTPtr = 0;
  if( isEmpty() ) return 0;
  else
    {
      tPtrListNode< NodeType > * temp = last;
      if( first == last ) first = last = 0;
      else
	{
	  last->prev->next = last->next;
	  if( first->prev == last ) first->prev = last->prev;
	  last = last->prev;
	}
      NTPtr = temp->Ptr;
      delete temp;
      --nNodes;
      return NTPtr;
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
#if defined(SUPPORT_DEPRECATED)
template< class NodeType >
inline int tPtrList< NodeType >::
removeNext( NodeType * /*NTPtr*/, tPtrListNode< NodeType > * ptr )
{
  NodeType * temp = removeNext( ptr );
  //NTPtr = temp;
  return temp != 0 ? 1:0;
}
#endif

template< class NodeType >
inline NodeType* tPtrList< NodeType >::
removeNext( tPtrListNode< NodeType >* ptr )
{
  NodeType* NTPtr = 0;
  if( ptr == 0 ) return NTPtr;
  if( ptr->next == 0 ) return NTPtr;
  if( ptr->next == last ) return removeFromBack();
  if( ptr->next == first ) return removeFromFront();
  tPtrListNode< NodeType > * temp = ptr->next;
  ptr->next = ptr->next->next;
  ptr->next->prev = ptr;
  NTPtr = temp->Ptr;
  delete temp;
  --nNodes;
  return NTPtr;
}

/**************************************************************************\
 **
 **  tPtrList::removePrev
 **
 **  Removes the item before _ptr_ on the list, returning the pointer in
 **  NTPtr.
 **
\**************************************************************************/
#if defined(SUPPORT_DEPRECATED)
template< class NodeType >
inline int tPtrList< NodeType >::
removePrev( NodeType * /*NTPtr*/, tPtrListNode< NodeType > * ptr )
{
  NodeType * temp = removePrev( ptr );
  //NTPtr = temp;
  return temp != 0 ? 1:0;
}
#endif

template< class NodeType >
inline NodeType* tPtrList< NodeType >::
removePrev( tPtrListNode< NodeType > * ptr )
{
  NodeType* NTPtr = 0;
  if( ptr == 0 ) return NTPtr;
  if( ptr->prev == 0 ) return NTPtr;
  if( ptr->prev == last ) return removeFromBack();
  if( ptr->prev == first ) return removeFromFront();
  tPtrListNode< NodeType >* temp = ptr->prev;
  ptr->prev = ptr->prev->prev;
  ptr->prev->next = ptr;
  NTPtr = temp->Ptr;
  delete temp;
  --nNodes;
  return NTPtr;
}

/**************************************************************************\
 **
 **  tPtrList::moveToFront and moveToBack
 **
 **  Moves the list item pointed to by mvnode to the front or back of
 **  the list, respectively.
 **
\**************************************************************************/

template< class NodeType >
inline void tPtrList< NodeType >::
moveToBack( tPtrListNode< NodeType > * mvnode )
{
  if( mvnode != last )
    {
      if( mvnode == first )
	{
	  first->next->prev = first->prev;
	  first = first->next;
	}
      else
	{
	  mvnode->prev->next = mvnode->next;
	  mvnode->next->prev = mvnode->prev;
	}
      mvnode->next = last->next;
      mvnode->prev = last;
      if( first->prev != 0 ) first->prev = mvnode;
      last->next = mvnode;
      last = mvnode;
    }
}

template< class NodeType >
inline void tPtrList< NodeType >::
moveToFront( tPtrListNode< NodeType > * mvnode )
{
  if( mvnode != first )
    {
      if( mvnode == last )
	{
	  last->prev->next = last->next;
	  last = last->prev;
	}
      else
	{
	  mvnode->prev->next = mvnode->next;
	  mvnode->next->prev = mvnode->prev;
	}
      mvnode->next = first;
      mvnode->prev = first->prev;
      if( last->next != 0 ) last->next = mvnode;
      first->prev = mvnode;
      first = mvnode;
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
template< class NodeType >
void tPtrList< NodeType >::
Flush()
{
  assert( this!=0 );
  if( !isEmpty() )
    {
      tPtrListNode<NodeType > *current = first, * temp;
      first = 0;
      while( last != 0 )
	{
	  temp = current;
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
template< class NodeType >
inline bool tPtrList< NodeType >::
isEmpty() const
{
  return BOOL( first == 0 );
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
template< class NodeType >
inline int tPtrList< NodeType >::
getSize() const { return nNodes;}

template< class NodeType >
inline tPtrListNode< NodeType > * tPtrList< NodeType >::
getFirstNC() {return first;}

template< class NodeType >
inline const tPtrListNode< NodeType > * tPtrList< NodeType >::
getFirst() const {return first;}

template< class NodeType >
inline const tPtrListNode< NodeType > * tPtrList< NodeType >::
getLast() const {return last;}

template< class NodeType >
inline const NodeType *tPtrList< NodeType >::
getIthPtr( int num ) const
{
  return getIth( num )->getPtr();
}

template< class NodeType >
inline NodeType *tPtrList< NodeType >::
getIthPtrNC( int num )
{
  return getIthNC( num )->getPtrNC();
}

template< class NodeType >
const tPtrListNode< NodeType >* tPtrList< NodeType >::
getIth( int num ) const
{
  int i;
  tPtrListNode< NodeType > const * curPtr;
  assert( num >= 0 && num < nNodes );
  if( num > nNodes / 2 )
    for( curPtr = last, i = nNodes-1; i>num; curPtr = curPtr->prev, --i );
  else
    for( curPtr = first, i = 0; i<num; curPtr = curPtr->next, ++i );
  return curPtr;
}

template< class NodeType >
tPtrListNode< NodeType >* tPtrList< NodeType >::
getIthNC( int num )
{
  int i;
  tPtrListNode< NodeType >* curPtr;
  assert( num >= 0 && num < nNodes );
  if( num > nNodes / 2 )
    for( curPtr = last, i = nNodes-1; i>num; curPtr = curPtr->prev, --i );
  else
    for( curPtr = first, i = 0; i<num; curPtr = curPtr->next, ++i );
  return curPtr;
}

/**************************************************************************\
 **
 **  tPtrList::makeCircular
 **
 **  Converts the list into a circular list by having the last item point
 **  to the first.
 **
\**************************************************************************/
template< class NodeType >
inline void tPtrList< NodeType >::
makeCircular()
{
  assert( first != 0 && last != 0 );
  last->next = first;
  first->prev = last;
}


//display list contents
template< class NodeType >
void tPtrList< NodeType >::
print() const
{
  if( isEmpty() )
    {
      std::cout<<"The list is empty\n"<<std::endl;
      return;
    }
  tPtrListNode< NodeType > * current = first;
  std::cout<<"The list is: ";
  while( current != 0 )
    {
      std::cout<<current->Ptr->getID() <<' ';
      current = current->next;
    }
  std::cout << '\n' <<std::endl;
}


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

template< class NodeType >
inline tPtrListIter< NodeType >::
tPtrListIter() :
  ptrlistPtr(0),
  curptrnode(0),
  counter(0)
{
  if (0) //DEBUG
    std::cout << "tPtrListIter()" << std::endl;
}

template< class NodeType >
inline tPtrListIter< NodeType >::
tPtrListIter(const tPtrListIter< NodeType >& orig) :
  ptrlistPtr(orig.ptrlistPtr),
  curptrnode(orig.curptrnode),
  counter(orig.counter)
{
  if (0) //DEBUG
    std::cout << "tPtrListIter(const tPtrListIter&)" << std::endl;
}

template< class NodeType >
inline tPtrListIter< NodeType >::
tPtrListIter( tPtrList< NodeType > &ptrlist ) :
  ptrlistPtr(&ptrlist),
  curptrnode(ptrlist.first),
  counter(0)
{
  if (0) //DEBUG
    std::cout << "tPtrListIter( ptrlist )" << std::endl;
}

template< class NodeType >
inline tPtrListIter< NodeType >::
tPtrListIter( tPtrList< NodeType > * ptrlist ) :
  ptrlistPtr(ptrlist),
  curptrnode(0),
  counter(0)
{
  assert( ptrlist != 0 );
  curptrnode = ptrlist->first;
}

template< class NodeType >
inline tPtrListIter< NodeType >::
~tPtrListIter()
{
  ptrlistPtr = 0;
  curptrnode = 0;
  if (0) //DEBUG
    std::cout << "    ~tPtrListIter()" << std::endl;
}

/**************************************************************************\
 **
 **  tPtrListIter::Reset
 **
 **  Points iterator at the 1st node on _ptrlist_ (provides a way of telling
 **  an iterator which list to work on).
 **
\**************************************************************************/
template< class NodeType >
inline void tPtrListIter< NodeType >::
Reset( tPtrList< NodeType > &ptrlist )
{
  assert( &ptrlist != 0 );
  ptrlistPtr = &ptrlist;
  assert( ptrlistPtr != 0 );
  curptrnode = ptrlistPtr->first;
  counter = 0;
}

/**************************************************************************\
 **
 **  tPtrListIter::First and Last
 **
 **  Move to the first or last item on the current list. Return TRUE if
 **  pointing to a valid list item, FALSE otherwise.
 **
\**************************************************************************/
template< class NodeType >
inline int tPtrListIter< NodeType >::
First()
{
  assert( ptrlistPtr != 0 );
  curptrnode = ptrlistPtr->first;
  counter = 0;
  if( curptrnode != 0 ) return 1;
  return 0;
}

template< class NodeType >
inline int tPtrListIter< NodeType >::
Last()
{
  assert( ptrlistPtr != 0 );
  curptrnode = ptrlistPtr->last;
  counter = ptrlistPtr->nNodes - 1;
  if( curptrnode != 0 ) return 1;
  return 0;
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
template< class NodeType >
int tPtrListIter< NodeType >::
Get( int num )
{
  assert( ptrlistPtr != 0 );
  tPtrListNode< NodeType > *tempnodeptr;
  for( tempnodeptr = ptrlistPtr->first, counter = 0;
       counter <= ptrlistPtr->getSize() && tempnodeptr != 0;
       tempnodeptr = tempnodeptr->next, ++counter )
    {
      if( tempnodeptr->getPtr()->getID() == num )
	break;
    }
  if( tempnodeptr == 0 ) return 0;
  if( tempnodeptr->getPtr()->getID() != num ) return 0;
  curptrnode = tempnodeptr;
  return 1;
}

template< class NodeType >
int tPtrListIter< NodeType >::
Get( const NodeType *desiredItemPtr )
{
  assert( ptrlistPtr != 0 );
  tPtrListNode<NodeType> *tempnodeptr;
  for( tempnodeptr=ptrlistPtr->first, counter = 0;
       counter <= ptrlistPtr->getSize() && tempnodeptr != 0;
       tempnodeptr=tempnodeptr->next, ++counter )
    {
      if( tempnodeptr->Ptr == desiredItemPtr )
	break;
    }
  if( tempnodeptr == 0 ) return 0;
  if( tempnodeptr->Ptr != desiredItemPtr ) return 0;
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
template< class NodeType >
inline int tPtrListIter< NodeType >::
Next()
{
  if( likely(curptrnode != 0) )
    {
      curptrnode = curptrnode->next;
      ++counter;
      return curptrnode != 0 ? 1:0;
    }
  // unlikely case
  return NextIfNoCurrent();
}

template< class NodeType >
int tPtrListIter< NodeType >::
NextIfNoCurrent()
{
  assert( ptrlistPtr != 0 );
  assert( curptrnode == 0 );
  curptrnode = ptrlistPtr->first;
  counter = 0;
  return curptrnode !=0 ? 1:0;
}

template< class NodeType >
inline int tPtrListIter< NodeType >::
Prev()
{
  if( likely(curptrnode != 0) )
    {
      curptrnode = curptrnode->prev;
      --counter;
      return curptrnode != 0 ? 1:0;
    }
  // unlikely case
  return PrevIfNoCurrent();
}

template< class NodeType >
int tPtrListIter< NodeType >::
PrevIfNoCurrent()
{
  assert( ptrlistPtr != 0 );
  assert( curptrnode == 0 );
  curptrnode = ptrlistPtr->last;
  counter = -1;
  return curptrnode != 0 ? 1:0;
}

/**************************************************************************\
 **
 **  tPtrListIter::Where
 **
 **  Returns the ID number of the current data item, or -1 if there is
 **  no current data item. Assumes data item has a getID() mbr function!
 **
\**************************************************************************/
template< class NodeType >
inline int tPtrListIter< NodeType >::
Where() const
{
  return ( curptrnode == 0 ) ?
    -1 : curptrnode->getPtr()->getID();
}


/**************************************************************************\
 **
 **  tPtrListIter::DatPtr
 **
 **  Returns copy of current data pointer.
 **
\**************************************************************************/
template< class NodeType >
inline NodeType *tPtrListIter< NodeType >::
DatPtr() const
{
  return ( curptrnode == 0 ) ?
    0 : curptrnode->Ptr;
}


/**************************************************************************\
 **
 **  tPtrListIter::NodePtr
 **
 **  Returns pointer to current list node.
 **
\**************************************************************************/
template< class NodeType >
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
template< class NodeType >
inline int tPtrListIter< NodeType >::
NextIsNotFirst() const
{
  assert( curptrnode != 0 );
  assert( ptrlistPtr != 0 );
  return ( curptrnode->next == ptrlistPtr->first ) ? 0:1;
}


/**************************************************************************\
 **
 **  tPtrListIter::FirstP and tPtrListIter::LastP
 **
 **  Move to the first or last item on the list and return a copy of the
 **  data pointer, or 0 if first/last item is empty.
 **
\**************************************************************************/
template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::
FirstP()
{
  return ( First() ) ?
    curptrnode->Ptr : 0;
}

template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::
LastP()
{
  return ( Last() ) ?
    curptrnode->Ptr : 0;
}


/**************************************************************************\
 **
 **  tPtrListIter::NextP and tPtrListIter::PrevP
 **
 **  Same as Next and Prev, except that the functions return a copy of the
 **  data pointer rather than a pointer to the list item.
 **
\**************************************************************************/
template< class NodeType >
inline NodeType *tPtrListIter< NodeType >::
PrevP()
{
  return ( Prev() ) ?
    curptrnode->Ptr : 0;
}

template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::
NextP()
{
  return ( Next() ) ?
    curptrnode->Ptr : 0;
}


/**************************************************************************\
 **
 **  tPtrListIter::GetP
 **
 **  Similar to Get, but returns a copy of the current data pointer rather
 **  than a pointer to the list item (or 0 if undefined).
 **
\**************************************************************************/
template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::
GetP( int num )
{
  return ( Get( num ) ) ?
    curptrnode->Ptr : 0;
}

template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::GetP( NodeType* nPtr )
{
   return ( Get( nPtr ) ) ?
       curptrnode->Ptr : 0;
}

/**************************************************************************\
 **
 **  tPtrListIter::ReportNextP and ReportPrevP
 **
 **  Returns a copy of the next or previous data pointer without actually
 **  moving to the next or previous item.
 **
\**************************************************************************/
template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::
ReportNextP() const
{
  assert( ptrlistPtr != 0 );
  if( curptrnode == 0 ) return 0;
  if( curptrnode->next == 0 ) return 0;
  return curptrnode->next->Ptr;
}

template< class NodeType >
inline NodeType * tPtrListIter< NodeType >::
ReportPrevP() const
{
  assert( ptrlistPtr != 0 );
  if( curptrnode == 0 ) return 0;
  if( curptrnode->prev == 0 ) return 0;
  return curptrnode->prev->Ptr;
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
template< class NodeType >
inline bool tPtrListIter< NodeType >::
AtEnd() const
{
  if( ptrlistPtr->last == 0 ) return true;
  if( ptrlistPtr->last->next == 0 ) return BOOL( curptrnode==0 );
  return BOOL( curptrnode == ptrlistPtr->first && counter != 0 );
}

#endif


