//-*-c++-*-

/**************************************************************************/
/**
 **  @file tList.h
 **  @brief Header file for classes tList, tListNodeXXX, and tListIter.
 **
 **  A tList is an object that implements a general linked list NodeType
 **  objects, where NodeType can be any type (double, int, other objects,
 **  etc). The one caveat is that tLists are not designed to be lists of
 **  pointers, which have some unique requirements and are thus handled
 **  by tPtrList objects. (Specifically, in a normal list you want to
 **  retrieve the actual item in the list, whereas in a pointer list you
 **  want to retrieve the item to which the list entry points).
 **
 **  Lists can be either linear or circular. The tList class provides
 **  a variety of methods for adding, moving, and retrieving list elements.
 **  For moving back and forth in a tList and retrieving items, it's often
 **  most useful to use a tListIter object (q.v.).
 **
 **  tListNode objects are the nodes on the list; each contains an instance
 **  of the given data type (double, int, class, etc) and a pointer to the
 **  next node in the list.
 **
 **  A tListIter is an iterator for the linked list tList objects (and their
 **  descendants). Its services include fetching data from the current entry
 **  on the list, advancing to the next or previous item on the list, etc.
 **
 **  See also tPtrList, tArray, tMatrix
 **
 **  Changes:
 **    - GT added currentItem member and routines FirstP and NextP to
 **      track position on list w/o an iterator, 1/22/99
 **    - moved all functions into .h file and inlined them (GT 1/20/00)
 **    - AD - March 2004: tListNode is a template argument.
 **
 **  $Id: tList.h,v 1.57 2004-06-16 13:37:33 childcvs Exp $
 */
/**************************************************************************/

#ifndef TLIST_H
#define TLIST_H

#include "tListFwd.h"
#include "../compiler.h"

/**************************************************************************/
/**
 ** @class tListNodeBasic
 **
 ** Class tListNodeBasic represents the items (or "nodes") on the list. Each
 ** tListNodeBasic object has two parts: the data (of type NodeType) and a
 ** pointer to the next item on the list. Capabilities include copy
 ** construction (from either another tListNodeBasic or a NodeType),
 ** returning a pointer or reference to either the data or the tListNodeBasic
 ** itself, and assignment and equality/inequality operations.
 **
 */
/**************************************************************************/
template< class NodeType >
class tListNodeBasic
{
  friend class tList< NodeType, tListNodeBasic< NodeType > >;
  friend class tMeshList< NodeType, tListNodeBasic< NodeType > >;
  friend class tListIter< NodeType, tListNodeBasic< NodeType > >;
  friend class tMeshListIter< NodeType, tListNodeBasic< NodeType > >;
public:
  inline tListNodeBasic();                                // default constructor
  inline tListNodeBasic( const tListNodeBasic< NodeType > & ); // copy constructor #1
  inline tListNodeBasic( const NodeType & );              // copy constructor #2
  ~tListNodeBasic() {next = prev = 0;}
  const tListNodeBasic< NodeType >
  &operator=( const tListNodeBasic< NodeType > & );           // assignment
  inline bool operator==( const tListNodeBasic< NodeType > & ) const; // equality
  inline bool operator!=( const tListNodeBasic< NodeType > & ) const; // inequality
  /*set*/
  inline NodeType getDataNC() const;               // returns copy of data item
  inline NodeType &getDataRefNC();                 // returns modifiable ref to data
  inline NodeType *getDataPtrNC();                 // returns modifiable ptr to data
  inline tListNodeBasic< NodeType > * getNextNC(); // returns ptr to next list node
  /*get*/
  inline NodeType getData() const;                     // returns const copy of data
  inline const NodeType &getDataRef() const;           // returns const ref to data
  inline const NodeType *getDataPtr() const;           // returns const ptr to data
  inline const tListNodeBasic< NodeType > * getNext() const;// returns const ptr to next
  inline const tListNodeBasic< NodeType > * getPrev() const;
  enum {
    isListable = false
  };
  static tListNodeBasic< NodeType > *getListPtr( NodeType const *ptr) {
    return 0;
  }

protected:
  NodeType data_;               // data item
  tListNodeBasic< NodeType > *next; // ptr to next node on list (=0 if end)
  tListNodeBasic< NodeType > *prev;
};


/**************************************************************************\
 **
 **         FUNCTIONS FOR CLASS tListNodeBasic< NodeType >
 **
 **         tListNodeBasics contain a data item and a pointer to the next
 **         tListNodeBasic in the tList (or tMeshList). The data item may
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
 **  tListNodeBasic constructors:
 **
 **  Default constructor: sets next to null
 **  Copy constructor #1: makes a copy of a given tListNodeBasic
 **  Copy constructor #2: fills in data item w/ copy of given NodeType
 **
\**************************************************************************/

//default constructor
template< class NodeType >
inline tListNodeBasic< NodeType >::
tListNodeBasic() :
  next(0),
  prev(0)
{}

//copy constructor with data reference
template< class NodeType >
inline tListNodeBasic< NodeType >::
tListNodeBasic( const tListNodeBasic< NodeType > &original ) :
  data_(original.data_),
  next(original.next),
  prev(original.prev)
{}

//value (by reference) constructor
template< class NodeType >
inline tListNodeBasic< NodeType >::
tListNodeBasic( const NodeType &info ) :
  data_(info),
  next(0),
  prev(0)
{}


/**************************************************************************\
 **
 **  tListNodeBasic overloaded operators:
 **
 **  Assignment: makes a copy (including next ptr)
 **  Equality: compares both data contents and next ptr
 **  Inequality: compares both data contents and next ptr
 **
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >
inline const tListNodeBasic< NodeType > &tListNodeBasic< NodeType >::
operator=( const tListNodeBasic< NodeType > &right )
{
  if( &right != this )
    {
      data_ = right.data_;
      next = right.next;
      prev = right.prev;
    }
  return *this;
}

//overloaded equality operator:
template< class NodeType >
inline bool tListNodeBasic< NodeType >::
operator==( const tListNodeBasic< NodeType > &right ) const
{
  if( next != right.next ) return false;
  if( prev != right.prev ) return false;
  if( &data_ != &(right.data_) ) return false;
  return true;
}

//overloaded inequality operator:
template< class NodeType >
inline bool tListNodeBasic< NodeType >::
operator!=( const tListNodeBasic< NodeType > &right ) const
{
  return ! operator==(right);
}


/**************************************************************************\
 **
 **  tListNodeBasic "get" functions:
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
template< class NodeType >
inline NodeType tListNodeBasic< NodeType >::
getDataNC() const {return data_;}

template< class NodeType >
inline NodeType &tListNodeBasic< NodeType >::
getDataRefNC() {return data_;}

template< class NodeType >
inline NodeType *tListNodeBasic< NodeType >::
getDataPtrNC() {return &data_;}

template< class NodeType >
inline tListNodeBasic< NodeType > *tListNodeBasic< NodeType >::
getNextNC() {return next;}

//return data by value
template< class NodeType >
inline NodeType tListNodeBasic< NodeType >::
getData() const {return data_;}

//return data by reference
template< class NodeType >
inline const NodeType &tListNodeBasic< NodeType >::
getDataRef() const {return data_;}

//return data by pointer
template< class NodeType >
inline const NodeType *tListNodeBasic< NodeType >::
getDataPtr() const {return &data_;}

//return next pointer
template< class NodeType >
inline const tListNodeBasic< NodeType > *tListNodeBasic< NodeType >::
getNext() const {return next;}

//return prev pointer
template< class NodeType >
inline const tListNodeBasic< NodeType > *tListNodeBasic< NodeType >::
getPrev() const {return prev;}


/**************************************************************************/
/**
 ** @class tListable
 **
 ** Used by classes T that can be used by tList<T,tListNodeListable<T> >.
 **
 */
/**************************************************************************/
class tListable
{
public:
  tListable() : listPtr(0) {}
  // listPtr is not copied and set to 0.
  tListable(tListable const &) : listPtr(0) {}
  ~tListable() {listPtr=0;}
  // listPtr is left identical.
  tListable& operator=(tListable const &) { return *this; }
  void setListPtr(void *ptr) { listPtr = ptr; }
  void *getListPtr() const { return listPtr; }
private:
  void *listPtr;  // Pointer to ListNode
};

/**************************************************************************/
/**
 ** @class tListNodeListable
 **
 ** Class tListNodeListable is similar to tListNode except that it sets
 ** a backpointer in NodeType.
 **
 */
/**************************************************************************/
template< class NodeType >
class tListNodeListable
{
  friend class tList< NodeType, tListNodeListable< NodeType > >;
  friend class tMeshList< NodeType, tListNodeListable< NodeType > >;
  friend class tListIter< NodeType, tListNodeListable< NodeType > >;
  friend class tMeshListIter< NodeType, tListNodeListable< NodeType > >;
public:
  inline tListNodeListable();                                // default constructor
  inline tListNodeListable( const tListNodeListable< NodeType > & ); // copy constructor #1
  inline tListNodeListable( const NodeType & );              // copy constructor #2
  ~tListNodeListable() {next = prev = 0;}
  const tListNodeListable< NodeType >
  &operator=( const tListNodeListable< NodeType > & );           // assignment
  inline bool operator==( const tListNodeListable< NodeType > & ) const; // equality
  inline bool operator!=( const tListNodeListable< NodeType > & ) const; // inequality
  /*set*/
  inline NodeType getDataNC() const;               // returns copy of data item
  inline NodeType &getDataRefNC();                 // returns modifiable ref to data
  inline NodeType *getDataPtrNC();                 // returns modifiable ptr to data
  inline tListNodeListable< NodeType > * getNextNC();// returns ptr to next list node
  /*get*/
  inline NodeType getData() const;                     // returns const copy of data
  inline const NodeType &getDataRef() const;           // returns const ref to data
  inline const NodeType *getDataPtr() const;           // returns const ptr to data
  inline const tListNodeListable< NodeType > * getNext() const;// returns const ptr to next
  inline const tListNodeListable< NodeType > * getPrev() const;
  enum {
    isListable = true
  };
  static tListNodeListable< NodeType > *getListPtr( NodeType const *dataPtr) {
    return static_cast<tListNodeListable< NodeType >*>(dataPtr->getListPtr());
  }

protected:
  NodeType data_;               // data item
  tListNodeListable< NodeType > *next; // ptr to next node on list (=0 if end)
  tListNodeListable< NodeType > *prev;
};


/**************************************************************************\
 **
 **         FUNCTIONS FOR CLASS tListNodeListable< NodeType >
 **
 **         tListNodeListables contain a data item and a pointer to the next
 **         tListNodeListable in the tList (or tMeshList). The data item may
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
 **  tListNodeListable constructors:
 **
 **  Default constructor: sets next to null
 **  Copy constructor #1: makes a copy of a given tListNodeListable
 **  Copy constructor #2: fills in data item w/ copy of given NodeType
 **
\**************************************************************************/

//default constructor
template< class NodeType >
inline tListNodeListable< NodeType >::
tListNodeListable() :
  next(0),
  prev(0)
{
  data_.setListPtr(this);
}

//copy constructor with data reference
template< class NodeType >
inline tListNodeListable< NodeType >::
tListNodeListable( const tListNodeListable< NodeType > &original ) :
  data_(original.data_),
  next(original.next),
  prev(original.prev)
{
  data_.setListPtr(this);
}

//value (by reference) constructor
template< class NodeType >
inline tListNodeListable< NodeType >::
tListNodeListable( const NodeType &info ) :
  data_(info),
  next(0),
  prev(0)
{
  data_.setListPtr(this);
}


/**************************************************************************\
 **
 **  tListNodeListable overloaded operators:
 **
 **  Assignment: makes a copy (including next ptr)
 **  Equality: compares both data contents and next ptr
 **  Inequality: compares both data contents and next ptr
 **
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >
inline const tListNodeListable< NodeType > &tListNodeListable< NodeType >::
operator=( const tListNodeListable< NodeType > &right )
{
  if( &right != this )
    {
      data_ = right.data_;
      next = right.next;
      prev = right.prev;
    }
  return *this;
}

//overloaded equality operator:
template< class NodeType >
inline bool tListNodeListable< NodeType >::
operator==( const tListNodeListable< NodeType > &right ) const
{
  if( next != right.next ) return false;
  if( prev != right.prev ) return false;
  if( &data_ != &(right.data_) ) return false;
  return true;
}

//overloaded inequality operator:
template< class NodeType >
inline bool tListNodeListable< NodeType >::
operator!=( const tListNodeListable< NodeType > &right ) const
{
  return ! operator==(right);
}


/**************************************************************************\
 **
 **  tListNodeListable "get" functions:
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
template< class NodeType >
inline NodeType tListNodeListable< NodeType >::
getDataNC() const {
  NodeType aData(data_);
  aData.setListPtr(0);
  return aData;
}

template< class NodeType >
inline NodeType &tListNodeListable< NodeType >::
getDataRefNC() {return data_;}

template< class NodeType >
inline NodeType *tListNodeListable< NodeType >::
getDataPtrNC() {return &data_;}

template< class NodeType >
inline tListNodeListable< NodeType > *tListNodeListable< NodeType >::
getNextNC() {return next;}

//return data by value
template< class NodeType >
inline NodeType tListNodeListable< NodeType >::
getData() const {
  NodeType aData(data_);
  aData.setListPtr(0);
  return aData;
}

//return data by reference
template< class NodeType >
inline const NodeType &tListNodeListable< NodeType >::
getDataRef() const {return data_;}

//return data by pointer
template< class NodeType >
inline const NodeType *tListNodeListable< NodeType >::
getDataPtr() const {return &data_;}

//return next pointer
template< class NodeType >
inline const tListNodeListable< NodeType > *tListNodeListable< NodeType >::
getNext() const {return next;}

//return prev pointer
template< class NodeType >
inline const tListNodeListable< NodeType > *tListNodeListable< NodeType >::
getPrev() const {return prev;}


/**************************************************************************/
/**
 ** @class tList
 **
 ** Class tList implements a linked list. The class includes pointers to
 ** the first, last, and current list nodes (arg "ListNodeType").
 **
 */
/**************************************************************************/
template< class NodeType, class ListNodeType >
class tList
{
  friend class tListIter< NodeType, ListNodeType >;
  friend class tMeshListIter< NodeType, ListNodeType >;
public:
  inline tList();                     // default constructor
  tList(const tList&);
  tList( const tList< NodeType, ListNodeType > * ); // copy constructor
  ~tList();                           // destructor
  const tList< NodeType, ListNodeType >
  &operator=( const tList< NodeType, ListNodeType > & );           // assignment
  inline bool operator==( const tList< NodeType, ListNodeType > & ) const; // equality
  inline bool operator!=( const tList< NodeType, ListNodeType > & ) const; // inequality
  inline void insertAtFront( const NodeType & ); // puts copy of item at list front
  inline void insertAtBack( const NodeType & );  // puts copy of item at list back
  inline void insertAtNext( const NodeType &, ListNodeType * );
  inline void insertAtPrev( const NodeType &, ListNodeType * );
  inline int removeFromFront( NodeType & ); // removes 1st item, puts it in ref
  inline int removeFromBack( NodeType & );  // removes last item, puts it in ref
  inline int removeNext( NodeType &, ListNodeType * );
  inline int removePrev( NodeType &, ListNodeType * );
  void Flush();        // clears and reinitializes list
  inline bool isEmpty() const; // returns true is list is empty, false otherwise
#ifndef NDEBUG
  void print() const;  // prints contents of list - DEBUG ONLY
#endif
  inline int getSize() const; // returns # of items on list
  inline ListNodeType const * getFirst() const; // returns ptr to 1st list node
  inline ListNodeType const * getLast() const;  // returns ptr to last list node
  inline ListNodeType * getFirstNC(); // returns ptr to 1st list node
  inline ListNodeType * getLastNC();  // returns ptr to last list node
  inline NodeType * FirstP();  // returns ptr to 1st data item & sets current to 1st
  inline NodeType * NextP();   // moves to next node and returns ptr to data item
  inline void moveToBack( ListNodeType *  );  // move given node to back
  inline void moveToFront( ListNodeType *  ); // move given node to front
  inline void moveToBefore( ListNodeType *, ListNodeType * );
  inline void moveToAfter( ListNodeType *, ListNodeType * );
  inline void makeCircular();   // makes list circular (last points to first)
  inline const NodeType getIthData( int ) const;     // rtns copy of given item #
  inline const NodeType *getIthDataPtr( int ) const; // rtns ptr to given item #
  inline const NodeType &getIthDataRef( int ) const; // rtns ref to given item #
  inline NodeType getIthDataNC( int ) const;     // rtns modifiable copy of item #
  inline NodeType *getIthDataPtrNC( int ); // rtns modifiable ptr to item #
  inline NodeType &getIthDataRefNC( int ); // rtns modifiable ref to item #
  ListNodeType * getListNode( const NodeType * ); // rtns ptr to node #
  const ListNodeType * getIthListNode( int ) const;
  ListNodeType * getIthListNodeNC( int );

#ifndef NDEBUG
  void DebugTellPtrs() const;
#endif

protected:
  int nNodes;                          // # of items on list
  ListNodeType * first;       // ptr to first node
  ListNodeType * last;        // ptr to last node
  ListNodeType * currentItem; // ptr to current item
  ListNodeType * getNewNode( const NodeType & ); // makes new node
};



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
template< class NodeType, class ListNodeType >
inline tList< NodeType, ListNodeType >::tList() :
  nNodes(0), first(0), last(0), currentItem(0)
{}

//copy constructor
template< class NodeType, class ListNodeType >
tList< NodeType, ListNodeType >::
tList( const tList< NodeType, ListNodeType > *original ) :
  nNodes(0), first(0), last(0), currentItem(0)
{
  int i;

  assert( original != 0 );
  ListNodeType * current = original->first;
  for( i=0; i<original->nNodes; ++i )
    {
      insertAtBack( current->getDataRef() );
      current = current->next;
    }
  assert( nNodes == original->nNodes );
  if (0) //DEBUG
    std::cout << "list copy instantiated" << first << std::endl;
  current = first;

}

//copy constructor
template< class NodeType, class ListNodeType >
tList< NodeType, ListNodeType >::
tList( const tList< NodeType, ListNodeType > &original ) :
  nNodes(0), first(0), last(0), currentItem(0)
{
  int i;

  assert( original != 0 );
  ListNodeType * current = original.first;
  for( i=0; i<original.nNodes; ++i )
    {
      insertAtBack( current->getDataRef() );
      current = current->next;
    }
  assert( nNodes == original.nNodes );
  if (0) //DEBUG
    std::cout << "list copy instantiated" << first << std::endl;
  current = first;

}

//destructor
template< class NodeType, class ListNodeType >
tList< NodeType, ListNodeType >::
~tList()
{
  if( !isEmpty() )
    {
      if (0) //DEBUG
	std::cout<<"Destroying nodes ... "<<std::endl;
      ListNodeType * current = first, * temp;
      while( current != 0 )
	{
	  temp = current;
	  current = current->next;
	  delete temp;
	}
    }
  first = 0;
  last = 0;
  currentItem = 0;
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
template< class NodeType, class ListNodeType >
const tList< NodeType, ListNodeType > &tList< NodeType, ListNodeType >::
operator=( const tList< NodeType, ListNodeType > &right )
{
  if( this != &right )
    {
      // create an equivalent object rather than pointers to the original
      Flush();
      ListNodeType *cn = right.first;
      if( cn != 0 )
	{
          insertAtBack( cn->getDataRef() );
          for( cn = cn->next; cn != last->next; cn = cn->next )
	    insertAtBack( cn->getDataRef() );
	  if( right.last->next == right.first ) last->next = first;
	  if( right.first->prev == right.last ) first->prev = last;
	}
      assert( nNodes == right.nNodes );
    }
  return *this;
}

//overloaded equality operator:
template< class NodeType, class ListNodeType >
inline bool tList< NodeType, ListNodeType >::
operator==( const tList< NodeType, ListNodeType > &right ) const
{
  if( nNodes != right.nNodes ) return false;
  if( first != right.first ) return false;
  if( last != right.last ) return false;
  return true;
}

//overloaded inequality operator:
template< class NodeType, class ListNodeType >
inline bool tList< NodeType, ListNodeType >::
operator!=( const tList< NodeType, ListNodeType > &right ) const
{
  return ! operator==(right);
}


/**************************************************************************\
 **
 **  tList::getNewNode
 **
 **  Creates a new ListNodeType and returns a pointer to it. Used by list
 **  insertion routines (see below); not publically accessible.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline ListNodeType * tList< NodeType, ListNodeType >::
getNewNode( const NodeType &value )
{
  ListNodeType * ptr =
    new ListNodeType( value );
  assert( ptr != 0 );
  ++nNodes;
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
template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
insertAtFront( const NodeType &value )
{
  if (0) //DEBUG
    std::cout << "ADD NEW NODE TO LIST AT FRONT" << std::endl;

  ListNodeType *newPtr = getNewNode( value );
  if( isEmpty() )
    first = last = currentItem = newPtr;
  else
    {
      newPtr->next = first;
      newPtr->prev = first->prev;
      first->prev = newPtr;
      if( last->next == first ) last->next = newPtr;
      first = newPtr;
    }
}

//insert at back
template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
insertAtBack( const NodeType &value )
{
  ListNodeType * newPtr = getNewNode( value );
  if( isEmpty() )
    {
      first = last = currentItem = newPtr;
    }
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

//insert at next spot in list
template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
insertAtNext( const NodeType &value, ListNodeType * prev )
{
  if( prev != 0 )
    {
      if( prev == last )
	{
	  insertAtBack( value );
	  return;
	}
      ListNodeType * newPtr = getNewNode( value );
      newPtr->next = prev->next;
      newPtr->prev = prev;
      prev->next->prev = newPtr;
      prev->next = newPtr;
    }
}

//insert at previous spot in list
template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
insertAtPrev( const NodeType &value, ListNodeType * node )
{
  if( node != 0 )
    {
      if( node == first )
	{
	  insertAtFront( value );
	  return;
	}
      ListNodeType * newPtr = getNewNode( value );
      newPtr->next = node;
      newPtr->prev = node->prev;
      node->prev->next = newPtr;
      node->prev = newPtr;
    }
}

//delete from front
template< class NodeType, class ListNodeType >
inline int tList< NodeType, ListNodeType >::
removeFromFront( NodeType &value )
{
  if( isEmpty() ) return 0;
  ListNodeType * temp = first;
  if( first == last ) first = last = currentItem = 0;
  else
    {
      first->next->prev = first->prev;
      if( last->next == first ) last->next = first->next;
      if( currentItem==first ) currentItem = first->next;
      first = first->next;
    }
  value = temp->getData();
  delete temp;
  --nNodes;
  return 1;
}

//delete from back
template< class NodeType, class ListNodeType >
inline int tList< NodeType, ListNodeType >::
removeFromBack( NodeType &value )
{
  if( isEmpty() ) return 0;
  ListNodeType * temp = last;
  if( first == last ) first = last = currentItem = 0;
  else
    {
      last->prev->next = last->next;
      if( first->prev == last ) first->prev = last->prev;
      if( currentItem==last ) currentItem = last->prev;
      last = last->prev;
    }
  value = temp->getData();
  delete temp;
  --nNodes;
  return 1;
}

//delete next node
template< class NodeType, class ListNodeType >
inline int tList< NodeType, ListNodeType >::
removeNext( NodeType &value, ListNodeType * ptr )
{
  if( ptr == 0 ) return 0;
  if( ptr->next == 0 ) return 0;
  if( ptr->next == last ) return removeFromBack( value );
  if( ptr->next == first ) return removeFromFront( value );
  ListNodeType * temp = ptr->next;
  if( currentItem == temp ) currentItem = ptr->next->next;
  ptr->next = ptr->next->next;
  ptr->next->prev = ptr;
  value = temp->getData();
  delete temp;
  --nNodes;
  return 1;
}

//delete previous node
template< class NodeType, class ListNodeType >
inline int tList< NodeType, ListNodeType >::
removePrev( NodeType &value, ListNodeType * ptr )
{
  if( ptr == 0 ) return 0;
  if( ptr->prev == 0 ) return 0;
  if( ptr->prev == last ) return removeFromBack( value );
  if( ptr->prev == first ) return removeFromFront( value );
  ListNodeType* temp = ptr->prev;
  if( currentItem == ptr->prev ) currentItem = ptr->prev->prev;
  ptr->prev = ptr->prev->prev;
  ptr->prev->next = ptr;
  value = temp->getData();
  delete temp;
  --nNodes;
  return 1;
}


/**************************************************************************\
 **
 **  tList::Flush
 **
 **  Deletes all nodes on list.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
void tList< NodeType, ListNodeType >::
Flush()
{
  ListNodeType * current = first, * temp;

  if( !isEmpty() )
    {
      while( current != 0 )
	{
	  temp = current;
	  current = current->next;
	  delete temp;
	}
      first = last = currentItem = 0;
    }
  assert( isEmpty() );
  nNodes = 0;
}


/**************************************************************************\
 **
 **  tList::isEmpty
 **
 **  Returns TRUE if first points to null; FALSE otherwise.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline bool tList< NodeType, ListNodeType >::
isEmpty() const
{
  return BOOL( first == 0 );
}

//display list contents -- for debugging only

template< class NodeType, class ListNodeType >
void tList< NodeType, ListNodeType >::
print() const
{
  if( isEmpty() )
    {
      std::cout<<"The list is empty\n"<<std::endl;
      return;
    }
  ListNodeType * current = first;
  std::cout<<"The list is: ";
  while( current != 0 )
    {
      current = current->next;
    }
  std::cout<< '\n' <<std::endl;
}


/**************************************************************************\
 **
 **  tList "get" functions:
 **
 **  getSize: returns # of items on list
 **  getFirst: returns const ptr to first ListNodeType
 **  getLast: returns const ptr to last ListNodeType
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

template< class NodeType, class ListNodeType >
inline int tList< NodeType, ListNodeType >::
getSize() const {return nNodes;}

template< class NodeType, class ListNodeType >
inline ListNodeType const * tList< NodeType, ListNodeType >::
getFirst() const {return first;}

template< class NodeType, class ListNodeType >
inline ListNodeType const * tList< NodeType, ListNodeType >::
getLast() const {return last;}

template< class NodeType, class ListNodeType >
inline ListNodeType * tList< NodeType, ListNodeType >::
getFirstNC() {return first;}

template< class NodeType, class ListNodeType >
inline ListNodeType * tList< NodeType, ListNodeType >::
getLastNC() {return last;}

// Added by gt 1/99
template< class NodeType, class ListNodeType >
inline NodeType * tList< NodeType, ListNodeType >::
FirstP()
{
  assert( first!=0 );
  currentItem = first;
  return first->getDataPtrNC();
}

// Added by gt 1/99
template< class NodeType, class ListNodeType >
inline NodeType * tList< NodeType, ListNodeType >::
NextP()
{
  assert( currentItem!=0 );
  currentItem = currentItem->next;
  if( currentItem!=0 )
    return currentItem->getDataPtrNC();
  else return 0;
}

template< class NodeType, class ListNodeType >
inline const NodeType tList< NodeType, ListNodeType >::
getIthData( int num ) const
{
  assert( num >= 0 && num < nNodes );
  return getIthListNode(num)->getData();
}

template< class NodeType, class ListNodeType >
inline const NodeType &tList< NodeType, ListNodeType >::
getIthDataRef( int num ) const
{
  assert( num >= 0 && num < nNodes );
  return getIthListNode(num)->getDataRef();
}

template< class NodeType, class ListNodeType >
inline const NodeType *tList< NodeType, ListNodeType >::
getIthDataPtr( int num ) const
{
  assert( num >= 0 && num < nNodes );
  return getIthListNode(num)->getDataPtr();
}

//set
template< class NodeType, class ListNodeType >
inline NodeType tList< NodeType, ListNodeType >::
getIthDataNC( int num ) const
{
  assert( num >= 0 && num < nNodes );
  return getIthListNode(num)->getData();
}

template< class NodeType, class ListNodeType >
inline NodeType &tList< NodeType, ListNodeType >::
getIthDataRefNC( int num )
{
  assert( num >= 0 && num < nNodes );
  return getIthListNodeNC(num)->getDataRefNC();
}

template< class NodeType, class ListNodeType >
inline NodeType *tList< NodeType, ListNodeType >::
getIthDataPtrNC( int num )
{
  assert( num >= 0 && num < nNodes );
  return getIthListNodeNC(num)->getDataPtrNC();
}


/**************************************************************************\
 **
 **  tList::getListNode
 **
 **  Finds and returns the list node containing the data pointed to by
 **  desiredDatPtr, or zero if not found.
 **
 **  Parameters:  desiredDatPtr -- pointer to the data item sought after
 **  Returns:  pointer to the ListNodeType containing desiredDatPtr, or zero
 **            if not found
 **  Notes: might be safer to implement with a const return type
 **  Created: 4/29/98 GT
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
ListNodeType * tList< NodeType, ListNodeType >::
getListNode( const NodeType * desiredDatPtr )
{
  ListNodeType * listnode = first;

  if( listnode==0 ) return 0;
  while( listnode->getDataPtr() != desiredDatPtr )
    {
      listnode = listnode->next;
      if( listnode==0 ) return 0;
    }
  return listnode;

}

/**************************************************************************\
 **
 **  tList::getIthListNode( i )
 **
 **  Finds and returns the list node indicated by i [0,nnodes-1),
 **  or zero if not found.
 **
 **  Parameters:  i -- number of item sought after
 **  Returns:  pointer to the ith ListNodeType, or zero
 **            if not found
 **  Notes: might be safer to implement with a const return type
 **  Created: 4/29/98 GT
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
const ListNodeType* tList< NodeType, ListNodeType >::
getIthListNode( int num ) const
{
  if( num < 0 || num >= nNodes )
    return 0;
  int i;
  ListNodeType const * curPtr;
  if( num > nNodes / 2 )
    for( curPtr = last, i = nNodes-1; i>num; curPtr = curPtr->prev, --i ) ;
  else
    for( curPtr = first, i = 0; i<num; curPtr = curPtr->next, ++i ) ;
  return curPtr;
}

template< class NodeType, class ListNodeType >
ListNodeType* tList< NodeType, ListNodeType >::
getIthListNodeNC( int num )
{
  if( num < 0 || num >= nNodes )
    return 0;
  int i;
  ListNodeType * curPtr;
  if( num > nNodes / 2 )
    for( curPtr = last, i = nNodes-1; i>num; curPtr = curPtr->prev, --i ) ;
  else
    for( curPtr = first, i = 0; i<num; curPtr = curPtr->next, ++i ) ;
  return curPtr;
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

template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
moveToBack( ListNodeType * mvnode )
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

template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
moveToFront( ListNodeType * mvnode )
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

template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
moveToBefore( ListNodeType* mvnode,
              ListNodeType* plcnode )
{
  if( mvnode == plcnode ) return;
  if( plcnode == first )
    {
      moveToFront( mvnode );
      return;
    }
  if( mvnode == last )
    {
      last->prev->next = last->next;
      last = last->prev;
    }
  else if( mvnode == first )
    {
      first->next->prev = first->prev;
      first = first->next;
    }
  else
    {
      mvnode->prev->next = mvnode->next;
      mvnode->next->prev = mvnode->prev;
    }
  mvnode->next = plcnode;
  mvnode->prev = plcnode->prev;
  plcnode->prev->next = mvnode;
  plcnode->prev = mvnode;
}

template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
moveToAfter( ListNodeType* mvnode,
             ListNodeType* plcnode )
{
  if( mvnode == plcnode ) return;
  if( plcnode->next == mvnode ) return;
  if( plcnode == last )
    {
      moveToBack( mvnode );
      return;
    }
  if( mvnode == last )
    {
      last->prev->next = last->next;
      last = last->prev;
    }
  else if( mvnode == first )
    {
      first->next->prev = first->prev;
      first = first->next;
    }
  else
    {
      mvnode->prev->next = mvnode->next;
      mvnode->next->prev = mvnode->prev;
    }
  mvnode->next = plcnode->next;
  mvnode->prev = plcnode;
  plcnode->next->prev = mvnode;
  plcnode->next = mvnode;
}

/**************************************************************************\
 **
 **  tList::makeCircular
 **
 **  Converts the list into a circular list by having the last item point
 **  to the first.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline void tList< NodeType, ListNodeType >::
makeCircular()
{
  assert( first != 0 && last != 0 );
  last->next = first;
  first->prev = last;
}


template< class NodeType, class ListNodeType >
void tList< NodeType, ListNodeType >::
DebugTellPtrs() const
{
  std::cout << "first: " << first << std::endl;
  std::cout << "last: " << last << std::endl;
}



/**************************************************************************/
/**
 ** @class tListIter
 **
 ** Helper class for tList. tListIters are "iterators" that walk up and
 ** down a tList, fetching items. Their chief advantage is that you can
 ** have multiple iterators on any given list at once, and thus multiple
 ** access points. Use of iterator classes is discussed by Deitel and
 ** Deitel, _C++ How to Program_, first edition, Prentice Hall, 1994.
 **
 ** Note that in the current implementation, list items are fetched by
 ** ID number, which presupposes that the list items have a member function
 ** getID. This restricts the generality of tList, and should be moved
 ** to tMeshList. (TODO)
 **
 */
/**************************************************************************/
//TO DO: make Get, Where, GetP, refer to place in list rather than use getID()
template< class NodeType, class ListNodeType >
class tListIter
{
  tListIter& operator=(const tListIter&);
  int NextIfNoCurrent();  // set 1st as current undefined
  int PrevIfNoCurrent();  // set last as current undefined
  void GetByPtrVerify( NodeType const *, ListNodeType const *);
public:
  inline tListIter();               // default constructor
  tListIter(const tListIter&);      // copy constructor
  tListIter( tList< NodeType, ListNodeType > & ); // constructor: reference to list
  tListIter( tList< NodeType, ListNodeType > * ); // constructor: ptr to list
  inline ~tListIter();   // destructor
  inline int First();    // sets position to 1st list item (rtns 0 on failure)
  inline int Last();     // sets position to last "    "     "
  int Get( int ); // use only if NodeType has member getID()!!
  int GetByPermID( int );  // use only if NodeType has member getPermID()!!
  inline int GetByPtr( NodeType const * );
  int GetByPtrSlow( NodeType const * );
  inline NodeType * GetByPtrP( NodeType const * );
  inline int Next();     // advances to next item (or 1st if current undefined)
  inline int Prev();     // moves to previous item (or last if current undef'd)
  inline int Where() const;  // use only if NodeType has member getID()!!
  inline bool AtEnd() const;  // returns true if at end of the list
  inline NodeType &DatRef();  // returns ref to current data item
  inline NodeType *DatPtr();  // returns ptr to current data item
  inline ListNodeType *NodePtr();  // returns ptr to current list node
  inline void Reset( tList< NodeType, ListNodeType > & ); // tells iterator to work on given list
  inline NodeType * FirstP();  // moves to 1st item and rtns ptr to it
  inline NodeType * LastP();   // moves to last  "   "
  inline NodeType * NextP();   // moves to next  "   "
  inline NodeType * PrevP();   // moves to previous " "
  inline NodeType * GetP( int num ); //use only if NodeType has member getID()!!
  inline NodeType * GetPByPermID( int num ); //use only if NodeType has member getPermID()!!

protected:
  ListNodeType * curnode;  // ptr to current list node
  tList< NodeType, ListNodeType > *listPtr;       // ptr to current list
  int counter;                      // current position on list (first=0)
};

/**************************************************************************\
 **
 **         FUNCTIONS FOR CLASS tListIter
 **
 **  A tListIter is an iterator for the linked list tList objects (and their
 **  descendants). Its services include fetching data from the current entry
 **  on the list, advancing to the next or previous item on the list, etc.
 **
 **  See also tMeshList.
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

template< class NodeType, class ListNodeType >
inline tListIter< NodeType, ListNodeType >::
tListIter() :
  curnode(0),
  listPtr(0),
  counter(0)
{
  if (0) //DEBUG
    std::cout << "tListIter()" << std::endl;
}

template< class NodeType, class ListNodeType >
inline tListIter< NodeType, ListNodeType >::
tListIter(const tListIter& c) :
  curnode(c.curnode),
  listPtr(c.listPtr),
  counter(c.counter)
{
  if (0) //DEBUG
    std::cout << "tListIter(const tListIter&)" << std::endl;
}

template< class NodeType, class ListNodeType >
inline tListIter< NodeType, ListNodeType >::
tListIter( tList< NodeType, ListNodeType > &list ) :
  curnode(list.first),
  listPtr(&list),
  counter(0)
{
  //assert( curnode != 0 );
  if (0) //DEBUG
    std::cout << "tListIter( list )" << std::endl;
}

template< class NodeType, class ListNodeType >
inline tListIter< NodeType, ListNodeType >::
tListIter( tList< NodeType, ListNodeType > *ptr ) :
  curnode(0),
  listPtr(ptr),
  counter(0)
{
  assert( ptr != 0 );
  curnode = ptr->first;
}

template< class NodeType, class ListNodeType >
inline tListIter< NodeType, ListNodeType >::
~tListIter()
{
  listPtr = 0;
  curnode = 0;
  if (0) //DEBUG
    std::cout << "~tListIter()" << std::endl;
}


/**************************************************************************\
 **
 **  tListIter::Reset
 **
 **  Points iterator at the 1st node on _list_ (provides a way of telling
 **  an iterator which list to work on).
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline void tListIter< NodeType, ListNodeType >::
Reset( tList< NodeType, ListNodeType > &list )
{
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
template< class NodeType, class ListNodeType >
inline int tListIter< NodeType, ListNodeType >::
First()
{
  assert( listPtr != 0 );
  curnode = listPtr->first;
  counter = 0;
  if( curnode != 0 ) return 1;
  return ( curnode == 0 && listPtr->isEmpty() ) ? 1:0;
}

template< class NodeType, class ListNodeType >
inline int tListIter< NodeType, ListNodeType >::
Last()
{
  assert( listPtr != 0 );
  curnode = listPtr->last;
  counter = -1;
  return ( curnode != 0 ) ? 1:0;
}


/**************************************************************************\
 **
 **  tListIter::Get
 **
 **  Move to list item with ID number _num_. Note: assumes that list items
 **  have a member function getID()! Returns 1 if found, 0 if not.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
int tListIter< NodeType, ListNodeType >::
Get( int num )
{
  assert( listPtr != 0 );
  if( num < 0 ) std::cout << "tListIter::Get(num): num < 0" << std::endl;
  ListNodeType *tempnodeptr;
  for( tempnodeptr = listPtr->first, counter = 0;
       counter <= listPtr->nNodes && tempnodeptr != 0;
       tempnodeptr = tempnodeptr->next, ++counter )
    {
      if( tempnodeptr->getDataPtr()->getID() == num ) break;
    }
  if( tempnodeptr == 0 ) return 0;
  if( tempnodeptr->getDataPtr()->getID() != num ) return 0;
  curnode = tempnodeptr;
  return 1;
}


/**************************************************************************\
 **
 **  tListIter::GetByPermID
 **
 **  Move to list item with permID number _num_. Note: assumes that list items
 **  have a member function getPermID()! Returns 1 if found, 0 if not.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
int tListIter< NodeType, ListNodeType >::
GetByPermID( int num )
{
  assert( listPtr != 0 );
  if( num < 0 ) std::cout << "tListIter::Get(num): num < 0" << std::endl;
  ListNodeType *tempnodeptr;
  for( tempnodeptr = listPtr->first, counter = 0;
       counter <= listPtr->nNodes && tempnodeptr != 0;
       tempnodeptr = tempnodeptr->next, ++counter )
    {
      if( tempnodeptr->getDataPtr()->getPermID() == num ) break;
    }
  if( tempnodeptr == 0 ) return 0;
  if( tempnodeptr->getDataPtr()->getPermID() != num ) return 0;
  curnode = tempnodeptr;
  return 1;
}


/**************************************************************************\
 **
 **  tListIter::GetByPtr
 **
 **  Move to list item with the same adress.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline int tListIter< NodeType, ListNodeType >::
GetByPtr( NodeType const *dataPtr )
{
  assert( listPtr != 0 );
  assert( dataPtr != 0 );

  if (ListNodeType::isListable) {
    // ListNodeType is "listable". Therefore, we use the back pointer.
    ListNodeType *tempnodeptr = ListNodeType::getListPtr( dataPtr );
    if (0) //DEBUG
      GetByPtrVerify( dataPtr, tempnodeptr);
    if (tempnodeptr == 0)
      return 0;
    curnode = tempnodeptr;
    counter = -2;
    return 1;
  } else
    // linear search.
    return GetByPtrSlow( dataPtr );
}

template< class NodeType, class ListNodeType >
void tListIter< NodeType, ListNodeType >::
GetByPtrVerify( NodeType const *dataPtr, ListNodeType const *tempnodeptr )
{
  if (Get(dataPtr->getID()) != 0)
    assert(tempnodeptr == curnode);
  else
    assert(tempnodeptr == 0);
}

template< class NodeType, class ListNodeType >
int tListIter< NodeType, ListNodeType >::
GetByPtrSlow( NodeType const *dataPtr )
{
  assert( listPtr != 0 );
  if( dataPtr == NULL )
    std::cout << "tListIter::GetByPtr(ptr): ptr < 0" << std::endl;

  // linear search.
  ListNodeType *tempnodeptr;
  for( tempnodeptr = listPtr->first, counter = 0;
       counter <= listPtr->nNodes && tempnodeptr != 0;
       tempnodeptr = tempnodeptr->next, ++counter )
    {
      if( tempnodeptr->getDataPtr() == dataPtr ) break;
    }
  if( tempnodeptr == 0 ) return 0;
  if( tempnodeptr->getDataPtr() != dataPtr ) return 0;
  curnode = tempnodeptr;
  return 1;
}

/**************************************************************************\
 **
 **  tListIter::GetByPtrP
 **
 **  Similar to Get, but returns a pointer to the current data item (or
 **  0 if undefined).
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline NodeType * tListIter< NodeType, ListNodeType >::
GetByPtrP( NodeType const *dataPtr )
{
  return ( GetByPtr( dataPtr ) ) ?
    curnode->getDataPtrNC() : 0;
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
template< class NodeType, class ListNodeType >
inline int tListIter< NodeType, ListNodeType >::
Next()
{
  if( likely(curnode != 0) )
    {
      curnode = curnode->next;
      ++counter;
      return curnode != 0 ? 1:0;
    }
  // unlikely case
  return NextIfNoCurrent();
}

template< class NodeType, class ListNodeType >
int tListIter< NodeType, ListNodeType >::
NextIfNoCurrent()
{
  assert( listPtr != 0 );
  assert( curnode == 0 );
  curnode = listPtr->first;
  counter = 0;
  return curnode != 0 ? 1:0;
}

template< class NodeType, class ListNodeType >
inline int tListIter< NodeType, ListNodeType >::
Prev()
{
  if( likely(curnode != 0) )
    {
      curnode = curnode->prev;
      --counter;
      return curnode != 0 ? 1:0;
    }
  // unlikely case
  return PrevIfNoCurrent();
}

template< class NodeType, class ListNodeType >
int tListIter< NodeType, ListNodeType >::
PrevIfNoCurrent()
{
  assert( listPtr != 0 );
  assert( curnode == 0 );
  curnode = listPtr->last;
  counter = -1;
  return curnode != 0 ? 1:0;
}

/**************************************************************************\
 **
 **  tListIter::FirstP and tListIter::LastP
 **
 **  Move to the first or last item on the list and return a pointer to the
 **  data, or 0 if first/last item is empty.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline NodeType * tListIter< NodeType, ListNodeType >::
FirstP()
{
  return ( First() ) ?
    curnode->getDataPtrNC() : 0;
}

template< class NodeType, class ListNodeType >
inline NodeType * tListIter< NodeType, ListNodeType >::
LastP()
{
  return ( Last() ) ?
    curnode->getDataPtrNC() : 0;
}

/**************************************************************************\
 **
 **  tListIter::NextP and tListIter::PrevP
 **
 **  Same as Next and Prev, except that the functions return a pointer to
 **  the current data item (or 0 if none exists).
 **
\**************************************************************************/

template< class NodeType, class ListNodeType >
inline NodeType * tListIter< NodeType, ListNodeType >::
NextP()
{
  return ( Next() ) ?
    curnode->getDataPtrNC() : 0;
}

template< class NodeType, class ListNodeType >
inline NodeType *tListIter< NodeType, ListNodeType >::
PrevP()
{
  return ( Prev() ) ?
    curnode->getDataPtrNC() : 0;
}


/**************************************************************************\
 **
 **  tListIter::GetP
 **
 **  Similar to Get, but returns a pointer to the current data item (or
 **  0 if undefined).
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline NodeType * tListIter< NodeType, ListNodeType >::
GetP( int num )
{
  return ( Get( num ) ) ?
    curnode->getDataPtrNC() : 0;
}


/**************************************************************************\
 **
 **  tListIter::GetPByPermID
 **
 **  Similar to GetByPermID, but returns a pointer to the current data item (or
 **  0 if undefined).
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline NodeType * tListIter< NodeType, ListNodeType >::
GetPByPermID( int num )
{
  return ( GetByPermID( num ) ) ?
    curnode->getDataPtrNC() : 0;
}


/**************************************************************************\
 **
 **  tListIter::Where
 **
 **  Returns the ID number of the current data item, or -1 if there is
 **  no current data item. Assumes data item has a getID() mbr function!
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline int tListIter< NodeType, ListNodeType >::
Where() const
{
  return ( curnode == 0 ) ?
    -1 : curnode->getDataPtr()->getID();
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
template< class NodeType, class ListNodeType >
inline bool tListIter< NodeType, ListNodeType >::
AtEnd() const
{
  if( listPtr->isEmpty() ) return true;
  if( listPtr->last->next == 0 ) return BOOL( curnode==0 );
  return BOOL( curnode == listPtr->first && counter != 0 );
}


/**************************************************************************\
 **
 **  tListIter::DatRef, DatPtr, and NodePtr
 **
 **  Returns a non-constant reference or pointer to the current data item
 **  or the current list node.
 **
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline NodeType &tListIter< NodeType, ListNodeType >::
DatRef()
{
  return curnode->getDataRefNC();
}

template< class NodeType, class ListNodeType >
inline NodeType *tListIter< NodeType, ListNodeType >::
DatPtr()
{
  return ( curnode == 0 ) ?
    0 : curnode->getDataPtrNC();
}

template< class NodeType, class ListNodeType >
inline ListNodeType *tListIter< NodeType, ListNodeType >::
NodePtr()
{
  return curnode;
}

#endif
