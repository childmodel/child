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
**  $Id: tPtrList.h,v 1.8 2000-01-13 19:19:54 gtucker Exp $
\**************************************************************************/

#ifndef TPTRLIST_H
#define TPTRLIST_H

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "../Classes.h" // TODO: include only needed stuff


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
    ~tPtrList();                              // destructor
    const tPtrList< NodeType > 
    &operator=( const tPtrList< NodeType > & );  // assignment
    void insertAtFront( NodeType * ); // puts ptr at list front
    void insertAtBack( NodeType * );  // puts ptr at list back
    void insertAtNext( NodeType *, tPtrListNode< NodeType > * );
    void insertAtPrev( NodeType *, tPtrListNode< NodeType > * );
    int removeFromFront( NodeType * );  // removes 1st item, puts in ptr
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
\**************************************************************************/
//TO DO: make Get, Where, GetP, refer to place in list rather than use getID()
template< class NodeType >
class tPtrListIter
{
  public:
   tPtrListIter();
   tPtrListIter( tPtrList< NodeType > & );
   ~tPtrListIter();
   int First();
   int Last();
   int Get( int );
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

#endif


