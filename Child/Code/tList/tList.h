/**************************************************************************\
**
**  tList.h: Header file for classes tList, tListNode, and tListIter.
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
**
**  $Id: tList.h,v 1.15 1999-04-05 15:33:15 gtucker Exp $
\**************************************************************************/

#ifndef TLIST_H
#define TLIST_H

#include "../Classes.h" // TODO: include only needed files


/**************************************************************************\
** class tListNode ********************************************************
**
** Class tListNode represents the items (or "nodes") on the list. Each
** tListNode object has two parts: the data (of type NodeType) and a
** pointer to the next item on the list. Capabilities include copy
** construction (from either another tListNode or a NodeType),
** returning a pointer or reference to either the data or the tListNode
** itself, and assignment and equality/inequality operations.
**
\**************************************************************************/
template< class NodeType >
class tListNode
{
    friend class tList< NodeType >;
    friend class tMeshList< NodeType >;
    friend class tListIter< NodeType >;
    friend class tMeshListIter< NodeType >;
public:
    tListNode();                                // default constructor
    tListNode( const tListNode< NodeType > & ); // copy constructor #1
    tListNode( const NodeType & );              // copy constructor #2
    const tListNode< NodeType >
    &operator=( const tListNode< NodeType > & );           // assignment
    int operator==( const tListNode< NodeType > & ) const; // equality
    int operator!=( const tListNode< NodeType > & ) const; // inequality
     /*set*/
    NodeType getDataNC();                     // returns copy of data item
    NodeType &getDataRefNC();                 // returns modifiable ref to data
    NodeType *getDataPtrNC();                 // returns modifiable ptr to data
    tListNode< NodeType > * getNextNC() const;// returns ptr to next list node
     /*get*/
    NodeType getData() const;                     // returns const copy of data
    const NodeType &getDataRef() const;           // returns const ref to data
    const NodeType *getDataPtr() const;           // returns const ptr to data
    const tListNode< NodeType > * getNext() const;// returns const ptr to next
   
protected:
    NodeType data;               // data item
    tListNode< NodeType > *next; // ptr to next node on list (=0 if end)
};


/**************************************************************************\
** class tList ************************************************************
**
** Class tList implements a linked list. The class includes pointers to
** the first, last, and current list nodes (see tListNode).
**
\**************************************************************************/
template< class NodeType >
class tList
{
    friend class tListIter< NodeType >;
    friend class tMeshListIter< NodeType >;
public:
    tList();                            // default constructor
    tList( const tList< NodeType > * ); // copy constructor
    ~tList();                           // destructor
    const tList< NodeType >
    &operator=( const tList< NodeType > & );           // assignment
    int operator==( const tList< NodeType > & ) const; // equality
    int operator!=( const tList< NodeType > & ) const; // inequality
    void insertAtFront( const NodeType & ); // puts copy of item at list front
    void insertAtBack( const NodeType & );  // puts copy of item at list back
    void insertAtNext( const NodeType &, tListNode< NodeType > * );
    void insertAtPrev( const NodeType &, tListNode< NodeType > * );
    int removeFromFront( NodeType & ); // removes 1st item, puts it in ref
    int removeFromBack( NodeType & );  // removes last item, puts it in ref
    int removeNext( NodeType &, tListNode< NodeType > * );
    int removePrev( NodeType &, tListNode< NodeType > * );
    void Flush();        // clears and reinitializes list
    int isEmpty() const; // returns 1 is list is empty, 0 otherwise
#ifndef NDEBUG
    void print() const;  // prints contents of list - DEBUG ONLY
#endif
    int getSize() const; // returns # of items on list
    tListNode< NodeType  > * getFirst() const; // returns ptr to 1st list node
    tListNode< NodeType  > * getLast() const;  // returns ptr to last list node
    NodeType * FirstP();  // returns ptr to 1st data item & sets current to 1st
    NodeType * NextP();   // moves to next node and returns ptr to data item
    void moveToBack( tListNode< NodeType > *  );  // move given node to back
    void moveToFront( tListNode< NodeType > *  ); // move given node to front
    void makeCircular();   // makes list circular (last points to first)
    //Xvoid setNNodes( int ); // ONLY sets nNodes -- danger Will Robinson!!
    const NodeType getIthData( int ) const;     // rtns copy of given item #
    const NodeType *getIthDataPtr( int ) const; // rtns ptr to given item #
    const NodeType &getIthDataRef( int ) const; // rtns ref to given item #
    NodeType getIthDataNC( int ) const;     // rtns modifiable copy of item #
    NodeType *getIthDataPtrNC( int ) const; // rtns modifiable ptr to item #
    NodeType &getIthDataRefNC( int ) const; // rtns modifiable ref to item #
    tListNode< NodeType > * getListNode( NodeType * ); // rtns ptr to node #
    
protected:
    int nNodes;                          // # of items on list
    tListNode< NodeType > * first;       // ptr to first node
    tListNode< NodeType > * last;        // ptr to last node
    tListNode< NodeType > * currentItem; // ptr to current item
    tListNode< NodeType > * getNewNode( const NodeType & ); // makes new node
};


/**************************************************************************\
** class tListIter *********************************************************
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
\**************************************************************************/
//TO DO: make Get, Where, GetP, refer to place in list rather than use getID()
template< class NodeType >
class tListIter
{
public:
    tListIter();                      // default constructor
    tListIter( tList< NodeType > & ); // constructor: reference to list
    tListIter( tList< NodeType > * ); // constructor: ptr to list
    ~tListIter();   // destructor
    int First();    // sets position to 1st list item (rtns 0 on failure)
    int Last();     // sets position to last "    "     "
    int Get( int ); // use only if NodeType has member getID()!!
    int Next();     // advances to next item (or 1st if current undefined)
    int Prev();     // moves to previous item (or last if current undef'd)
    int Where();    // use only if NodeType has member getID()!!
    int AtEnd();    // returns 1 if at end of the list
    NodeType &DatRef();  // returns ref to current data item
    NodeType *DatPtr();  // returns ptr to current data item
    tListNode< NodeType > *NodePtr();  // returns ptr to current list node
    void Reset( tList< NodeType > & ); // tells iterator to work on given list
    NodeType * FirstP();  // moves to 1st item and rtns ptr to it
    NodeType * LastP();   // moves to last  "   "  
    NodeType * NextP();   // moves to next  "   "
    NodeType * PrevP();   // moves to previous " "
    NodeType * GetP( int num ); //use only if NodeType has member getID()!!
    
protected:
    tListNode< NodeType > * curnode;  // ptr to current list node
    tList< NodeType > *listPtr;       // ptr to current list
    int counter;                      // current position on list (first=0)
};

#endif
