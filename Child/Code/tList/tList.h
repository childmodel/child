/**************************************************************************\
**
**  tList.h: Header file for classes tList, tListNode, and tListIter.
**
**  A tList is an object that implements a general linked list NodeType
**  objects, where NodeType can be any type (float, int, other objects,
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
**  of the given data type (float, int, class, etc) and a pointer to the
**  next node in the list.
**
**  A tListIter is an iterator for the linked list tList objects (and their
**  descendants). Its services include fetching data from the current entry
**  on the list, advancing to the next or previous item on the list, etc.
**
**  $Id: tList.h,v 1.4 1998-01-21 22:09:29 gtucker Exp $
\**************************************************************************/

#ifndef TLIST_H
#define TLIST_H

template<class NodeType> class tGridList<NodeType>;
template<class NodeType> class tGridListIter<NodeType>;


/** class tListNode  ********************************************************/
template< class NodeType >
class tListNode
{
   friend class tList< NodeType >;
   friend class tGridList< NodeType >;
   friend class tListIter< NodeType >;
   friend class tGridListIter< NodeType >;
  public:
   tListNode();
   tListNode( const tListNode< NodeType > & );
   tListNode( const NodeType & );
   const tListNode< NodeType >
       &operator=( const tListNode< NodeType > & );
   int operator==( const tListNode< NodeType > & ) const;
   int operator!=( const tListNode< NodeType > & ) const;
     /*set*/
   NodeType getDataNC();
   NodeType &getDataRefNC();
   NodeType *getDataPtrNC();
   tListNode< NodeType > * getNextNC() const;
     /*get*/
   NodeType getData() const;
   const NodeType &getDataRef() const;
   const NodeType *getDataPtr() const;
   const tListNode< NodeType > * getNext() const;
   
  protected:
   NodeType data;
   tListNode< NodeType > *next;
};


/** class tList ************************************************************/
template< class NodeType >
class tList
{
   friend class tListIter< NodeType >;
   friend class tGridListIter< NodeType >;
     /*friend class tListIter< NodeType, tNode >;
   friend class tListIter< NodeType, tLNode >;
   friend class tListIter< NodeType, tGrid >;
   friend class tListIter< NodeType, tGrid >;*/
  public:
   tList();
   tList( const tList< NodeType > * );
   ~tList();
   const tList< NodeType >
       &operator=( const tList< NodeType > & );
   int operator==( const tList< NodeType > & ) const;
   int operator!=( const tList< NodeType > & ) const;
   void insertAtFront( const NodeType & );
   void insertAtBack( const NodeType & );
   void insertAtNext( const NodeType &, tListNode< NodeType > * );
   void insertAtPrev( const NodeType &, tListNode< NodeType > * );
   int removeFromFront( NodeType & );
   int removeFromBack( NodeType & );
   int removeNext( NodeType &, tListNode< NodeType > * );
   int removePrev( NodeType &, tListNode< NodeType > * );
   void Flush();
   int isEmpty() const;
   void print() const;
   /*void input( int );*/
   int getSize() const;
   tListNode< NodeType  > * getFirst() const;
   tListNode< NodeType  > * getLast() const;
   void moveToBack( tListNode< NodeType > *  );
   void moveToFront( tListNode< NodeType > *  );
   void makeCircular();
   void setNNodes( int );
   const NodeType getIthData( int ) const;
   const NodeType *getIthDataPtr( int ) const;
   const NodeType &getIthDataRef( int ) const;
   NodeType getIthDataNC( int ) const;
   NodeType *getIthDataPtrNC( int ) const;
   NodeType &getIthDataRefNC( int ) const;
   
  protected:
   int nNodes;
   tListNode< NodeType > * first;
   tListNode< NodeType > * last;
   tListNode< NodeType > * getNewNode( const NodeType & );
};



/** class tListIter ********************************************************/
template< class NodeType >
class tListIter
{
  public:
   tListIter();
   tListIter( tList< NodeType > & );
   tListIter( tList< NodeType > * );
   ~tListIter();
   int First();
   int Last();
   int Get( int );
   int Next();
   int Prev();
   int Where();
   int AtEnd();
   NodeType &DatRef();
   NodeType *DatPtr();
   tListNode< NodeType > *NodePtr();
   void Reset( tList< NodeType > & );
   NodeType * FirstP();
   NodeType * LastP();
   NodeType * NextP();
   NodeType * PrevP();
   NodeType * GetP( int num );
   
  protected:
   tListNode< NodeType > * curnode;
   tList< NodeType > *listPtr;   
};

#endif
