/**************************************************************************\
**
**  tList.h: Header file for class tList< NodeType >
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
**  $Id: tList.h,v 1.1 1998-01-14 20:33:55 gtucker Exp $
\**************************************************************************/
#ifndef TLIST_H
#define TLIST_H

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

#endif
