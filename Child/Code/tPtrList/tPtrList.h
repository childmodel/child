/**************************************************************************\
**
**  tPtrList.h: Header file for tPtrList and tPtrListNode objects.
**
**  $Id: tPtrList.h,v 1.1 1998-01-21 00:49:19 gtucker Exp $
\**************************************************************************/


#ifndef TPTRLIST_H
#define TPTRLIST_H

/** class tPtrListNode ******************************************************/
template< class NodeType >
class tPtrListNode
{
   friend class tPtrList< NodeType >;
   friend class tPtrListIter< NodeType >;
  public:
   tPtrListNode();
   tPtrListNode( const tPtrListNode< NodeType > & );
   tPtrListNode( NodeType * );
   ~tPtrListNode();
   const tPtrListNode< NodeType >
       &operator=( const tPtrListNode< NodeType > & );
   int operator==( const tPtrListNode< NodeType > & ) const;
   int operator!=( const tPtrListNode< NodeType > & ) const;
   NodeType *getPtrNC();
   const NodeType *getPtr() const;
   tPtrListNode< NodeType > *getNextNC();
   const tPtrListNode< NodeType > * getNext() const;
  private:
   NodeType * Ptr;
   tPtrListNode< NodeType > * next;
};


/** class tPtrList **********************************************************/
template< class NodeType >
class tPtrList
{
   friend class tPtrListIter< NodeType >;
  public:
   tPtrList();
   tPtrList( const tPtrList< NodeType > & );
   ~tPtrList();
   const tPtrList< NodeType > &operator=( const tPtrList< NodeType > & );
   void insertAtFront( NodeType * );
   void insertAtBack( NodeType * );
   void insertAtNext( NodeType *, tPtrListNode< NodeType > * );
   void insertAtPrev( NodeType *, tPtrListNode< NodeType > * );
   int removeFromFront( NodeType * );
   int removeFromBack( NodeType * );
   int removeNext( NodeType *, tPtrListNode< NodeType > * );
   int removePrev( NodeType *, tPtrListNode< NodeType > * );
   void Flush();
   int isEmpty() const;
   void print() const;
   /*void input( int, tList< NodeType > * );*/
   int getSize() const;
   tPtrListNode< NodeType > * getFirstNC();
   const tPtrListNode< NodeType > * getFirst() const;
   void moveToBack( tPtrListNode< NodeType > *  );
   void moveToFront( tPtrListNode< NodeType > *  );
   tPtrListNode< NodeType > * getLast() const;
   void makeCircular();
   const NodeType *getIthPtr( int ) const;
   NodeType *getIthPtrNC( int ) const;
  private:
   int nNodes;
   tPtrListNode< NodeType > * first;
   tPtrListNode< NodeType > * last;
   tPtrListNode< NodeType > * getNewNode( NodeType * );
};

#endif
