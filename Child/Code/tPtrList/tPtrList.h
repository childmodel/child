/**************************************************************************\
**
**  tPtrList.h: Header file for tPtrList, tPtrListNode, and tPtrListIter
**              objects.
**
**  $Id: tPtrList.h,v 1.5 1998-02-03 00:48:30 stlancas Exp $
\**************************************************************************/


#ifndef TPTRLIST_H
#define TPTRLIST_H

#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include "../Classes.h"
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



/** class tPtrListIter ******************************************************/
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

