/**************************************************************************\
**
**  tGridList.h
**
**  Header file for derived classes tGridList and tGridListIter.
**
**  A tGridList is derived from the generic linked list class tList.
**  It is used in CHILD to store lists of grid elements (nodes and edges),
**  and differs from a generic list in being divided into two parts:
**  (1) an "active" part, representing elements that are not part of the
**  mesh boundary and are therefore subject to active processes (whatever
**  those may be; in CHILD the processes are runoff, erosion, and
**  sedimentation); and (2) a "boundary" part, containing elements along
**  the mesh boundary.
**
**  A tGridListIter is an iterator for a tGridList. It has the same services
**  as a tListIterator (getting the current node, moving to the first, last,
**  next, or previous node on the list, etc). It also will move to the
**  last "active" (non-boundary) node on a grid list, or to the first
**  boundary node on the list. It adds special functions FirstP, NextP
**  that are identical to the tListIter functions First and Next except
**  that they return a pointer to the data portion of the node (or zero if
**  the end of the list is reached, or the current node is null for some
**  other reason).
**
**  Modifications:
**   - added "MoveToActiveBack()" function, 12/97 GT
**
**  $Id: tMeshList.h,v 1.5 1999-01-05 22:46:27 stlancas Exp $
\**************************************************************************/

#ifndef TGRIDLIST_H
#define TGRIDLIST_H

#include "../Classes.h"
#include "../tList/tList.h"

/** class tGridList ********************************************************/
template< class NodeType >
class tGridList : public tList< NodeType >
{
   friend class tListIter< NodeType  >;
   friend class tGridListIter< NodeType  >;
  public:
   tGridList();
   tGridList( const tGridList< NodeType > * );
   ~tGridList();
   const tGridList< NodeType >
       &operator=( const tGridList< NodeType > & );
   int operator==( const tGridList< NodeType > & ) const;
   int operator!=( const tGridList< NodeType > & ) const;
   int getActiveSize() const;
   tListNode< NodeType  > * getLastActive() const;
   int isActiveEmpty() const;
   int isBoundEmpty() const;
   void insertAtBoundFront( const NodeType & );
   int removeFromBoundFront( NodeType & );
   void insertAtActiveBack( const NodeType & );
   int removeFromActiveBack( NodeType & );
   void setNActiveNodes( int );
   int removeNext( NodeType &value, tListNode< NodeType > * );
   int removePrev( NodeType &value, tListNode< NodeType > * );
   void moveToBack( tListNode< NodeType > * );
   void moveToFront( tListNode< NodeType > * );
   void moveToActiveBack( tListNode< NodeType > * );
   void moveToBoundFront( tListNode< NodeType > * );
   void moveToBack( NodeType * );
   void insertAtFront( const NodeType & );
   int removeFromFront( NodeType & );
   int InActiveList( tListNode< NodeType > * );
   void Flush();
   
  protected:
   int nActiveNodes;
   tListNode< NodeType > * lastactive;
};


/** class tGridListIter *****************************************************/
template< class NodeType >
class tGridListIter
                : public tListIter< NodeType >
{
  public:
   tGridListIter();
   tGridListIter( tGridList< NodeType > & );
   tGridListIter( tGridList< NodeType > * );
   ~tGridListIter();
   int LastActive();
   int FirstBoundary();
   int IsActive();
   NodeType * LastActiveP();
   NodeType * FirstBoundaryP();
//   NodeType * FirstP();
//   NodeType * NextP();
  //private:
   //tGridList< NodeType > *gridlistPtr; 
};


#endif
