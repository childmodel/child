/**************************************************************************\
**
**  tMeshList.h
**
**  Header file for derived classes tMeshList and tMeshListIter.
**  (formerly tGridList)
**
**  A tMeshList is derived from the generic linked list class tList.
**  It is used in CHILD to store lists of grid elements (nodes and edges),
**  and differs from a generic list in being divided into two parts:
**  (1) an "active" part, representing elements that are not part of the
**  mesh boundary and are therefore subject to active processes (whatever
**  those may be; in CHILD the processes are runoff, erosion, and
**  sedimentation); and (2) a "boundary" part, containing elements along
**  the mesh boundary.
**
**  A tMeshListIter is an iterator for a tMeshList. It has the same services
**  as a tListIterator (getting the current node, moving to the first, last,
**  next, or previous node on the list, etc). It also will move to the
**  last "active" (non-boundary) node on a grid list, or to the first
**  boundary node on the list. It adds special functions FirstP, NextP
**  that are identical to the tListIter functions First and Next except
**  that they return a pointer to the data portion of the node (or zero if
**  the end of the list is reached, or the current node is null for some
**  other reason).
**
**  See also tList, tArray, tMatrix, tMesh
**
**  Modifications:
**   - added "MoveToActiveBack()" function, 12/97 GT
**
**  $Id: tMeshList.h,v 1.8 2000-12-07 12:08:21 gtucker Exp $
\**************************************************************************/

#ifndef TMESHLIST_H
#define TMESHLIST_H

#include "../Classes.h"
#include "../tList/tList.h"


/**************************************************************************\
** class tMeshList ********************************************************
**
** Class tMeshList implements a linked list that is divided into two
** parts, an "active" (front) and "inactive" (back) part. It is derived
** from tList.
**
\**************************************************************************/
template< class NodeType >
class tMeshList : public tList< NodeType >
{
   friend class tListIter< NodeType  >;
   friend class tMeshListIter< NodeType  >;
  public:
   tMeshList();
   tMeshList( const tMeshList< NodeType > * );
   ~tMeshList();
   const tMeshList< NodeType >
       &operator=( const tMeshList< NodeType > & );
   int operator==( const tMeshList< NodeType > & ) const;
   int operator!=( const tMeshList< NodeType > & ) const;
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
   int nActiveNodes;                    // # of active nodes on list
   tListNode< NodeType > * lastactive;  // ptr to last active node
};


/**************************************************************************\
** class tMeshListIter *****************************************************
**
** Helper class for tMeshList, derived from tListIter ("iterators" that
** walk up and down a tList, fetching items -- see tList.h/.cpp). 
** In addition to tListIter capabilities, tMeshListIter adds methods to
** move to and/or fetch the last active or first boundary (inactive)
** items, and to indicate whether it is on currently on the active portion
** of the list.
**
\**************************************************************************/
template< class NodeType >
class tMeshListIter
                : public tListIter< NodeType >
{
  public:
   tMeshListIter();
   tMeshListIter( tMeshList< NodeType > & );
   tMeshListIter( tMeshList< NodeType > * );
   ~tMeshListIter();
   int LastActive();
   int FirstBoundary();
   int IsActive();
   NodeType * LastActiveP();
   NodeType * FirstBoundaryP();
//   NodeType * FirstP();
//   NodeType * NextP();
  //private:
   //tMeshList< NodeType > *meshlistPtr; 
};

/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#ifdef __GNUC__
#include "tMeshList.cpp"
#endif

#endif
