/**************************************************************************\
**
**  tGridList.h
**
**  Header file for tGridList derived objects.
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
**  Modifications:
**   - added "MoveToActiveBack()" function, 12/97 GT
**
**  $Id: tMeshList.h,v 1.1 1998-01-14 20:19:32 gtucker Exp $
\**************************************************************************/
#ifndef TGRIDLIST_H
#define TGRIDLIST_H

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
   void insertAtFront( const NodeType & );
   int removeFromFront( NodeType & );
   void Flush();
   
  protected:
   int nActiveNodes;
   tListNode< NodeType > * lastactive;
};

#endif
