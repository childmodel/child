//-*-c++-*- 

/**************************************************************************\
**
**  tMeshList.h
**
**  Header file for derived classes tMeshList and tMeshListIter.
**
**  A tMeshList is derived from the generic linked list class tList.
**  It is used in CHILD to store lists of mesh elements (nodes and edges),
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
**  last "active" (non-boundary) node on a mesh list, or to the first
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
**   - 09-2002 AD: Merge some of Stephen's bidirectional list patches
**
**  $Id: tMeshList.h,v 1.14 2002-09-06 16:42:52 arnaud Exp $
\**************************************************************************/

#ifndef TMESHLIST_H
#define TMESHLIST_H

#include "../Classes.h"
#include "../tList/tList.h"
#include "../tAssert.h"

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

   tMeshList(const tMeshList< NodeType > &);
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
**  FUNCTIONS FOR CLASS tMeshList
\**************************************************************************/

/**************************************************************************\
**
**  Constructors
**
**  Default: sets values to zero (empty list)
**  Copy constructor: creates a copy of _original_
**
\**************************************************************************/
template< class NodeType >                     //tMeshtList
tMeshList< NodeType >::
tMeshList() :
  nActiveNodes(0),
  lastactive(0)
{
  if (0) //DEBUG
    cout << "                  from tMeshList()" << first << endl;
}

template< class NodeType >                     //tMeshtList
tMeshList< NodeType >::
tMeshList( const tMeshList< NodeType > *original ) :
  tList< NodeType >( original ),
  nActiveNodes(original->nActiveNodes),
  lastactive(0)
{
  if (0) //DEBUG
    cout << "                  from tMeshList( original )" << first << endl;
  if ( nActiveNodes > 0 )
    lastactive = tList< NodeType >::getIthListNode( nActiveNodes - 1 );
}

template< class NodeType >                     //tMeshtList
tMeshList< NodeType >::
~tMeshList()
{
  if (0) //DEBUG
    cout << "                  from ~tMeshList()" << first << endl;
}


/**************************************************************************\
**
**  tMeshList overloaded operators
**
**  Assignment: creates a copy of right-hand list
**  Equality and inequality: adds test of # of active items and lastactive
**                           to basic tList operations
**
**  Modifications:
**    - assignment op: GT modified 12/7/99 to correct error in handling
**      lastactive pointer
**
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >                     //tMeshList
const tMeshList< NodeType > &tMeshList< NodeType >::
operator=( const tMeshList< NodeType > &right )
{
   if( this != &right )
   {
      tList< NodeType >::operator=( right );
      nActiveNodes = right.nActiveNodes;
      lastactive = 0;
      if ( nActiveNodes > 0 )
	lastactive = tList< NodeType >::getIthListNode( nActiveNodes - 1 );
   }
   return *this;
}

//overloaded equality operator:
template< class NodeType >                      //tMeshList
int tMeshList< NodeType >::
operator==( const tMeshList< NodeType > &right ) const
{
   if( tList< NodeType >::operator!=( right ) ) return 0;
   if( nActiveNodes != right.nActiveNodes ) return 0;
   if( lastactive != right.lastactive ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                      //tMeshList
int tMeshList< NodeType >::
operator!=( const tMeshList< NodeType > &right ) const
{
   return ! operator==(right);
}

/**************************************************************************\
**
**  tMeshList "get" functions
**
\**************************************************************************/

template< class NodeType >                      //tMeshList
inline int tMeshList< NodeType >::
getActiveSize() const {return nActiveNodes;}

template< class NodeType >                      //tMeshList
inline tListNode< NodeType > *
tMeshList< NodeType >::
getLastActive() const {return lastactive;}

template< class NodeType >                     //tMeshtList
inline void tMeshList< NodeType >::
setNActiveNodes( int val ) {nActiveNodes = ( val >= 0 ) ? val : 0;}

template< class NodeType >                     //tMeshtList
inline int tMeshList< NodeType >::
isActiveEmpty() const
{
  if( lastactive == 0 ) return 1;
  else return 0;
}

template< class NodeType >                     //tMeshtList
inline int tMeshList< NodeType >::
isBoundEmpty() const
{
   if( lastactive == last ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tMeshList insertion and removal functions
**
**  Adds and removes items to/from the list. Supplements tList
**  functionality by adding capability to add items to front of
**  "boundary" section or rear of "active" section. Updates
**  nActiveNodes as appropriate.
**
\**************************************************************************/

template< class NodeType >                         //tList
void tMeshList< NodeType >::
insertAtFront( const NodeType &value )
{
   tList< NodeType >::insertAtFront( value );
   if( value.getBoundaryFlag() == kNonBoundary )
   {
     if( isActiveEmpty() ) lastactive = first;
     ++nActiveNodes;
   }
}

template< class NodeType >                     //tMeshtList
void tMeshList< NodeType >::
insertAtBoundFront( const NodeType &value )
{
   // Case list empty or active part of list empty:
   if( isEmpty() || lastactive==0 )
   {
      tList< NodeType >::insertAtFront( value );
      return;
   }
   // Usual case: list and active part of list NOT empty:
   insertAtNext( value, lastactive );
}


template< class NodeType >                     //tMeshtList
int tMeshList< NodeType >::
removeFromBoundFront( NodeType &value )
{
   if( lastactive == 0 ) return removeFromFront( value );
   return removeNext( value, lastactive );
}
   

template< class NodeType >                     //tMeshtList
void tMeshList< NodeType >::
insertAtActiveBack( const NodeType &value )
{
   // Case list empty or active part of list empty:
   if( lastactive==0 )
   {
      insertAtFront( value );
      return;
   }
   // Usual case: list and active part of list NOT empty:
   tList< NodeType >::insertAtNext( value, lastactive );
   lastactive = lastactive->next;
   ++nActiveNodes;
}

template< class NodeType >                     //tMeshtList
int tMeshList< NodeType >::
removeFromActiveBack( NodeType &value )
{
   if( lastactive == 0 ) return 0;
   if( first == lastactive ) return removeFromFront( value );
   return removeNext( value, lastactive->prev );
}

template< class NodeType >                         //tList
int tMeshList< NodeType >::
removeFromFront( NodeType &value )
{
   if( !isActiveEmpty() )
   {
      --nActiveNodes;
      if( lastactive == first ) lastactive = 0;
   }
   return tList< NodeType >::removeFromFront( value );
}

//delete next node
template< class NodeType >                         //tList
int tMeshList< NodeType >::
removeNext( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr->next == 0 ) return 0;
   if( value.getBoundaryFlag() == kNonBoundary )
   {
      --nActiveNodes;
      if( ptr->next == lastactive )
          lastactive = ptr;
   }
   return tList< NodeType >::removeNext( value, ptr );
}

//delete previous node
template< class NodeType >                         //tList
int tMeshList< NodeType >::
removePrev( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr->prev == 0 ) return 0;
   if( value.getBoundaryFlag() == kNonBoundary )
   {
      --nActiveNodes;
      if( ptr->prev == lastactive )
          lastactive = lastactive->prev;
   }
   return tList< NodeType >::removePrev( value, ptr );
}


/**************************************************************************\
**
**  tMeshList::moveToBack ( tListNode * )
**
**  Moves mvnode to the back of the list (the boundary portion).
**  Handles case of moved node being the last active node, in which case
**  _lastactive_ needs to be updated.
**
**  Modifications:
**    - if moved node is active, nActiveNodes is now decremented 4/98 GT
**        (note: does not properly handle the case of list w/ only one node
**      that's active -- in this case, node is not moved (it's already last)
**      and nActiveNodes isn't updated. TODO) 
**
\**************************************************************************/
template< class NodeType >                         //tList
void tMeshList< NodeType >::
moveToBack( tListNode< NodeType > * mvnode ) 
{
   if (0) //DEBUG
     cout << "moveToBack( tListNode )\n";
   
   assert( mvnode>0 );
   if( mvnode != last )
   {
      //if( InActiveList( mvnode ) ) nActiveNodes--;
      if( mvnode->getDataPtr()->getBoundaryFlag() == kNonBoundary )
          --nActiveNodes;
      if( mvnode == lastactive )
      {
         if( mvnode != first )
         {
            lastactive = mvnode->prev;
         }
         else lastactive = 0;
      }
      tList< NodeType >::moveToBack( mvnode );
   }
}


/**************************************************************************\
**
**  tMeshList::moveToBack ( NodeType * )
**
**  Finds the ListNode whose data are identical to mvnodedata and calls
**  moveToBack( tListNode ) to move it to the back of the list.
**
**  Parameters: mvnodedata -- ptr to data in node to be moved
**  Assumes: mvnodedata valid and contained in the list
**
\**************************************************************************/
template< class NodeType >                         //tList
void tMeshList< NodeType >::
moveToBack( NodeType * mvnodedata ) 
{
   assert( getListNode( mvnodedata )!=0 );  // failure: null or not on list
   moveToBack( getListNode( mvnodedata ) );
}


/**************************************************************************\
**
**  tMeshList::moveToFront
**
**  Moves mvnode to the front of the list, taking care to handle the case
**  in which the node being moved is the last on the active section
**  (doesn't check whether node is active or inactive however, and thus
**  doesn't update nActiveNodes...TODO)
**
\**************************************************************************/
template< class NodeType >                         //tList
void tMeshList< NodeType >::
moveToFront( tListNode< NodeType > * mvnode ) 
{
   if( mvnode != first )
   {
      if( mvnode == lastactive )
      {
         lastactive = mvnode->prev;
      }
      tList< NodeType >::moveToFront( mvnode );
   }
}


/**************************************************************************\
**
**  tMeshList::moveToActiveBack
**
**  Moves mvnode to the back of the "active" portion of the list
**  (does not update nActiveNodes if the node happens to be inactive!)
**
\**************************************************************************/
template< class NodeType >                         //tList
void tMeshList< NodeType >::
moveToActiveBack( tListNode< NodeType > * mvnode ) 
{
   if( !lastactive )
   {
      moveToBack( mvnode );
      lastactive = mvnode;
      return;
   }
   
   if( mvnode != lastactive )
   {
      // if node was in active part of list, decrement nActiveNodes:
      if( mvnode->getDataPtr()->getBoundaryFlag() != kNonBoundary )
          ++nActiveNodes;
      // Detach mvnode from its position on the list:
      if( mvnode == first )
      {
         //if( first->prev != 0 )
         //    first->prev->next = first->next;
         //first->next->prev = first->prev;
         if( first->next ) first = first->next;
      }
      //else
      //{
      if( mvnode->prev ) mvnode->prev->next = mvnode->next;
      if( mvnode->next ) mvnode->next->prev = mvnode->prev;
         //}
      // Insert it at the end of the active part of the list:
      mvnode->next = lastactive->next;
      mvnode->prev = lastactive;
      if( lastactive->next ) lastactive->next->prev = mvnode;
      lastactive->next = mvnode;
      if( lastactive == last )
      {
         last = mvnode;
         // If it's a circular list, make sure to preserve circularity:
         if( last->next != 0 ) last->next = first;
      }
      lastactive = mvnode;
   }
}


/**************************************************************************\
**
**  tMeshList::moveToBoundFront
**
**  Moves mvnode to the front of the "boundary" portion of the list,
**  making sure to update nActiveNodes is the node was previously on
**  the active portion of the list.
**
\**************************************************************************/
template< class NodeType >                         //tList
void tMeshList< NodeType >::
moveToBoundFront( tListNode< NodeType > * mvnode ) 
{
   if( !lastactive )
   {
      moveToFront( mvnode );
      return;
   }
   
   if( mvnode != lastactive->next )
   {
      // if node was in active part of list, decrement nActiveNodes:
      //if( InActiveList( mvnode ) ) --nActiveNodes;
      if( mvnode->getDataPtr()->getBoundaryFlag() == kNonBoundary )
          --nActiveNodes;
      // Detach mvnode from its position on the list:
      if( mvnode == first )
      {
         //if( first->next ) first->next->prev = first->prev;
         if( first->next ) first = first->next;
      }
      if( mvnode == lastactive )
      {
         lastactive = lastactive->prev;
         return;
      }
      //else
      //{
      if( mvnode->prev ) mvnode->prev->next = mvnode->next;
      if( mvnode->next ) mvnode->next->prev = mvnode->prev;
         //}
      // Insert it after the end of the active part of the list:
      mvnode->next = lastactive->next;
      mvnode->prev = lastactive;
      if( lastactive->next ) lastactive->next->prev = mvnode;
      lastactive->next = mvnode;
      if( lastactive == last )
      {
         last = mvnode;
         // If it's a circular list, make sure to preserve circularity:
         if( last->next != 0 ) last->next = first;
      }
   }
}


/**************************************************************************\
**
**  tMeshList::Flush
**
**  Also reinitializes lastactive and nActiveNodes
**
\**************************************************************************/
template< class NodeType >                         //tList
void tMeshList< NodeType >::
Flush()
{
   tList< NodeType >::Flush();
   lastactive = 0;
   nActiveNodes = 0;
}


/**************************************************************************\
**
**  tMeshList::InActiveList
**
**  Reports whether a given list node is in the active portion of the list.
**
**  Parameters:  theNode -- list node to test
**  Returns:  1 if theNode is present in the active portion of the list,
**            0 otherwise.
**  Created:  4/29/98 GT
**
\**************************************************************************/
template< class NodeType >                         //tList
int tMeshList< NodeType >::
InActiveList( tListNode< NodeType > * theNode )
{
   tListNode< NodeType > * listnode = first;

   if( nActiveNodes==0 ) return 0;
   while( listnode!=lastactive && listnode!=theNode )
       listnode = listnode->next;
   if( listnode==theNode ) return 1;
   else return 0;
   
}

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

/**************************************************************************\
**     FUNCTIONS FOR DERIVED CLASS tMeshListIter
\**************************************************************************/

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
tMeshListIter()
{
  if (0) //DEBUG
    cout << "    from tMeshListIter()" << endl;
}

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
tMeshListIter( tMeshList< NodeType > &list )
        : tListIter< NodeType >( list )
{
   curnode = listPtr->first;
}

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
tMeshListIter( tMeshList< NodeType > *ptr )
        : tListIter< NodeType >( ptr )
{
   assert( ptr != 0 );
   curnode = listPtr->first;
   assert( curnode != 0 );
}

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
~tMeshListIter()
{}


/**************************************************************************\
**
**  tMeshListIter::LastActive
**
**  Moves the iterator to the last active node.
**
\**************************************************************************/
template< class NodeType >   //tMeshListIter
int tMeshListIter< NodeType >::
LastActive()
{
   tMeshList< NodeType > *meshlistPtr = 
     static_cast< tMeshList< NodeType > * >(listPtr);
   assert( meshlistPtr != 0 );
   curnode = meshlistPtr->lastactive;
   counter = meshlistPtr->nActiveNodes-1;
   if( curnode != 0 ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tMeshListIter::FirstBoundary
**
**  Moves the iterator to the first boundary node.
**
\**************************************************************************/
template< class NodeType >   //tMeshListIter
int tMeshListIter< NodeType >::
FirstBoundary()
{
   tMeshList< NodeType > *meshlistPtr = 
     static_cast< tMeshList< NodeType > * >(listPtr);
   assert( meshlistPtr != 0 );
   if( meshlistPtr->isActiveEmpty() ) curnode = listPtr->first;
   else if( meshlistPtr->isBoundEmpty() ) curnode = 0;
   else curnode = meshlistPtr->lastactive->next;
   counter = meshlistPtr->nActiveNodes;
   if( curnode != 0 ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tMeshListIter::FirstBoundaryP
**
**  Moves the iterator to the first boundary node and returns a pointer
**  to the data at that location.
**
\**************************************************************************/
template< class NodeType >   //tMeshListIter
NodeType* tMeshListIter< NodeType >::
FirstBoundaryP()
{
  if ( FirstBoundary() ) return curnode->getDataPtrNC();
  else return 0;
}


/**************************************************************************\
**
**  tMeshListIter::LastActiveP
**
**  Moves the iterator to the last active node and returns a pointer
**  to the data at that location.
**
\**************************************************************************/
template< class NodeType >   //tMeshListIter
NodeType *tMeshListIter< NodeType >::
LastActiveP()
{
   if ( LastActive() ) return curnode->getDataPtrNC();
   else return 0;
}


/**************************************************************************\
**
**  tMeshListIter::IsActive
**
**  Indicates whether the current item is on the active portion of the
**  list, returning 1 if so, 0 if not. Assumes NodeType has a member
**  function getBoundaryFlag.
**
\**************************************************************************/
template< class NodeType >   //tMeshListIter
int tMeshListIter< NodeType >::
IsActive()
{
   if( curnode!=0 )
   {
      assert( curnode->getDataPtr()!=0 );
      return 
	curnode->getDataRef().getBoundaryFlag() == kNonBoundary;
   }
   return 0;
}

#endif
