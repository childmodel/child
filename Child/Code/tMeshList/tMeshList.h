//-*-c++-*-

/**************************************************************************/
/**
**  @file tMeshList.h
**  @brief Header file for derived classes tMeshList and tMeshListIter.
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
**  $Id: tMeshList.h,v 1.30 2004-06-16 13:37:38 childcvs Exp $
*/
/**************************************************************************/

#ifndef TMESHLIST_H
#define TMESHLIST_H

#include "../Classes.h"
#include "../tList/tList.h"
#include <assert.h>

/**************************************************************************/
/**
** @class tMeshList
**
** Class tMeshList implements a linked list that is divided into two
** parts, an "active" (front) and "inactive" (back) part. It is derived
** from tList.
**
*/
/**************************************************************************/
template< class NodeType, class ListNodeType >
class tMeshList : public tList< NodeType, ListNodeType >
{
   friend class tListIter< NodeType, ListNodeType  >;
   friend class tMeshListIter< NodeType, ListNodeType  >;

  public:
   tMeshList();
   tMeshList( const tMeshList< NodeType, ListNodeType > * );
   tMeshList(const tMeshList< NodeType, ListNodeType > &);
   ~tMeshList();
   const tMeshList< NodeType, ListNodeType >
       &operator=( const tMeshList< NodeType, ListNodeType > & );
   bool operator==( const tMeshList< NodeType, ListNodeType > & ) const;
   bool operator!=( const tMeshList< NodeType, ListNodeType > & ) const;
   inline int getActiveSize() const;
   inline ListNodeType * getLastActive();
   int isActiveEmpty() const;
   int isBoundEmpty() const;
   void insertAtBoundFront( const NodeType & );
   int removeFromBoundFront( NodeType & );
   void insertAtActiveBack( const NodeType & );
   int removeFromActiveBack( NodeType & );
   inline void setNActiveNodes( int );
   int removeNext( NodeType &value, ListNodeType * );
   int removePrev( NodeType &value, ListNodeType * );
   void moveToBack( ListNodeType * );
   void moveToFront( ListNodeType * );
   void moveToActiveBack( ListNodeType * );
   void moveToBoundFront( ListNodeType * );
   void moveToBack( NodeType const * );
   void insertAtFront( const NodeType & );
   int removeFromFront( NodeType & );
   inline void moveToBefore( ListNodeType*, ListNodeType* );
   inline void moveToAfter( ListNodeType*, ListNodeType* );
   int InActiveList( ListNodeType const * );
   void Flush();
   int CheckConsistency( const char * );

  protected:
   int nActiveNodes;                    // # of active nodes on list
   ListNodeType * lastactive;  // ptr to last active node
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
template< class NodeType, class ListNodeType >
tMeshList< NodeType, ListNodeType >::
tMeshList() :
  nActiveNodes(0),
  lastactive(0)
{
  if (0) //DEBUG
    std::cout << "                  from tMeshList()" << this->first << std::endl;
}

template< class NodeType, class ListNodeType >
tMeshList< NodeType, ListNodeType >::
tMeshList( const tMeshList< NodeType, ListNodeType > *original ) :
  tList< NodeType, ListNodeType >( original ),
  nActiveNodes(original->nActiveNodes),
  lastactive(0)
{
  if (0) //DEBUG
    std::cout << "                  from tMeshList( original )" << this->first
	 << std::endl;
  if ( nActiveNodes > 0 )
    lastactive = tList< NodeType, ListNodeType >::getIthListNodeNC( nActiveNodes - 1 );
}

template< class NodeType, class ListNodeType >
tMeshList< NodeType, ListNodeType >::
tMeshList( const tMeshList< NodeType, ListNodeType > &original ) :
  tList< NodeType, ListNodeType >( original ),
  nActiveNodes(original.nActiveNodes),
  lastactive(0)
{
  if (0) //DEBUG
    std::cout << "                  from tMeshList( original )" << this->first
	 << std::endl;
  if ( nActiveNodes > 0 )
    lastactive = tList< NodeType, ListNodeType >::getIthListNodeNC( nActiveNodes - 1 );
}

template< class NodeType, class ListNodeType >
tMeshList< NodeType, ListNodeType >::
~tMeshList()
{
  if (0) //DEBUG
    std::cout << "                  from ~tMeshList()" << this->first
	 << std::endl;
  lastactive = 0;
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
template< class NodeType, class ListNodeType >
const tMeshList< NodeType, ListNodeType > &tMeshList< NodeType, ListNodeType >::
operator=( const tMeshList< NodeType, ListNodeType > &right )
{
   if( this != &right )
   {
      tList< NodeType, ListNodeType >::operator=( right );
      nActiveNodes = right.nActiveNodes;
      lastactive = 0;
      if ( nActiveNodes > 0 )
	lastactive = tList< NodeType, ListNodeType >::getIthListNode( nActiveNodes - 1 );
   }
   return *this;
}

//overloaded equality operator:
template< class NodeType, class ListNodeType >
bool tMeshList< NodeType, ListNodeType >::
operator==( const tMeshList< NodeType, ListNodeType > &right ) const
{
   if( tList< NodeType, ListNodeType >::operator!=( right ) ) return false;
   if( nActiveNodes != right.nActiveNodes ) return false;
   if( lastactive != right.lastactive ) return false;
   return true;
}

//overloaded inequality operator:
template< class NodeType, class ListNodeType >
bool tMeshList< NodeType, ListNodeType >::
operator!=( const tMeshList< NodeType, ListNodeType > &right ) const
{
   return ! operator==(right);
}

/**************************************************************************\
**
**  tMeshList "get" functions
**
\**************************************************************************/

template< class NodeType, class ListNodeType >
inline int tMeshList< NodeType, ListNodeType >::
getActiveSize() const {return nActiveNodes;}

template< class NodeType, class ListNodeType >
inline ListNodeType *
tMeshList< NodeType, ListNodeType >::
getLastActive() {return lastactive;}

template< class NodeType, class ListNodeType >
inline void tMeshList< NodeType, ListNodeType >::
setNActiveNodes( int val ) {nActiveNodes = ( val >= 0 ) ? val : 0;}

template< class NodeType, class ListNodeType >
inline int tMeshList< NodeType, ListNodeType >::
isActiveEmpty() const
{
  if( lastactive == 0 ) return 1;
  else return 0;
}

template< class NodeType, class ListNodeType >
inline int tMeshList< NodeType, ListNodeType >::
isBoundEmpty() const
{
   if( lastactive == this->last ) return 1;
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

template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
insertAtFront( const NodeType &value )
{
   tList< NodeType, ListNodeType >::insertAtFront( value );
   if( value.isNonBoundary() )
   {
     if( isActiveEmpty() ) lastactive = this->first;
     ++nActiveNodes;
   }
}

template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
insertAtBoundFront( const NodeType &value )
{
   // Case list empty or active part of list empty:
   if( this->isEmpty() || lastactive==0 )
   {
      tList< NodeType, ListNodeType >::insertAtFront( value );
      return;
   }
   // Usual case: list and active part of list NOT empty:
   insertAtNext( value, lastactive );
}


template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
removeFromBoundFront( NodeType &value )
{
   if( lastactive == 0 ) return removeFromFront( value );
   return removeNext( value, lastactive );
}


template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
insertAtActiveBack( const NodeType &value )
{
   // Case list empty or active part of list empty:
   if( lastactive==0 )
   {
      insertAtFront( value );
      return;
   }
   // Usual case: list and active part of list NOT empty:
   tList< NodeType, ListNodeType >::insertAtNext( value, lastactive );
   lastactive = lastactive->next;
   ++nActiveNodes;
}

template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
removeFromActiveBack( NodeType &value )
{
   if( lastactive == 0 ) return 0;
   if( this->first == lastactive ) return removeFromFront( value );
   return removeNext( value, lastactive->prev );
}

template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
removeFromFront( NodeType &value )
{
   if( !isActiveEmpty() )
   {
      --nActiveNodes;
      if( lastactive == this->first ) lastactive = 0;
   }
   return tList< NodeType, ListNodeType >::removeFromFront( value );
}

//delete next node
template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
removeNext( NodeType &value, ListNodeType * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr->next == 0 ) return 0;
   if( value.isNonBoundary() )
   {
      --nActiveNodes;
      if( ptr->next == lastactive )
          lastactive = ptr;
   }
   return tList< NodeType, ListNodeType >::removeNext( value, ptr );
}

//delete previous node
template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
removePrev( NodeType &value, ListNodeType * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr->prev == 0 ) return 0;
   if( value.isNonBoundary() )
   {
      --nActiveNodes;
      if( ptr->prev == lastactive )
          lastactive = lastactive->prev;
   }
   return tList< NodeType, ListNodeType >::removePrev( value, ptr );
}


/**************************************************************************\
**
**  tMeshList::moveToBack ( ListNodeType * )
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
template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
moveToBack( ListNodeType * mvnode )
{
   if (0) //DEBUG
     std::cout << "moveToBack( ListNodeType )\n";

   assert( mvnode!=0 );
   if( mvnode != this->last )
   {
      //if( InActiveList( mvnode ) ) nActiveNodes--;
      if( mvnode->getDataPtr()->isNonBoundary() )
          --nActiveNodes;
      if( mvnode == lastactive )
      {
         if( mvnode != this->first )
         {
            lastactive = mvnode->prev;
         }
         else lastactive = 0;
      }
      tList< NodeType, ListNodeType >::moveToBack( mvnode );
   }
}


/**************************************************************************\
**
**  tMeshList::moveToBack ( NodeType * )
**
**  Finds the ListNode whose data are identical to mvnodedata and calls
**  moveToBack( ListNodeType ) to move it to the back of the list.
**
**  Parameters: mvnodedata -- ptr to data in node to be moved
**  Assumes: mvnodedata valid and contained in the list
**
\**************************************************************************/
template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
moveToBack( NodeType const * mvnodedata )
{
   ListNodeType * mvnode = getListNode( mvnodedata );
   assert( mvnode!=0 );  // failure: null or not on list
   moveToBack( mvnode );
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
template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
moveToFront( ListNodeType * mvnode )
{
   if( mvnode != this->first )
   {
      if( mvnode == lastactive )
      {
         lastactive = mvnode->prev;
      }
      tList< NodeType, ListNodeType >::moveToFront( mvnode );
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
template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
moveToActiveBack( ListNodeType * mvnode )
{
   if( !lastactive )
   {
      moveToFront( mvnode );
      lastactive = mvnode;
      if( ! mvnode->getDataPtr()->isNonBoundary() )
	nActiveNodes = 1;
      return;
   }

   if( mvnode != lastactive )
   {
      // if node was not in active part of list, increment nActiveNodes:
      if( ! mvnode->getDataPtr()->isNonBoundary() )
          ++nActiveNodes;
      // Detach mvnode from its position on the list:
      if( mvnode == this->first ) {
	assert( this->first->next );
	this->first =this-> first->next;
      }
      if( mvnode == this->last ) {
	assert( this->last->prev );
	this->last = this->last->prev;
      }
      if( mvnode->prev ) mvnode->prev->next = mvnode->next;
      if( mvnode->next ) mvnode->next->prev = mvnode->prev;
      // Insert it at the end of the active part of the list:
      mvnode->next = lastactive->next;
      mvnode->prev = lastactive;
      if( lastactive->next ) lastactive->next->prev = mvnode;
      lastactive->next = mvnode;
      // set lastactive and, if necessary last
      if( lastactive == this->last ) this->last = mvnode;
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
template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
moveToBoundFront( ListNodeType * mvnode )
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
      if( mvnode->getDataPtr()->isNonBoundary() )
          --nActiveNodes;
      if( mvnode == lastactive ) 
	{
	  assert( lastactive->prev );
	  lastactive = lastactive->prev;
	  return;
	}

      // Detach mvnode from its position on the list:
      if( mvnode == this->first ) {
         assert( this->first->next );
	 this->first = this->first->next;
      }
      if( mvnode == this->last ) {
         assert( this->last->prev );
	 this->last = this->last->prev;
      }
      if( mvnode->prev ) mvnode->prev->next = mvnode->next;
      if( mvnode->next ) mvnode->next->prev = mvnode->prev;
      // Insert it after the end of the active part of the list:
      mvnode->next = lastactive->next;
      mvnode->prev = lastactive;
      if( lastactive->next ) lastactive->next->prev = mvnode;
      lastactive->next = mvnode;
      if( lastactive == this->last ) this->last = mvnode;
   }
}


template< class NodeType, class ListNodeType >
inline void tMeshList< NodeType, ListNodeType >::
moveToBefore( ListNodeType* mvnode,
              ListNodeType* plcnode )
{
  tList< NodeType, ListNodeType >::moveToBefore( mvnode, plcnode );
  if( plcnode == lastactive->next
      && mvnode->getDataPtr()->isNonBoundary() )
    lastactive = mvnode;
}

template< class NodeType, class ListNodeType >
inline void tMeshList< NodeType, ListNodeType >::
moveToAfter( ListNodeType* mvnode,
	     ListNodeType* plcnode )
{
  tList< NodeType, ListNodeType >::moveToAfter( mvnode, plcnode );
  if( plcnode == lastactive
      && mvnode->getDataPtr()->isNonBoundary() )
    lastactive = mvnode;
}

/**************************************************************************\
**
**  tMeshList::Flush
**
**  Also reinitializes lastactive and nActiveNodes
**
\**************************************************************************/
template< class NodeType, class ListNodeType >
void tMeshList< NodeType, ListNodeType >::
Flush()
{
   tList< NodeType, ListNodeType >::Flush();
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
template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
InActiveList( ListNodeType const * theNode )
{
   ListNodeType * listnode = this->first;

   if( nActiveNodes==0 ) return 0;
   while( listnode!=lastactive && listnode!=theNode )
       listnode = listnode->next;
   if( listnode==theNode ) return 1;
   else return 0;

}

/**************************************************************************\
**
**  tMeshList::CheckConsistency
**
**  Internal consistency check
**
\**************************************************************************/
template< class NodeType, class ListNodeType >
int tMeshList< NodeType, ListNodeType >::
CheckConsistency( const char *ListName ){
  NodeType *cl;
  tMeshListIter<NodeType, ListNodeType> Iter( this );
  int nactive = 0;
  for( cl=Iter.FirstP(); ; cl=Iter.NextP() ){
    if (!Iter.IsActive()){
      std::cerr << "Element #" << cl->getID()
		<< " is not active but within the active part of the "
		<< ListName << " list.\n";
      return 1;
    }
    ++nactive;
    if (Iter.NodePtr() == this->getLastActive()) break;
    if (Iter.AtEnd()) {
      assert(0); /*NOTREACHED*/
      abort();
    }
  }
  if (nactive != this->getActiveSize()){
    std::cerr << "The " << ListName << " list contains " << nactive
	      << " elements but 'getActiveSize()' gives "
	      << this->getActiveSize() << ".\n";
    return 1;
  }
  for( cl=Iter.FirstBoundaryP(); !(Iter.AtEnd()); cl=Iter.NextP() ){
    if (Iter.IsActive()){
      std::cerr << "Element #" << cl->getID()
		<< " is active but within the boundary part of "
		<< ListName << " list.\n";
      return 1;
    }
  }
  return 0;
}

/**************************************************************************/
/**
** @class tMeshListIter
**
** Helper class for tMeshList, derived from tListIter ("iterators" that
** walk up and down a tList, fetching items -- see tList.h/.cpp).
** In addition to tListIter capabilities, tMeshListIter adds methods to
** move to and/or fetch the last active or first boundary (inactive)
** items, and to indicate whether it is on currently on the active portion
** of the list.
**
*/
/**************************************************************************/
template< class NodeType, class ListNodeType >
class tMeshListIter
                : public tListIter< NodeType, ListNodeType >
{
  public:
   tMeshListIter();
   tMeshListIter( tMeshList< NodeType, ListNodeType > & );
   tMeshListIter( tMeshList< NodeType, ListNodeType > * );
   ~tMeshListIter();
   int LastActive();
   int FirstBoundary();
   inline bool IsActive() const;
   NodeType * LastActiveP();
   NodeType * FirstBoundaryP();
   NodeType * GetNodePtrByPermID( int pid );
};

/**************************************************************************\
**     FUNCTIONS FOR DERIVED CLASS tMeshListIter
\**************************************************************************/

template< class NodeType, class ListNodeType >
tMeshListIter< NodeType, ListNodeType >::
tMeshListIter()
{
  if (0) //DEBUG
    std::cout << "    from tMeshListIter()" << std::endl;
}

template< class NodeType, class ListNodeType >
tMeshListIter< NodeType, ListNodeType >::
tMeshListIter( tMeshList< NodeType, ListNodeType > &list )
        : tListIter< NodeType, ListNodeType >( list )
{
   this->curnode = this->listPtr->first;
}

template< class NodeType, class ListNodeType >
tMeshListIter< NodeType, ListNodeType >::
tMeshListIter( tMeshList< NodeType, ListNodeType > *ptr )
        : tListIter< NodeType, ListNodeType >( ptr )
{
   assert( ptr != 0 );
   this->curnode = this->listPtr->first;
   assert( this->curnode != 0 );
}

template< class NodeType, class ListNodeType >
tMeshListIter< NodeType, ListNodeType >::
~tMeshListIter()
{}


/**************************************************************************\
**
**  tMeshListIter::LastActive
**
**  Moves the iterator to the last active node.
**
\**************************************************************************/
template< class NodeType, class ListNodeType >
int tMeshListIter< NodeType, ListNodeType >::
LastActive()
{
   tMeshList< NodeType, ListNodeType > *meshlistPtr =
     static_cast< tMeshList< NodeType, ListNodeType > * >(this->listPtr);
   assert( meshlistPtr != 0 );
   this->curnode = meshlistPtr->lastactive;
   this->counter = meshlistPtr->nActiveNodes-1;
   if( this->curnode != 0 ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tMeshListIter::FirstBoundary
**
**  Moves the iterator to the first boundary node.
**
\**************************************************************************/
template< class NodeType, class ListNodeType >
int tMeshListIter< NodeType, ListNodeType >::
FirstBoundary()
{
   tMeshList< NodeType, ListNodeType > *meshlistPtr =
     static_cast< tMeshList< NodeType, ListNodeType > * >(this->listPtr);
   assert( meshlistPtr != 0 );
   if( meshlistPtr->isActiveEmpty() ) this->curnode = this->listPtr->first;
   else if( meshlistPtr->isBoundEmpty() ) this->curnode = 0;
   else this->curnode = meshlistPtr->lastactive->next;
   this->counter = meshlistPtr->nActiveNodes;
   if( this->curnode != 0 ) return 1;
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
template< class NodeType, class ListNodeType >
NodeType* tMeshListIter< NodeType, ListNodeType >::
FirstBoundaryP()
{
  if ( FirstBoundary() ) return this->curnode->getDataPtrNC();
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
template< class NodeType, class ListNodeType >
NodeType *tMeshListIter< NodeType, ListNodeType >::
LastActiveP()
{
   if ( LastActive() ) return this->curnode->getDataPtrNC();
   else return 0;
}


/**************************************************************************\
**
**  tMeshListIter::IsActive
**
**  Indicates whether the current item is on the active portion of the
**  list, returning true if so, false if not. Assumes NodeType has a member
**  function isNonBoundary().
**
\**************************************************************************/
template< class NodeType, class ListNodeType >
inline bool tMeshListIter< NodeType, ListNodeType >::
IsActive() const
{
   if( this->curnode!=0 )
   {
      assert( this->curnode->getDataPtr()!=0 );
      return
	this->curnode->getDataRef().isNonBoundary();
   }
   return false;
}


#endif
