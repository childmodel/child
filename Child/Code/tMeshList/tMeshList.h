//-*-c++-*- 

/**************************************************************************\
**
**  tMeshList.h
**
**  Header file for derived classes tMeshList and tMeshListIter.
**  (formerly tGridList)
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
**
**  $Id: tMeshList.h,v 1.11 2002-09-04 16:39:04 arnaud Exp $
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
   tMeshList(const tMeshList&);
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
tMeshList()
  :
  nActiveNodes(0),
  lastactive(0)
{
     //cout << "                  from tMeshList()" << first << endl;
}

template< class NodeType >                     //tMeshtList
tMeshList< NodeType >::
tMeshList( const tMeshList< NodeType > *original )
  : tList< NodeType >( original ),
  nActiveNodes(original->nActiveNodes),
  lastactive(original->lastactive)
{
     //cout << "                  from tMeshList( original )" << first << endl;
}

template< class NodeType >                     //tMeshtList
tMeshList< NodeType >::
~tMeshList()
{
     //cout << "                  from ~tMeshList()" << first << endl;
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
      lastactive = first;
      int i;
      for( i=1; i<nActiveNodes; i++ ) lastactive = lastactive->next;
   }
     //cout << "tMeshList assigned" << first << endl;
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
int tMeshList< NodeType >::
getActiveSize() const {return nActiveNodes;}

template< class NodeType >                      //tMeshList
tListNode< NodeType > *
tMeshList< NodeType >::
getLastActive() const {return lastactive;}

template< class NodeType >                     //tMeshtList
void tMeshList< NodeType >::
setNActiveNodes( int val ) {nActiveNodes = ( val >= 0 ) ? val : 0;}

template< class NodeType >                     //tMeshtList
int tMeshList< NodeType >::
isActiveEmpty() const
{
     //cout << "checking if tMeshList empty of active nodes" << endl;
   if( lastactive == 0 )
   {
        //cout << "tMeshList is empty of active nodes" << endl;
      return 1;
   }
   else
   {
        //cout << "tMeshList is not empty of active nodes" << endl;
      return 0;
   }
}

template< class NodeType >                     //tMeshtList
int tMeshList< NodeType >::
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
   if( isActiveEmpty() ) lastactive = first;
   nActiveNodes++;
}

template< class NodeType >                     //tMeshtList
void tMeshList< NodeType >::
insertAtBoundFront( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node at boundary front of tMeshList" << endl;
   assert( newPtr>0 );
   assert( this != 0 );
   
   if( isEmpty() )  // Case list empty
   {
      first = last = newPtr;
   }
   else if( lastactive==0 ) // Case active part of list empty
   {
        //cout << "Case no active nodes" << endl;
      newPtr->next = first;
      first = newPtr;
   }
   else  // Usual case: list and active part of list NOT empty
   {
      newPtr->next = lastactive->next;
      lastactive->next = newPtr;
      if( lastactive==last ) last = newPtr; // Case: new is last (only) bdy
   }
   
     //nNodes++;  // why commented out? why not increment nnodes & nactive?
     //cout << "added node to front of boundary tMeshList" << endl;
    
}


template< class NodeType >                     //tMeshtList
int tMeshList< NodeType >::
removeFromBoundFront( NodeType &value )
{
   assert( &value != 0 );
   if( isEmpty() ) return 0;
   else if( last == lastactive ) return 0;
   else
   {
      tListNode< NodeType > * temp = lastactive->next;
      if( first == last ) first = last = 0;
      else lastactive->next = lastactive->next->next;
      value = temp->data;
      delete temp;
      nNodes--;
      return 1;
   }
}
   

template< class NodeType >                     //tMeshtList
void tMeshList< NodeType >::
insertAtActiveBack( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node at active back of tMeshList" << endl;
   assert( this != 0 );
   //cout << " isActiveEmpty() = " << isActiveEmpty() << endl;
   if( isEmpty() )
   {
      first = last = lastactive = newPtr;
   }
   else if( !( isEmpty() ) && isActiveEmpty() && !( isBoundEmpty() ) )
   {
      lastactive = newPtr;
      lastactive->next = first;
      first = lastactive;
      //cout << "first = lastactive: " << first->getDataRef().getID()
      //   << " " << lastactive->getDataRef().getID() << endl;
      //newPtr->next = first;
      //first = lastactive = newPtr;
   }
   else if( !( isEmpty() ) && isBoundEmpty() )
   {
      newPtr->next = lastactive->next;
      lastactive->next = newPtr;
      lastactive = newPtr;
      last = lastactive;
   }
   else
   {
      newPtr->next = lastactive->next;
      lastactive->next = newPtr;
      lastactive = newPtr;
   }
   if( isBoundEmpty() ) last = lastactive;
     //nNodes++;
   nActiveNodes++;
     //cout << "added node to back of active tMeshList" << endl;
}

template< class NodeType >                     //tMeshtList
int tMeshList< NodeType >::
removeFromActiveBack( NodeType &value )
{
   if( isEmpty() ) return 0;
   else
   {
      tListNode< NodeType > * temp = lastactive;
      if( first == lastactive ) lastactive = 0;
      if( first == last ) first = last = 0;
      else
      {
         tListNode< NodeType > * current = first;
         while( current->next != lastactive ) current = current->next;
         current->next = lastactive->next;
         lastactive->next = 0;
         lastactive = current;
      }
      value = temp->data;
      delete temp;
      nNodes--;
      nActiveNodes--;
      return 1;
   }
}

template< class NodeType >                         //tList
int tMeshList< NodeType >::
removeFromFront( NodeType &value )
{
   if( !( isActiveEmpty() ) )
   {
      nActiveNodes--;
      if( lastactive == first ) lastactive = 0;
   }
   return tList< NodeType >::removeFromFront( value );
}

//delete next node
template< class NodeType >                         //tList
int tMeshList< NodeType >::
removeNext( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr->next == 0 ) return 0;
   if( ptr == 0 ) return 0;
   if( ptr->next == lastactive ) return removeFromActiveBack( value );
   if( ptr == lastactive ) return removeFromBoundFront( value );
   if( tList< NodeType >::removeNext( value, ptr ) )
   {
      if( !( value.getBoundaryFlag() ) ) nActiveNodes--;
      return 1;
   }
   return 0;
}

//delete previous node
template< class NodeType >                         //tList
int tMeshList< NodeType >::
removePrev( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr == first && last->next == 0 ) return 0;
   if( lastactive->next == ptr ) return removeFromActiveBack( value );
   if( tList< NodeType >::removePrev( value, ptr ) )
   {
      if( !( value.getBoundaryFlag() ) ) nActiveNodes--;
      return 1;
   }
   return 0;
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
   //cout << "moveToBack( tListNode )\n";
   
   assert( mvnode>0 );
   tListNode< NodeType > * prev;
   if( mvnode != last )
   {
      if( InActiveList( mvnode ) )
         nActiveNodes--;
      if( mvnode == lastactive )
      {
         if( mvnode != first )
         {
            for( prev = first; prev->next != mvnode; prev = prev->next );
            lastactive = prev;
         }
         else
            lastactive = 0;
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
   tListNode< NodeType > *prev;
   if( mvnode != first )
   {
      if( mvnode == lastactive )
      {
         for( prev = first; prev->next != mvnode; prev = prev->next );
         lastactive = prev;
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
   tListNode< NodeType > * prev;
   if( mvnode != lastactive )
   {
      // Detach mvnode from its position on the list
      if( mvnode == first ) first = first->next;
      else
      {
         prev = first;
         while( prev->next != mvnode ) prev = prev->next;
         prev->next = mvnode->next;
      }

      // Insert it at the end of the active part of the list
      mvnode->next = lastactive->next;
      lastactive->next = mvnode;
      if( lastactive == last )
      {
         last = mvnode;
         // If it's a circular list, make sure to preserve circularity
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
   tListNode< NodeType > * prev;
   if( mvnode != lastactive->next )
   {
      // if node was in active part of list, decrement nActiveNodes
      if( InActiveList( mvnode ) ) --nActiveNodes;
      // Detach mvnode from its position on the list
      if( mvnode == first ) first = first->next;
      else
      {
         prev = first;
         while( prev->next != mvnode ) prev = prev->next;
         prev->next = mvnode->next;
      }

      // Insert it after the end of the active part of the list
      mvnode->next = lastactive->next;
      lastactive->next = mvnode;
      if( lastactive == last )
      {
         last = mvnode;
         // If it's a circular list, make sure to preserve circularity
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
     //cout << "    from tMeshListIter()" << endl;
}

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
tMeshListIter( tMeshList< NodeType > &list )
        : tListIter< NodeType >( list )
{
   assert( &list != 0 );
     //gridlistPtr = &list;
   curnode = /*grid*/listPtr->first;
     //if( listPtr->first != 0 ) assert( curnode != 0 );
     //cout << "    from tMeshListIter( list )" << endl;
}

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
tMeshListIter( tMeshList< NodeType > *ptr )
        : tListIter< NodeType >( ptr )
{
   assert( ptr != 0 );
     //gridlistPtr = &list;
   curnode = /*grid*/listPtr->first;
   assert( curnode != 0 );
     //cout << "    from tMeshListIter( ptr )" << endl;
}

template< class NodeType >   //tMeshListIter
tMeshListIter< NodeType >::
~tMeshListIter()
{
     //cout << "    from ~tMeshListIter()" << endl;
}


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
   tMeshList< NodeType > *meshlistPtr;

   meshlistPtr = static_cast< tMeshList< NodeType > * >(listPtr);
   assert( meshlistPtr != 0 );
   curnode = meshlistPtr->lastactive;
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
   tMeshList< NodeType > *meshlistPtr;
   meshlistPtr = static_cast< tMeshList< NodeType > * >(listPtr);
   assert( meshlistPtr != 0 );
   if( meshlistPtr->isActiveEmpty() ) curnode = listPtr->first;
   else if( meshlistPtr->isBoundEmpty() ) curnode = 0;
   else curnode = meshlistPtr->lastactive->next;
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
   tMeshList< NodeType > *meshlistPtr;
   meshlistPtr = static_cast< tMeshList< NodeType > * >(listPtr);
   assert( meshlistPtr != 0 );
   if( meshlistPtr->isActiveEmpty() ) curnode = listPtr->first;
   else if( meshlistPtr->isBoundEmpty() ) curnode = 0;
   else curnode = meshlistPtr->lastactive->next;
   if( curnode != 0 ) return curnode->getDataPtrNC();
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
   tMeshList< NodeType > *meshlistPtr;
   meshlistPtr = static_cast< tMeshList< NodeType > * >(listPtr);
   assert( meshlistPtr != 0 );
   curnode = meshlistPtr->lastactive;
   if( curnode != 0 ) return curnode->getDataPtrNC();
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
