/**************************************************************************\
**
**  tGridList.cpp
**
**  Functions for derived classes tGridList and tGridListIter. The classes
**  are declared in tGridList.h (q.v.).
**
**  Modifications:
**   - added "MoveToActiveBack()" function, 12/97 GT
**
**  $Id: tMeshList.cpp,v 1.7 1999-04-02 22:17:38 gtucker Exp $
\**************************************************************************/

#include <assert.h>
#include "tGridList.h"


/**************************************************************************\
**  FUNCTIONS FOR CLASS tGridList
\**************************************************************************/

/**************************************************************************\
**
**  Constructors
**
**  Default: sets values to zero (empty list)
**  Copy constructor: creates a copy of _original_
**
\**************************************************************************/
template< class NodeType >                     //tGridtList
tGridList< NodeType >::
tGridList()
{
   nActiveNodes = 0;
   lastactive = 0;
     //cout << "                  from tGridList()" << first << endl;
}

template< class NodeType >                     //tGridtList
tGridList< NodeType >::
tGridList( const tGridList< NodeType > *original )
        : tList< NodeType >( original )
{
   nActiveNodes = original->nActiveNodes;
   if( original->lastactive != 0 ) lastactive = original->lastactive;
   else lastactive = 0;
     //cout << "                  from tGridList( original )" << first << endl;
}

template< class NodeType >                     //tGridtList
tGridList< NodeType >::
~tGridList()
{
     //cout << "                  from ~tGridList()" << first << endl;
}


/**************************************************************************\
**
**  tGridList overloaded operators
**
**  Assignment: creates a copy of right-hand list
**  Equality and inequality: adds test of # of active items and lastactive
**                           to basic tList operations
**
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >                     //tGridList
const tGridList< NodeType > &tGridList< NodeType >::
operator=( const tGridList< NodeType > &right )
{
   if( this != &right )
   {
      tList< NodeType >::operator=( right );
      lastactive = right.lastactive;
      nActiveNodes = right.nActiveNodes;
   }
     //cout << "tGridList assigned" << first << endl;
   return *this;
}

//overloaded equality operator:
template< class NodeType >                      //tGridList
int tGridList< NodeType >::
operator==( const tGridList< NodeType > &right ) const
{
   if( tList< NodeType >::operator!=( right ) ) return 0;
   if( nActiveNodes != right.nActiveNodes ) return 0;
   if( lastactive != right.lastactive ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                      //tGridList
int tGridList< NodeType >::
operator!=( const tGridList< NodeType > &right ) const
{
   if( tList< NodeType >::operator!=( right ) ) return 1;
   if( nActiveNodes != right.nActiveNodes ) return 1;
   if( lastactive != right.lastactive ) return 1;
   return 0;
}

/**************************************************************************\
**
**  tGridList "get" functions
**
\**************************************************************************/

template< class NodeType >                      //tGridList
int tGridList< NodeType >::
getActiveSize() const {return nActiveNodes;}

template< class NodeType >                      //tGridList
tListNode< NodeType > *
tGridList< NodeType >::
getLastActive() const {return lastactive;}

template< class NodeType >                     //tGridtList
void tGridList< NodeType >::
setNActiveNodes( int val ) {nActiveNodes = ( val >= 0 ) ? val : 0;}

template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
isActiveEmpty() const
{
     //cout << "checking if tGridList empty of active nodes" << endl;
   if( lastactive == 0 )
   {
        //cout << "tGridList is empty of active nodes" << endl;
      return 1;
   }
   else
   {
        //cout << "tGridList is not empty of active nodes" << endl;
      return 0;
   }
}

template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
isBoundEmpty() const
{
   if( lastactive == last ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tGridList insertion and removal functions
**
**  Adds and removes items to/from the list. Supplements tList
**  functionality by adding capability to add items to front of
**  "boundary" section or rear of "active" section. Updates
**  nActiveNodes as appropriate.
**
\**************************************************************************/

template< class NodeType >                         //tList
void tGridList< NodeType >::
insertAtFront( const NodeType &value )
{
   tList< NodeType >::insertAtFront( value );
   if( isActiveEmpty() ) lastactive = first;
   nActiveNodes++;
}

template< class NodeType >                     //tGridtList
void tGridList< NodeType >::
insertAtBoundFront( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node at boundary front of tGridList" << endl;
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
     //cout << "added node to front of boundary tGridList" << endl;
    
}


template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
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
   

template< class NodeType >                     //tGridtList
void tGridList< NodeType >::
insertAtActiveBack( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
     //cout << "add new node at active back of tGridList" << endl;
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
     //cout << "added node to back of active tGridList" << endl;
}

template< class NodeType >                     //tGridtList
int tGridList< NodeType >::
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
int tGridList< NodeType >::
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
int tGridList< NodeType >::
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
int tGridList< NodeType >::
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
**  tGridList::moveToBack ( tListNode * )
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
void tGridList< NodeType >::
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
**  tGridList::moveToBack ( NodeType * )
**
**  Finds the ListNode whose data are identical to mvnodedata and calls
**  moveToBack( tListNode ) to move it to the back of the list.
**
**  Parameters: mvnodedata -- ptr to data in node to be moved
**  Assumes: mvnodedata valid and contained in the list
**
\**************************************************************************/
template< class NodeType >                         //tList
void tGridList< NodeType >::
moveToBack( NodeType * mvnodedata ) 
{
   assert( getListNode( mvnodedata )!=0 );  // failure: null or not on list
   moveToBack( getListNode( mvnodedata ) );
}


/**************************************************************************\
**
**  tGridList::moveToFront
**
**  Moves mvnode to the front of the list, taking care to handle the case
**  in which the node being moved is the last on the active section
**  (doesn't check whether node is active or inactive however, and thus
**  doesn't update nActiveNodes...TODO)
**
\**************************************************************************/
template< class NodeType >                         //tList
void tGridList< NodeType >::
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
**  tGridList::moveToActiveBack
**
**  Moves mvnode to the back of the "active" portion of the list
**  (does not update nActiveNodes if the node happens to be inactive!)
**
\**************************************************************************/
template< class NodeType >                         //tList
void tGridList< NodeType >::
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
**  tGridList::moveToBoundFront
**
**  Moves mvnode to the front of the "boundary" portion of the list,
**  making sure to update nActiveNodes is the node was previously on
**  the active portion of the list.
**
\**************************************************************************/
template< class NodeType >                         //tList
void tGridList< NodeType >::
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
**  tGridList::Flush
**
**  Also reinitializes lastactive and nActiveNodes
**
\**************************************************************************/
template< class NodeType >                         //tList
void tGridList< NodeType >::
Flush()
{
   tList< NodeType >::Flush();
   lastactive = 0;
   nActiveNodes = 0;
}


/**************************************************************************\
**
**  tGridList::InActiveList
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
int tGridList< NodeType >::
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
**     FUNCTIONS FOR DERIVED CLASS tGridListIter
\**************************************************************************/

template< class NodeType >   //tGridListIter
tGridListIter< NodeType >::
tGridListIter()
{
     //cout << "    from tGridListIter()" << endl;
}

template< class NodeType >   //tGridListIter
tGridListIter< NodeType >::
tGridListIter( tGridList< NodeType > &list )
        : tListIter< NodeType >( list )
{
   assert( &list != 0 );
     //gridlistPtr = &list;
   curnode = /*grid*/listPtr->first;
     //if( listPtr->first != 0 ) assert( curnode != 0 );
     //cout << "    from tGridListIter( list )" << endl;
}

template< class NodeType >   //tGridListIter
tGridListIter< NodeType >::
tGridListIter( tGridList< NodeType > *ptr )
        : tListIter< NodeType >( ptr )
{
   assert( ptr != 0 );
     //gridlistPtr = &list;
   curnode = /*grid*/listPtr->first;
   assert( curnode != 0 );
     //cout << "    from tGridListIter( ptr )" << endl;
}

template< class NodeType >   //tGridListIter
tGridListIter< NodeType >::
~tGridListIter()
{
     //cout << "    from ~tGridListIter()" << endl;
}


/**************************************************************************\
**
**  tGridListIter::LastActive
**
**  Moves the iterator to the last active node.
**
\**************************************************************************/
template< class NodeType >   //tGridListIter
int tGridListIter< NodeType >::
LastActive()
{
   tGridList< NodeType > *gridlistPtr;

   gridlistPtr = ( tGridList< NodeType > * ) listPtr;
   assert( gridlistPtr != 0 );
   curnode = gridlistPtr->lastactive;
   if( curnode != 0 ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tGridListIter::FirstBoundary
**
**  Moves the iterator to the first boundary node.
**
\**************************************************************************/
template< class NodeType >   //tGridListIter
int tGridListIter< NodeType >::
FirstBoundary()
{
   tGridList< NodeType > *gridlistPtr;
   gridlistPtr = ( tGridList< NodeType > * ) listPtr;
   assert( gridlistPtr != 0 );
   if( gridlistPtr->isActiveEmpty() ) curnode = listPtr->first;
   else if( gridlistPtr->isBoundEmpty() ) curnode = 0;
   else curnode = gridlistPtr->lastactive->next;
   if( curnode != 0 ) return 1;
   else return 0;
}


/**************************************************************************\
**
**  tGridListIter::FirstBoundaryP
**
**  Moves the iterator to the first boundary node and returns a pointer
**  to the data at that location.
**
\**************************************************************************/
template< class NodeType >   //tGridListIter
NodeType* tGridListIter< NodeType >::
FirstBoundaryP()
{
   tGridList< NodeType > *gridlistPtr;
   gridlistPtr = ( tGridList< NodeType > * ) listPtr;
   assert( gridlistPtr != 0 );
   if( gridlistPtr->isActiveEmpty() ) curnode = listPtr->first;
   else if( gridlistPtr->isBoundEmpty() ) curnode = 0;
   else curnode = gridlistPtr->lastactive->next;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}


/**************************************************************************\
**
**  tGridListIter::LastActiveP
**
**  Moves the iterator to the last active node and returns a pointer
**  to the data at that location.
**
\**************************************************************************/
template< class NodeType >   //tGridListIter
NodeType *tGridListIter< NodeType >::
LastActiveP()
{
   tGridList< NodeType > *gridlistPtr;
   gridlistPtr = ( tGridList< NodeType > * ) listPtr;
   assert( gridlistPtr != 0 );
   curnode = gridlistPtr->lastactive;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}


/**************************************************************************\
**
**  tGridListIter::IsActive
**
**  Indicates whether the current item is on the active portion of the
**  list, returning 1 if so, 0 if not. Assumes NodeType has a member
**  function getBoundaryFlag.
**
\**************************************************************************/
template< class NodeType >   //tGridListIter
int tGridListIter< NodeType >::
IsActive()
{
   int act;
   if( curnode!=0 )
   {
      assert( curnode->getDataPtr()!=0 );
      act = curnode->getDataRef().getBoundaryFlag();
      if( act == kNonBoundary ) return 1;
   }
   return 0;
}

