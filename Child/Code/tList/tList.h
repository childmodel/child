//-*-c++-*- 

/**************************************************************************\
**
**  tList.h: Header file for classes tList, tListNode, and tListIter.
**
**  A tList is an object that implements a general linked list NodeType
**  objects, where NodeType can be any type (double, int, other objects,
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
**  tListNode objects are the nodes on the list; each contains an instance
**  of the given data type (double, int, class, etc) and a pointer to the
**  next node in the list.
**
**  A tListIter is an iterator for the linked list tList objects (and their
**  descendants). Its services include fetching data from the current entry
**  on the list, advancing to the next or previous item on the list, etc.
**
**  See also tPtrList, tArray, tMatrix
**
**  Changes:
**    - GT added currentItem member and routines FirstP and NextP to
**      track position on list w/o an iterator, 1/22/99
**    - moved all functions into .h file and inlined them (GT 1/20/00)
**
**  $Id: tList.h,v 1.29 2002-08-12 15:41:16 arnaud Exp $
\**************************************************************************/

#ifndef TLIST_H
#define TLIST_H

template< class NodeType > class tList;
template< class NodeType > class tListIter;
template< class NodeType > class tMeshList;
template< class NodeType > class tMeshListIter;

/**************************************************************************\
** class tListNode ********************************************************
**
** Class tListNode represents the items (or "nodes") on the list. Each
** tListNode object has two parts: the data (of type NodeType) and a
** pointer to the next item on the list. Capabilities include copy
** construction (from either another tListNode or a NodeType),
** returning a pointer or reference to either the data or the tListNode
** itself, and assignment and equality/inequality operations.
**
\**************************************************************************/
template< class NodeType >
class tListNode
{
    friend class tList< NodeType >;
    friend class tMeshList< NodeType >;
    friend class tListIter< NodeType >;
    friend class tMeshListIter< NodeType >;
public:
    tListNode();                                // default constructor
    tListNode( const tListNode< NodeType > & ); // copy constructor #1
    tListNode( const NodeType & );              // copy constructor #2
    const tListNode< NodeType >
    &operator=( const tListNode< NodeType > & );           // assignment
    int operator==( const tListNode< NodeType > & ) const; // equality
    int operator!=( const tListNode< NodeType > & ) const; // inequality
     /*set*/
    NodeType getDataNC();                     // returns copy of data item
    NodeType &getDataRefNC();                 // returns modifiable ref to data
    NodeType *getDataPtrNC();                 // returns modifiable ptr to data
    tListNode< NodeType > * getNextNC() const;// returns ptr to next list node
     /*get*/
    NodeType getData() const;                     // returns const copy of data
    const NodeType &getDataRef() const;           // returns const ref to data
    const NodeType *getDataPtr() const;           // returns const ptr to data
    const tListNode< NodeType > * getNext() const;// returns const ptr to next
   
protected:
    NodeType data;               // data item
    tListNode< NodeType > *next; // ptr to next node on list (=0 if end)
};


/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tListNode< NodeType >
**
**         tListNodes contain a data item and a pointer to the next
**         tListNode in the tList (or tGridList). The data item may
**         be of any type, specified in the template brackets.
**
**         Some of the functions for retrieving the data are duplicated
**         by tListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/

/**************************************************************************\
**
**  tListNode constructors:
**
**  Default constructor: sets next to null
**  Copy constructor #1: makes a copy of a given tListNode
**  Copy constructor #2: fills in data item w/ copy of given NodeType
**
\**************************************************************************/

//default constructor
template< class NodeType >                     //tListNode
inline tListNode< NodeType >::
tListNode() :
  next(0)
{
     //data = 0;
}

//copy constructor with data reference
template< class NodeType >                     //tListNode
inline tListNode< NodeType >::
tListNode( const tListNode< NodeType > &original ) :
  data(original.data),
  next(original.next)
{}

//value (by reference) constructor 
template< class NodeType >                     //tListNode
inline tListNode< NodeType >::
tListNode( const NodeType &info ) :
  data(info),
  next(0)
{}


/**************************************************************************\
**
**  tListNode overloaded operators:
**
**  Assignment: makes a copy (including next ptr)
**  Equality: compares both data contents and next ptr
**  Inequality: compares both data contents and next ptr
**
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >                     //tListNode
inline const tListNode< NodeType > &tListNode< NodeType >::
operator=( const tListNode< NodeType > &right )
{
   if( &right != this )
   {
        //delete data;
        //data = new NodeType;
      assert( &data != 0 );
      data = right.data;
      next = right.next;
   }
   return *this;
}

//overloaded equality operator:
template< class NodeType >                     //tListNode
inline int tListNode< NodeType >::
operator==( const tListNode< NodeType > &right ) const
{
   if( next != right.next ) return 0;
   if( &data != &(right.data) ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                     //tListNode
inline int tListNode< NodeType >::
operator!=( const tListNode< NodeType > &right ) const
{
   return ! operator==(right);
}


/**************************************************************************\
**
**  tListNode "get" functions:
**  (note: to "set" an item, use non-const "get")
**
**  getDataNC: returns a non-const (modifiable) copy of data
**  getDataRefNC: returns a non-const (modifiable) reference to data
**  getDataPtrNC: returns a non-const (modifiable) pointer to data
**  getNextNC: returns non-const ptr to next item on list
**  getData: returns const copy of data
**  getDataRef: returns const reference to data
**  getDataPtr: returns const ptr to data
**  getNext: returns const ptr to next item on list
**
\**************************************************************************/

//set data by returning non-const
template< class NodeType >                     //tListNode
inline NodeType tListNode< NodeType >::
getDataNC() {return data;}

template< class NodeType >                     //tListNode
inline NodeType &tListNode< NodeType >::
getDataRefNC() {return data;}

template< class NodeType >                     //tListNode
inline NodeType *tListNode< NodeType >::
getDataPtrNC() {return &data;}

template< class NodeType >                     //tListNode
inline tListNode< NodeType > * tListNode< NodeType >::
getNextNC() const {return next;}

//return data by value
template< class NodeType >                     //tListNode
inline NodeType tListNode< NodeType >::
getData() const {return data;}

//return data by reference
template< class NodeType >                     //tListNode
inline const NodeType &tListNode< NodeType >::
getDataRef() const {return data;}

//return data by pointer
template< class NodeType >                     //tListNode
inline const NodeType *tListNode< NodeType >::
getDataPtr() const {return &data;}

//return next pointer
template< class NodeType >                     //tListNode
inline const tListNode< NodeType > * tListNode< NodeType >::
getNext() const {return next;}



/**************************************************************************\
** class tList ************************************************************
**
** Class tList implements a linked list. The class includes pointers to
** the first, last, and current list nodes (see tListNode).
**
\**************************************************************************/
template< class NodeType >
class tList
{
    friend class tListIter< NodeType >;
    friend class tMeshListIter< NodeType >;
    tList(const tList&);
public:
    tList();                            // default constructor
    tList( const tList< NodeType > * ); // copy constructor
    ~tList();                           // destructor
    const tList< NodeType >
    &operator=( const tList< NodeType > & );           // assignment
    int operator==( const tList< NodeType > & ) const; // equality
    int operator!=( const tList< NodeType > & ) const; // inequality
    void insertAtFront( const NodeType & ); // puts copy of item at list front
    void insertAtBack( const NodeType & );  // puts copy of item at list back
    void insertAtNext( const NodeType &, tListNode< NodeType > * );
    void insertAtPrev( const NodeType &, tListNode< NodeType > * );
    int removeFromFront( NodeType & ); // removes 1st item, puts it in ref
    int removeFromBack( NodeType & );  // removes last item, puts it in ref
    int removeNext( NodeType &, tListNode< NodeType > * );
    int removePrev( NodeType &, tListNode< NodeType > * );
    void Flush();        // clears and reinitializes list
    int isEmpty() const; // returns 1 is list is empty, 0 otherwise
#ifndef NDEBUG
    void print() const;  // prints contents of list - DEBUG ONLY
#endif
    int getSize() const; // returns # of items on list
    tListNode< NodeType  > * getFirst() const; // returns ptr to 1st list node
    tListNode< NodeType  > * getLast() const;  // returns ptr to last list node
    NodeType * FirstP();  // returns ptr to 1st data item & sets current to 1st
    NodeType * NextP();   // moves to next node and returns ptr to data item
    void moveToBack( tListNode< NodeType > *  );  // move given node to back
    void moveToFront( tListNode< NodeType > *  ); // move given node to front
    void makeCircular();   // makes list circular (last points to first)
    NodeType getIthData( int ) const;     // rtns copy of given item #
    const NodeType *getIthDataPtr( int ) const; // rtns ptr to given item #
    const NodeType &getIthDataRef( int ) const; // rtns ref to given item #
    NodeType getIthDataNC( int ) const;     // rtns modifiable copy of item #
    NodeType *getIthDataPtrNC( int ) const; // rtns modifiable ptr to item #
    NodeType &getIthDataRefNC( int ) const; // rtns modifiable ref to item #
    tListNode< NodeType > * getListNode( NodeType * ); // rtns ptr to node #

#ifndef NDEBUG
    void DebugTellPtrs();
#endif
    
protected:
    int nNodes;                          // # of items on list
    tListNode< NodeType > * first;       // ptr to first node
    tListNode< NodeType > * last;        // ptr to last node
    tListNode< NodeType > * currentItem; // ptr to current item
    tListNode< NodeType > * getNewNode( const NodeType & ); // makes new node
};



/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tList< NodeType >
**
**         tList is a linked list of tListNodes. Functions to add and
**         remove items to/from list, etc.
**
**         Again, some functions are duplicated by tListIter.
**
**         Created: fall, 97, SL
**
\**************************************************************************/

/**************************************************************************\
**
**  tList constructors & destructor:
**
**  Default constructor: initializes all values to 0 (empty list)
**  Copy constructor: makes a complete copy of another tList
**  Destructor: deletes all nodes on list
**
\**************************************************************************/

//default constructor
template< class NodeType >                         //tList
inline tList< NodeType >::tList() :
  nNodes(0), first(0), last(0), currentItem(0)
{
     //cout << "list instantiated" << first << endl;
}

//copy constructor
template< class NodeType >                         //tList
inline tList< NodeType >::
tList( const tList< NodeType > *original ) :
  nNodes(0)
{
   int i;

   assert( original != 0 );
     //nNodes = original->nNodes;
   tListNode<NodeType> * current = original->first;
   for( i=0; i<original->nNodes; i++ )
   {
      insertAtBack( current->data );
      current = current->next;
   }
   assert( nNodes == original->nNodes );
   cout << "list copy instantiated" << first << endl;
   current = first;
   
}

//destructor
template< class NodeType >                         //tList
inline tList< NodeType >::
~tList()
{
   if( !isEmpty() )
   {
        //cout<<"Destroying nodes ... "<<endl;
      tListNode<NodeType > * current = first, * temp;
      while( current != 0 )
      {
         temp = current;
         //cout<<temp->data<<endl;
         current = current->next;
         delete temp;
      }
   }
     //cout<<"All nodes destroyed"<<endl<<endl;
}

/**************************************************************************\
**
**  tList overloaded operators:
**
**  Assignment: clears the list and makes a copy of the right-hand list
**  Equality: returns TRUE if first & last pointers are identical and
**            nNodes is the same. Note that two lists with identical
**            contents are still not considered equal! (TODO--makes sense?)
**  Inequality: opposite of equality
**
\**************************************************************************/

//overloaded assignment operator
template< class NodeType >                         //tList
inline const tList< NodeType > &tList< NodeType >::
operator=( const tList< NodeType > &right )
{
   if( this != &right )
   {
      // create an equivalent object rather than pointers to the original
      Flush();
      tListNode< NodeType > *cn = right.first;
      if( cn != 0 )
      {
          insertAtBack( cn->data );
          for( cn = cn->next; cn != last->next; cn = cn->next )
             insertAtBack( cn->data );
      }
        /*first = right.first;
      last = right.last;*/
      assert( nNodes == right.nNodes );
   }
     //cout << "list assigned" << first << endl;
   return *this;
}

//overloaded equality operator:
template< class NodeType >                         //tList
inline int tList< NodeType >::
operator==( const tList< NodeType > &right ) const
{
   if( nNodes != right.nNodes ) return 0;
   if( first != right.first ) return 0;
   if( last != right.last ) return 0;
   return 1;
}

//overloaded inequality operator:
template< class NodeType >                         //tList
inline int tList< NodeType >::
operator!=( const tList< NodeType > &right ) const
{
   if( nNodes != right.nNodes ) return 1;
   if( first != right.first ) return 1;
   if( last != right.last ) return 1;
   return 0;
}


/**************************************************************************\
**
**  tList::getNewNode
**
**  Creates a new tListNode and returns a pointer to it. Used by list
**  insertion routines (see below); not publically accessible.
**
\**************************************************************************/
template< class NodeType >                         //tList
inline tListNode< NodeType > * tList< NodeType >::
getNewNode( const NodeType &value )
{
   tListNode< NodeType > * ptr =
       new tListNode< NodeType >( value );
   assert( ptr != 0 );
   //cout << "new node created" << endl;
   nNodes++;
   return ptr;
}


/**************************************************************************\
**
**  tList: list insertion and removal routines
**
**  A collection of routines to add nodes to the list, remove them from
**  the list, and change the position of nodes on the list. In the case
**  of insertion routines, a new copy of the given node is created first.
**  The list removal routines take a reference to a NodeType variable;
**  on return, this variable contains a copy of the data item removed.
**  List removal routines return TRUE if there was an item to remove;
**  FALSE otherwise.
**
**    insertAtFront: make a new node containing _value_ and put at top
**    insertAtBack: make a new node containing _value_ and put at bottom
**    insertAtNext: make new node containing _value_ and place it on the
**                  list after _prev_
**    insertAtPrev: make new node containing _value_ and place it on the
**                  list before _node_
**    removeFromFront: remove 1st item on list and place a copy in _value_
**    removeFromBack: remove last item on list and place a copy in _value_
**    removeNext: remove the node following node _ptr_ and place a copy
**                in _value_
**    removePrev: remove the node before node _ptr_ and place a copy
**                in _value_
**
\**************************************************************************/

//insert at front
template< class NodeType >                         //tList
inline void tList< NodeType >::
insertAtFront( const NodeType &value )
{
   //cout << "ADD NEW NODE TO LIST AT FRONT" << endl;
   
   tListNode< NodeType > *newPtr = getNewNode( value );
   if( isEmpty() ) first = last = currentItem = newPtr;
   else
   {
      newPtr->next = first;
      if( last->next == first ) last->next = newPtr;
      first = newPtr;
   }
     //nNodes++; now handled by getNewNode
}

//insert at back
template< class NodeType >                         //tList
inline void tList< NodeType >::
insertAtBack( const NodeType &value )
{
   tListNode< NodeType > * newPtr = getNewNode( value );
   //cout << "add new node to list in back" << endl;
   if( isEmpty() )
   {
      first = last = currentItem = newPtr;
   }
   else
   {
      newPtr->next = last->next;
      last->next = newPtr;
      last = newPtr;
   }

   //cout << "data is " << newPtr->data << endl;
}

//insert at next spot in list
template< class NodeType >                         //tList
inline void tList< NodeType >::
insertAtNext( const NodeType &value, tListNode< NodeType > * prev )
{
   if( prev != 0 )
   {
      if( prev == last )
      {
         insertAtBack( value );
         return;
      }
      tListNode< NodeType > * newPtr = getNewNode( value );
      newPtr->next = prev->next;
      prev->next = newPtr;
   }
}

//insert at previous spot in list
template< class NodeType >                         //tList
inline void tList< NodeType >::
insertAtPrev( const NodeType &value, tListNode< NodeType > * node )
{
   tListNode< NodeType > * prev;
   if( node != 0 )
   {
      if( node == first )
      {
         insertAtFront( value );
         return;
      }
      tListNode< NodeType > * newPtr = getNewNode( value );
      for( prev = first; prev->next != node; prev = prev->next ); 
      newPtr->next = prev->next;
      prev->next = newPtr;
   }
}

//delete from front
template< class NodeType >                         //tList
inline int tList< NodeType >::
removeFromFront( NodeType &value )
{
   if( isEmpty() ) return 0;
   else
   {
      tListNode< NodeType > * temp = first;
      if( first == last ) first = last = currentItem = 0;
      else
      {
         if( last->next == first ) last->next = first->next;
         if( currentItem==first ) currentItem = first->next;
         first = first->next;
      }
      value = temp->data;
      delete temp;
      nNodes--;
      return 1;
   }
}

//delete from back
template< class NodeType >                         //tList
inline int tList< NodeType >::
removeFromBack( NodeType &value )
{
   if( isEmpty() ) return 0;
   else
   {
      tListNode< NodeType > * temp = last;
      if( first == last ) first = last = currentItem = 0;
      else
      {
         tListNode< NodeType > * nexttolast = first;
         while( nexttolast->next != last ) nexttolast = nexttolast->next;
         nexttolast->next = last->next;
         if( currentItem==last ) currentItem = nexttolast;
         last = nexttolast;
      }
      value = temp->data;
      delete temp;
      nNodes--;
      return 1;
   }
}

//delete next node
template< class NodeType >                         //tList
inline int tList< NodeType >::
removeNext( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr->next == 0 ) return 0;
   if( ptr == 0 ) return 0;
   if( ptr->next == last ) return removeFromBack( value );
   else if( ptr->next == first ) return removeFromFront( value );
   //if( ptr == last ) return 0;
   tListNode< NodeType > * temp = ptr->next;
   if( currentItem == temp ) currentItem = ptr;
   ptr->next = ptr->next->next;
   value = temp->data;
   delete temp;
   nNodes--;
   return 1;
}

//delete previous node
template< class NodeType >                         //tList
inline int tList< NodeType >::
removePrev( NodeType &value, tListNode< NodeType > * ptr )
{
   if( ptr == 0 ) return 0;
   if( ptr == first && last->next == 0 ) return 0;
   if( ptr == first ) return removeFromBack( value );
   tListNode< NodeType > * temp, *prev;
   for( temp = first; temp->next->next != ptr; temp = temp->next );
   prev = temp;
   temp = temp->next;
   if( temp == first ) return removeFromFront( value );
   prev->next = prev->next->next;
   value = temp->data;
   if( currentItem == temp ) currentItem = prev;
   delete temp;
   nNodes--;
   return 1;
}


/**************************************************************************\
**
**  tList::Flush
**
**  Deletes all nodes on list.
**
\**************************************************************************/
template< class NodeType >                         //tList
inline void tList< NodeType >::
Flush()
{
   tListNode<NodeType > * current = first, * temp;

   if( !isEmpty() )
   {
      //cout<<"Destroying nodes ... "<<endl;
      //while( removeFromBack( temp ) );
      while( current != 0 )
      {
         temp = current;
         //cout<<temp->data<<endl;
         current = current->next;
         delete temp;
      }
      first = last = currentItem = 0;
   }
   assert( isEmpty() );
   nNodes = 0;
     //cout<<"All nodes destroyed"<<endl<<endl;
}


/**************************************************************************\
**
**  tList::isEmpty
**
**  Returns TRUE if first points to null; FALSE otherwise.
**
\**************************************************************************/
//empty?
template< class NodeType >                         //tList
inline int tList< NodeType >::
isEmpty() const
{
   if( first == 0 )
   {
        //cout << "list is empty" << endl;
      return 1;
   }
   else
   {
        //cout << "list is not empty" << endl;
      return 0;
   }
}

//display list contents -- for debugging only
/*
template< class NodeType >                         //tList
void tList< NodeType >::
print() const
{
   if( isEmpty() )
   {
      cout<<"The list is empty"<<endl<<endl;
      return;
   }
   tListNode< NodeType > * current = first;
   cout<<"The list is: ";
   while( current != 0 )
   {
      //cout<< current->data <<' ';
      current = current->next;
   }
   cout<<endl<<endl;
}*/


/**************************************************************************\
**
**  tList "get" functions:
**
**  getSize: returns # of items on list
**  getFirst: returns const ptr to first tListNode
**  getLast: returns const ptr to last tListNode
**  FirstP: returns ptr to first data item and sets currentItem to first
**  NextP: returns ptr to data item following currentItem and advances
**         currentItem to the next one on the list
**  getIthData: returns const copy of list data item number _num_
**  getIthDataRef: returns const reference to list data item number _num_
**  getIthDataRef: returns const ptr to list data item number _num_
**  getIthDataNC: returns non-const copy of list data item number _num_
**  getIthDataRefNC: returns non-const ref to list data item number _num_
**  getIthDataRefNC: returns non-const ptr to list data item number _num_
**  (see also getListNode, below)
**
\**************************************************************************/

template< class NodeType >                         //tList
inline int tList< NodeType >::
getSize() const {return nNodes;}

template< class NodeType >                         //tList
inline tListNode< NodeType > * tList< NodeType >::
getFirst() const {return first;}

template< class NodeType >                         //tList
inline tListNode< NodeType > * tList< NodeType >::
getLast() const {return last;}

// Added by gt 1/99
template< class NodeType >                         //tList
inline NodeType * tList< NodeType >::
FirstP() 
{
   assert( first!=0 );
   currentItem = first;
   return &first->data;
}

// Added by gt 1/99
template< class NodeType >                         //tList
inline NodeType * tList< NodeType >::
NextP() 
{
   assert( currentItem!=0 );
   currentItem = currentItem->next;
   if( currentItem!=0 )
       return &currentItem->data;
   else return 0;
}

template< class NodeType >                         //tList
inline NodeType tList< NodeType >::
getIthData( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
     if(num>= nNodes)
         {
            cout<<"using an index which is too large"<<endl;
            cout<<"you have "<<nNodes<<endl;
            cout<<"you wanted list member number "<<num<<endl;
         }
//     if(num<0)
//     {
//        cout<<"using a negative index"<<endl;
//        cout<<"you have "<<nNodes<<endl;
//        cout<<"you wanted list member number "<<num<<endl;
//     }
   //cout<<"num="<<num<<"nNodes="<<nNodes<<endl;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   assert( curPtr!=0 );
   return curPtr->getData();
}

template< class NodeType >                         //tList
inline const NodeType &tList< NodeType >::
getIthDataRef( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataRef();
}

template< class NodeType >                         //tList
inline const NodeType *tList< NodeType >::
getIthDataPtr( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0;i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataPtr();
}

//set
template< class NodeType >                         //tList
inline NodeType tList< NodeType >::
getIthDataNC( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getData();
}

template< class NodeType >                         //tList
inline NodeType &tList< NodeType >::
getIthDataRefNC( int num ) const
{
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataRefNC();
}

template< class NodeType >                         //tList
inline NodeType *tList< NodeType >::
getIthDataPtrNC( int num ) const
 {
   int i;
   tListNode< NodeType > * curPtr;
   assert( num >= 0 && num < nNodes );
   for( curPtr = first, i = 0; i<num; i++ )
   {
      curPtr = curPtr->next;
   }
   return curPtr->getDataPtrNC();
}


/**************************************************************************\
**
**  tList::getListNode
**
**  Finds and returns the list node containing the data pointed to by
**  desiredDatPtr, or zero if not found.
**
**  Parameters:  desiredDatPtr -- pointer to the data item sought after
**  Returns:  pointer to the tListNode containing desiredDatPtr, or zero
**            if not found
**  Notes: might be safer to implement with a const return type
**  Created: 4/29/98 GT
**
\**************************************************************************/
template< class NodeType >
inline tListNode< NodeType > * tList< NodeType >::
getListNode( NodeType * desiredDatPtr )
{
   tListNode< NodeType > * listnode = first;

   if( listnode==0 ) return 0;
   while( &(listnode->data) != desiredDatPtr )
   {
      listnode = listnode->next;
      if( listnode==0 ) return 0;
   }
   return listnode;
   
}


/**************************************************************************\
**
**  tList "move" functions:
**
**  moveToBack: looks for _mvnode_ on the list and moves it to the back
**              of the list. Assumes that _mvnode_ IS an item on the list.
**  moveToFront: same thing, but moves it to the front
**  
\**************************************************************************/

template< class NodeType >                         //tList
inline void tList< NodeType >::
moveToBack( tListNode< NodeType > * mvnode ) 
{
   tListNode< NodeType > * prev;
   if( mvnode != last )
   {  
      if( mvnode == first ) first = first->next;
      else
      {
         for( prev = first; prev->next != mvnode; prev = prev->next );
         prev->next = mvnode->next;
      }
      mvnode->next = last->next;
      last->next = mvnode;
      last = mvnode;
      if( last->next != 0 ) last->next = first;
   }
}

template< class NodeType >                         //tList
inline void tList< NodeType >::
moveToFront( tListNode< NodeType > * mvnode ) 
{
   tListNode< NodeType > * prev;
   if( mvnode != first )
   {
      for( prev = first; prev->next != mvnode; prev = prev->next );
      prev->next = mvnode->next;
      mvnode->next = first;
      first = mvnode;
      if( last == mvnode ) last = prev;
      if( last->next != 0 ) last->next = first;
   }
}



/* GT commented this out 1/99 -- doesn't seem to be used anywhere, and
   is a bit dangerous bcs the actual # of items on list is untouched
template< class NodeType >                         //tList
void tList< NodeType >::
setNNodes( int val ) {nNodes = ( val >= 0 ) ? val : 0;}*/

/**************************************************************************\
**
**  tList::makeCircular
**
**  Converts the list into a circular list by having the last item point
**  to the first.
**
\**************************************************************************/
template< class NodeType >                         //tList
inline void tList< NodeType >::
makeCircular() {last->next = first;}


template< class NodeType >
inline void tList< NodeType >::
DebugTellPtrs()
{
   cout << "first: " << first << endl;
   cout << "last: " << last << endl << flush;
}



/**************************************************************************\
** class tListIter *********************************************************
**
** Helper class for tList. tListIters are "iterators" that walk up and
** down a tList, fetching items. Their chief advantage is that you can
** have multiple iterators on any given list at once, and thus multiple
** access points. Use of iterator classes is discussed by Deitel and
** Deitel, _C++ How to Program_, first edition, Prentice Hall, 1994.
**
** Note that in the current implementation, list items are fetched by
** ID number, which presupposes that the list items have a member function
** getID. This restricts the generality of tList, and should be moved
** to tMeshList. (TODO)
**
\**************************************************************************/
//TO DO: make Get, Where, GetP, refer to place in list rather than use getID()
template< class NodeType >
class tListIter
{
  tListIter& operator=(const tListIter&);
 public:
    tListIter();                      // default constructor
    tListIter(const tListIter&);      // copy constructor
    tListIter( tList< NodeType > & ); // constructor: reference to list
    tListIter( tList< NodeType > * ); // constructor: ptr to list
    ~tListIter();   // destructor
    int First();    // sets position to 1st list item (rtns 0 on failure)
    int Last();     // sets position to last "    "     "
    int Get( int ); // use only if NodeType has member getID()!!
    int Next();     // advances to next item (or 1st if current undefined)
    int Prev();     // moves to previous item (or last if current undef'd)
    int Where();    // use only if NodeType has member getID()!!
    int AtEnd();    // returns 1 if at end of the list
    NodeType &DatRef();  // returns ref to current data item
    NodeType *DatPtr();  // returns ptr to current data item
    tListNode< NodeType > *NodePtr();  // returns ptr to current list node
    void Reset( tList< NodeType > & ); // tells iterator to work on given list
    NodeType * FirstP();  // moves to 1st item and rtns ptr to it
    NodeType * LastP();   // moves to last  "   "  
    NodeType * NextP();   // moves to next  "   "
    NodeType * PrevP();   // moves to previous " "
    NodeType * GetP( int num ); //use only if NodeType has member getID()!!
    
protected:
    tListNode< NodeType > * curnode;  // ptr to current list node
    tList< NodeType > *listPtr;       // ptr to current list
    int counter;                      // current position on list (first=0)
};

/**************************************************************************\
**
**         FUNCTIONS FOR CLASS tListIter
**
**  A tListIter is an iterator for the linked list tList objects (and their
**  descendants). Its services include fetching data from the current entry
**  on the list, advancing to the next or previous item on the list, etc.
**
**  See also tGridList.
**
**  Created: fall, 97, SL.
**  Modifications:  added an "AtEnd" function that signals whether the
**  the iterator has fallen off the end of the list, 11/17/97 gt.
**
\**************************************************************************/

/**************************************************************************\
**
**  tListIter constructors & destructor:
**
**  Default constructor: initializes all values to 0
**  Constructor (reference version): sets listPtr to point to _list_ and
**                                   points curnode to 1st node on list
**  Constructor (pointer version): same thing, but takes a pointer to a
**                                 tList as an argument
**
\**************************************************************************/

template< class NodeType >        //tListIter
inline tListIter< NodeType >::
tListIter() :
  curnode(0),
  listPtr(0),
  counter(0)
{
   //cout << "tListIter()" << endl;
}

template< class NodeType >        //tListIter
inline tListIter< NodeType >::
tListIter(const tListIter& c) :
  curnode(c.curnode),
  listPtr(c.listPtr),
  counter(c.counter)
{
  if (0) //DEBUG
    cout << "tListIter(const tListIter&)" << endl;
}

template< class NodeType >        //tListIter
inline tListIter< NodeType >::
tListIter( tList< NodeType > &list ) :
  curnode(list.first),
  listPtr(&list),
  counter(0)
{
   //assert( curnode != 0 );
   //cout << "tListIter( list )" << endl;
}

template< class NodeType >        //tListIter
inline tListIter< NodeType >::
tListIter( tList< NodeType > *ptr ) :
  curnode(0),
  listPtr(0),
  counter(0)
{
   assert( ptr != 0 );
   listPtr = ptr;
   curnode = ptr->first;
   //assert( curnode != 0 );
     //cout << "tListIter( ptr )" << endl;
   
}

template< class NodeType >        //tListIter
inline tListIter< NodeType >::
~tListIter()
{
   listPtr = 0;
   curnode = 0;
     //cout << "~tListIter()" << endl;
}


/**************************************************************************\
**
**  tListIter::Reset
**
**  Points iterator at the 1st node on _list_ (provides a way of telling
**  an iterator which list to work on).
**
\**************************************************************************/
template< class NodeType >        //tListIter
inline void tListIter< NodeType >::
Reset( tList< NodeType > &list )
{
   assert( &list != 0 );
   listPtr = &list;
   curnode = list.first;
   counter = 0;
}


/**************************************************************************\
**
**  tListIter::First and tListIter::Last
**
**  Move to the first or last item on the current list. Return TRUE if
**  pointing to a valid list item, FALSE otherwise.
**
\**************************************************************************/
template< class NodeType >        //tListIter
inline int tListIter< NodeType >::
First()
{
   assert( listPtr != 0 );
   curnode = listPtr->first;
   counter = 0;
   if( curnode != 0 ) return 1;
   else if( curnode == 0 && listPtr->isEmpty() ) return 1;
   return 0;
}

template< class NodeType >       //tListIter
inline int tListIter< NodeType >::
Last()
{
   assert( listPtr != 0 );
   curnode = listPtr->last;
   counter = -1;
   if( curnode != 0 ) return 1;
     //else if( curnode == 0 && listPtr->isEmpty() ) return 1;
   return 0;
}


/**************************************************************************\
**
**  tListIter::Get
**
**  Move to list item with ID number _num_. Note: assumes that list items
**  have a member function getID()! Returns 1 if found, 0 if not.
**
\**************************************************************************/
template< class NodeType >     //tListIter
inline int tListIter< NodeType >::
Get( int num )
{
   assert( listPtr != 0 );
//     if( num < 0 )
//     {
//        cout << "tListIter::Get(num): num < 0" << endl;
//        //return 0;
//     }
   tListNode< NodeType > *tempnodeptr;
   for( tempnodeptr = listPtr->first, counter = 0;
        counter <= listPtr->nNodes && tempnodeptr != 0;
        tempnodeptr = tempnodeptr->next, counter++ )
   {
      if( tempnodeptr->getDataPtr()->getID() == num ) break;
   }
   if( tempnodeptr == 0 )
   {
      //cout << "tListIter::Get(num): tempnodeptr == 0" << endl;
      return 0;
   }
   if( tempnodeptr->getDataPtr()->getID() != num )
   {
      //cout << "tListIter::Get(num): tempnodeptr->data.getID() != num" << endl;
      return 0;
   }
   curnode = tempnodeptr;
   return 1;
}

   
/**************************************************************************\
**
**  tListIter::Next and tListIter::Prev
**
**  Move to the next or previous item on the current list. Return TRUE if
**  pointing to a valid list item, FALSE otherwise. If we're not 
**  initially pointing to any item, then move to the first or last item,
**  respectively. Both assume we're working on a valid list.
**
\**************************************************************************/
template< class NodeType >        //tListIter
inline int tListIter< NodeType >::
Next()
{
   assert( listPtr != 0 );

   // if current position undefined, move to first node...
   if( curnode == 0 )
   {
      curnode = listPtr->first;
      counter = 0;
      if( curnode != 0 ) return 1;
      else return 0;
   }

   // ...otherwise just move to the next one
   curnode = curnode->next;
   counter++;
   if( curnode != 0 ) return 1;
   else return 0;
}

template< class NodeType >       //tListIter
inline int tListIter< NodeType >::
Prev()
{
   assert( listPtr != 0 );

   // if current position undefined, move to the last one
   if( curnode == 0 )
   {
      curnode = listPtr->last;
      counter = -1; // why -1 and not nNodes? TODO
      if( curnode != 0 ) return 1;
      else return 0;
   }

   // if we're at the first node, the previous one is only defined if we're
   // a circular list, in which case last points to first -- so move to last
   if( curnode == listPtr->first )
   {
      if( listPtr->last->next == 0 ) return 0;
      else
      {
         assert( curnode == listPtr->last->next );
         curnode = listPtr->last;
         counter = -1; // why -1?
         return 1;
      }
   }

   // general case: search through the list until we reach the one before
   // curnode, then set curnode to that one
   tListNode< NodeType > *tempnode;
   //int id = curnode->data.getID();
   for( tempnode = listPtr->first;
        tempnode->next != curnode; //tempnode->next->data.getID() != id;
        tempnode = tempnode->next );
   curnode = tempnode;
   counter--;
   assert( curnode != 0 );
   return 1;
}


/**************************************************************************\
**
**  tListIter::FirstP and tListIter::LastP
**
**  Move to the first or last item on the list and return a pointer to the
**  data, or 0 if first/last item is empty.
**
\**************************************************************************/
template< class NodeType >        //tListIter
inline NodeType * tListIter< NodeType >::
FirstP()
{
   assert( listPtr != 0 );
   curnode = listPtr->first;
   counter = 0;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}
   
template< class NodeType >        //tListIter
inline NodeType * tListIter< NodeType >::
LastP()
{
   assert( listPtr != 0 );
   curnode = listPtr->last;
   counter = 0;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}
   

/**************************************************************************\
**
**  tListIter::NextP and tListIter::PrevP
**
**  Same as Next and Prev, except that the functions return a pointer to
**  the current data item (or 0 if none exists).
**
\**************************************************************************/

template< class NodeType >        //tListIter
inline NodeType * tListIter< NodeType >::
NextP()
{
   assert( listPtr != 0 );
   if( curnode == 0 )
   {
      curnode = listPtr->first;
      counter = 0;
      if( curnode != 0 ) return curnode->getDataPtrNC();
      else return 0;
   }
   curnode = curnode->next;
   counter++;
   if( curnode != 0 ) return curnode->getDataPtrNC();
   else return 0;
}

template< class NodeType >       //tListIter
inline NodeType *tListIter< NodeType >::
PrevP()
{
   assert( listPtr != 0 );
   if( curnode == 0 )
   {
      curnode = listPtr->last;
      counter = -1;
      if( curnode != 0 ) return curnode->getDataPtrNC();
      else return 0;
   }
   if( curnode == listPtr->first )
   {
      if( listPtr->last->next == 0 ) return 0;
      else
      {
         assert( curnode == listPtr->last->next );
         curnode = listPtr->last;
         counter = -1;
         return curnode->getDataPtrNC();
      }
   }
   tListNode< NodeType > *tempnode;
   //int id = curnode->data.getID();
   for( tempnode = listPtr->first;
        tempnode->next != curnode;  //tempnode->next->data.getID() != id;
        tempnode = tempnode->next );
   curnode = tempnode;
   assert( curnode != 0 );
   counter--;
   return curnode->getDataPtrNC();
}


/**************************************************************************\
**
**  tListIter::GetP
**
**  Similar to Get, but returns a pointer to the current data item (or
**  0 if undefined).
**
\**************************************************************************/
template< class NodeType >       //tListIter
inline NodeType * tListIter< NodeType >::
GetP( int num )
{
   assert( listPtr != 0 );
   if( num < 0 ) return 0;
   //cout << "Get: num " << num << "; ";

   return
     (Get(num) != 0 ? curnode->getDataPtrNC() : 0 );
}


/**************************************************************************\
**
**  tListIter::Where
**
**  Returns the ID number of the current data item, or -1 if there is
**  no current data item. Assumes data item has a getID() mbr function!
**
\**************************************************************************/
template< class NodeType >       //tListIter
inline int tListIter< NodeType >::
Where()
{
   if( curnode == 0 ) return -1;
   return curnode->getDataPtr()->getID();
}


/**************************************************************************\
**
**  tListIter::AtEnd
**
**  Returns TRUE if:
**   - the list is empty
**   - the list is non-circular and the current item is null
**   - the list is circular, the current item is the first, and the 
**     counter is nonzero (meaning we've gone all the way through the
**     list and come back to the start)
**
\**************************************************************************/
template< class NodeType >       //tListIter
inline int tListIter< NodeType >::
AtEnd()
{
   if( listPtr->isEmpty() ) return 1;
   if( listPtr->last->next == 0 ) return ( curnode==0 );
   else return ( curnode == listPtr->first && counter != 0 );
   //return curnode==0;
}


/**************************************************************************\
**
**  tListIter::DatRef, DatPtr, and NodePtr
**
**  Returns a non-constant reference or pointer to the current data item
**  or the current list node.
**
\**************************************************************************/
template< class NodeType >       //tListIter
inline NodeType &tListIter< NodeType >::
DatRef()
{
     //if( curnode == 0 ) return 0;
   return curnode->getDataRefNC();
}

template< class NodeType >       //tListIter
inline NodeType *tListIter< NodeType >::
DatPtr()
{
   if( curnode == 0 ) return 0;
   return curnode->getDataPtrNC();
}

template< class NodeType >       //tListIter
inline tListNode< NodeType > *tListIter< NodeType >::
NodePtr()
{
   //assert( curnode != 0 );
   return curnode;
}

#endif
