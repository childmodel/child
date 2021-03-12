//-*-c++-*-

/**************************************************************************/
/**
 **  @file
 **  @brief Foward declarations for tList.h
 **
 **  AD - March 2004
 **
 **  $Id: tListFwd.h,v 1.2 2004-03-25 17:27:48 childcvs Exp $
 */
/**************************************************************************/

#ifndef TLISTFW_H
#define TLISTFW_H

template< class NodeType > class tListNodeBasic;
template< class NodeType, class ListNodeType = tListNodeBasic<NodeType> > class tList;
template< class NodeType, class ListNodeType = tListNodeBasic<NodeType> > class tListIter;
template< class NodeType, class ListNodeType = tListNodeBasic<NodeType> > class tMeshList;
template< class NodeType, class ListNodeType = tListNodeBasic<NodeType> > class tMeshListIter;

#endif
