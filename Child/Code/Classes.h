//-*-c++-*- 

/****************************************************************************/
/**
**   @file Classes.h
**   @brief Declarations of classes used in CHILD model
**
**   $Id: Classes.h,v 1.8 2004-03-24 15:20:44 childcvs Exp $
*/
/****************************************************************************/

#ifndef CLASSES_H
#define CLASSES_H

/** CLASSES ********************************************************/
class tTriangle;
class tNode;
class tLNode;
class tEdge;
template< class tSubNode > class tMesh;
class tChannel;
class tRegolith;
class tErode;
class tDeposit;
class tMeander;
class tBedrock;
class tSurface;

template< class T > class tArray;
template< class NodeType > class tListNode;
template< class NodeType > class tList;
template< class NodeType > class tMeshList;
template< class NodeType > class tPtrListNode;
template< class NodeType > class tPtrList;
template< class NodeType > class tListIter;
template< class NodeType > class tMeshListIter;
template< class NodeType > class tPtrListIter;

class tInputFile;
template< class tSubNode > class tListInputDataMesh;

class tEquilibCheck;
class tSedTransPwrLaw;
class tBedErodePwrLaw;
class tStreamNet;
class tStreamMeander;
class tStorm;
class tRunTimer;

#endif
