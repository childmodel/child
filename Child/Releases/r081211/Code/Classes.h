//-*-c++-*- 

/****************************************************************************/
/**
**   @file Classes.h
**   @brief Declarations of classes used in CHILD model
**
**   $Id: Classes.h,v 1.9 2004/03/25 12:14:47 childcvs Exp $
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
#include "tList/tListFwd.h"
template< class NodeType > class tPtrListNode;
template< class NodeType > class tPtrList;
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
