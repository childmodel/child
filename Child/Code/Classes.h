/****************************************************************************\
**
**   Classes.h: Declarations of classes used in CHILD model
**
**   $Id: Classes.h,v 1.2 2000-03-24 16:53:13 gtucker Exp $
\****************************************************************************/

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
template< class tSubNode > class tListInputData;
class tListIFStreams;
template< class tSubNode > class tListOutputData;
class tListOFStreams;

class tEquilibCheck;
class tSedTransPwrLaw;
class tBedErodePwrLaw;
class tStreamNet;
class tStreamMeander;
class tStorm;
class tRunTimer;

#endif
