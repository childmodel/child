//-*-c++-*-

/*************************************************************************/
/**
 **  @file tOutput.h
 **  @brief Header file for classes tOutput and tLOutput
 **
 **  The tOutput class handles output of triangulated mesh data to
 **  files. The class handles only output of mesh data (nodes,
 **  edges, and triangles); output of additional data (e.g., water or
 **  sediment flow) can be handled by classes derived from tOutput.
 **
 **  tOutput provides functions to open and initialize output files and
 **  write output at a specified time in a simulation. The class is
 **  templated in order to allow for a pointer to a templated
 **  tMesh object.
 **
 **  To handle output of application-specific data, one can create a
 **  class inherited from tOutput and overload its virtual
 **  WriteNodeData function to output the additional data.
 **
 **  Note that in the present version, the files tOutput.h/.cpp
 **  contain the inherited class tLOutput which handles output for
 **  the CHILD model. In the future, such inherited classes will be
 **  kept in separate files to preserve the generality of tOutput.
 **
 **  Recent modifications:
 **    - 1/00: GT added vegofs for output of vegetation cover
 **    - 6/01: GT added chanwidthofs for output of channel widths
 **      (only when non-regime hydraulic geometry model used)
 **    - 7/03: AD added tOutputBase and tTSOutputImp
 **    - 8/03: AD Random number generator handling
 **
 **  $Id: tOutput.h,v 1.48 2003-10-15 09:19:28 childcvs Exp $
 */
/*************************************************************************/

#ifndef TOUTPUT_H
#define TOUTPUT_H

#if !defined(HAVE_NO_NAMESPACE)
# include <fstream>
using namespace std;
#else
# include <fstream.h>
#endif
#include "../MeshElements/meshElements.h"
#include "../tInputFile/tInputFile.h"
#include "../tMesh/tMesh.h"


/**************************************************************************/
/**
 ** @class tOutputBase
 **
 ** Class tOutputBase is use to as base class. It contains common utilities.
 ** The constructor is protected so that the class cannot be used directly.
 **
 */
/**************************************************************************/
template< class tSubNode >
class tOutputBase
{
  tOutputBase(const tOutputBase&);
  tOutputBase& operator=(const tOutputBase&);
protected:
  tOutputBase( tMesh<tSubNode> * meshPtr, const tInputFile &infile );
  virtual ~tOutputBase() { m = 0; }

protected:
  enum{ kMaxNameSize = 80 };
  tMesh<tSubNode> * m;          // ptr to mesh (for access to nodes, etc)
  char baseName[kMaxNameSize];  // name of output files

  void CreateAndOpenFile( ofstream * theOFStream, const char * extension ) const;
  // write time/number of element
  static void WriteTimeNumberElements( ofstream &, double, int );
};

/**************************************************************************/
/**
 ** @class tOutput
 **
 ** Class tOutput handles output of mesh data (nodes, edges, and
 ** triangles). The constructor opens and initializes the files;
 ** the WriteOutput function writes basic mesh data and calls the
 ** virtual function WriteNodeData to write any application-specific
 ** data.
 **
 */
/**************************************************************************/
template< class tSubNode >
class tOutput : public tOutputBase<tSubNode>
{
  tOutput(const tOutput&);
  tOutput& operator=(const tOutput&);
public:
  tOutput( tMesh<tSubNode> * meshPtr, const tInputFile &infile);
  void WriteOutput( double time );

private:
  ofstream nodeofs;             // output file for node data
  ofstream edgofs;              // output file for edge data
  ofstream triofs;              // output file for triangle data
  ofstream zofs;                // output file for node "z" data
  ofstream vaofs;               // output file for Voronoi areas

protected:
  bool CanonicalNumbering;      // Output in canonical order

  virtual void WriteNodeData( double time );

private:
  // renumber in list order
  void RenumberIDInListOrder();
  // use to ensure a canonical ordering
  static int orderRNode(const void*, const void*);
  static int orderREdge(const void*, const void*);
  static int orderRTriangle(const void*, const void*);
  void RenumberIDCanonically();
  // write an individual record
  inline void WriteNodeRecord( tNode * );
  inline void WriteEdgeRecord( tEdge * );
  inline void WriteTriangleRecord( tTriangle const *);
};


template< class tSubNode >
inline void tOutput<tSubNode>::WriteNodeRecord( tNode *cn )
{
  nodeofs << cn->getX() << ' ' << cn->getY() << ' '
	  << cn->getEdg()->getID() << ' '
	  << cn->getBoundaryFlag() << '\n';
  zofs << cn->getZ() << '\n';
  vaofs << cn->getVArea() << '\n';
}

template< class tSubNode >
inline void tOutput<tSubNode>::WriteEdgeRecord( tEdge *ce )
{
  edgofs << ce->getOriginPtrNC()->getID() << ' '
	 << ce->getDestinationPtrNC()->getID() << ' '
	 << ce->getCCWEdg()->getID() << '\n';
}

template< class tSubNode >
inline void tOutput<tSubNode>::WriteTriangleRecord( tTriangle const *ct)
{
  const size_t index[] = {ct->index()[0], ct->index()[1], ct->index()[2]};
  triofs
    << ct->pPtr(index[0])->getID() << ' '
    << ct->pPtr(index[1])->getID() << ' '
    << ct->pPtr(index[2])->getID() << ' '
    << (ct->tPtr(index[0]) ? ct->tPtr(index[0])->getID() : -1) << ' '
    << (ct->tPtr(index[1]) ? ct->tPtr(index[1])->getID() : -1) << ' '
    << (ct->tPtr(index[2]) ? ct->tPtr(index[2])->getID() : -1) << ' '
    << ct->ePtr(index[0])->getID() << ' '
    << ct->ePtr(index[1])->getID() << ' '
    << ct->ePtr(index[2])->getID() << '\n';
}

/**************************************************************************/
/**
 ** @class tLOutput
 **
 ** Class tLOutput handles application-specific data for the CHILD model.
 ** The constructor creates additional output files, and the overloaded
 ** WriteNodeData function writes the data to files.
 ** (TODO: move to separate file)
 **
 ** Modifications:
 **  - 2/02 added output streams tauofs and qsofs for shear stress and
 **    sed flux, resp. (GT)
 **  - 7/03 added tTSOutputImp to separate time series output (AD)
 **
 */
/**************************************************************************/
template< class tSubNode > class tTSOutputImp;

template< class tSubNode >
class tLOutput : public tOutput<tSubNode>
{
  tLOutput(const tLOutput&);
  tLOutput& operator=(const tLOutput&);
public:
  tLOutput( tMesh<tSubNode> * meshPtr, const tInputFile &infile, tRand *);
  virtual ~tLOutput();
  void WriteTSOutput();
  bool OptTSOutput() const;

protected:
  virtual void WriteNodeData( double time );
private:
  ofstream randomofs;  // Random number generator state
  ofstream drareaofs;  // Drainage areas
  ofstream netofs;     // Downstream neighbor IDs
  ofstream slpofs;     // Slopes in the direction of flow
  ofstream qofs;       // Discharge
  ofstream layofs;     // Layer info
  ofstream texofs;     // Texture info
  ofstream vegofs;     // Vegetation cover %
  ofstream flowdepofs; // Flow depth
  ofstream chanwidthofs; // Channel width
  ofstream flowpathlenofs;  // Flow path length
  ofstream tauofs;     // Shear stress
  ofstream qsofs;      // Sed flux

  tTSOutputImp<tSubNode> *TSOutput;  // Time Series output
  tRand *rand;

  int counter;

  inline void WriteActiveNodeData( tSubNode * );
  inline void WriteAllNodeData( tSubNode * );
};

// Write data, including layer info
template< class tSubNode >
inline void tLOutput<tSubNode>::WriteActiveNodeData( tSubNode *cn )
{
  assert( cn!=0 );
  drareaofs << cn->getDrArea() << '\n';
  if( cn->getDownstrmNbr() )
    netofs << cn->getDownstrmNbr()->getID() << '\n';
  layofs << ' ' << cn->getNumLayer() << '\n';
  int i=0;
  while(i<cn->getNumLayer()){
    layofs << cn->getLayerCtime(i) << ' ' << cn->getLayerRtime(i) << ' '
	   << cn->getLayerEtime(i) << '\n'
	   << cn->getLayerDepth(i) << ' ' << cn->getLayerErody(i) << ' '
	   << cn->getLayerSed(i) << '\n';
    int j=0;
    while(j<cn->getNumg()){
      layofs << cn->getLayerDgrade(i,j) << ' ';
      j++;
    }
    layofs << '\n';
    i++;
  }
}

// Write discharge, vegetation, & texture data, etc.
template< class tSubNode >
inline void tLOutput<tSubNode>::WriteAllNodeData( tSubNode *cn )
{
  slpofs << (cn->getBoundaryFlag() == kNonBoundary ? cn->calcSlope():0.) << '\n';
  qofs << cn->getQ() << '\n';
  if( vegofs.good() ) vegofs << cn->getVegCover().getVeg() << '\n';
  if( flowdepofs.good() )
    flowdepofs << cn->getHydrDepth() << '\n';
  if( chanwidthofs.good() )
    chanwidthofs << cn->getHydrWidth() << '\n';
  if( cn->getNumg()>1 ) // temporary hack TODO
    {
      texofs << cn->getLayerDgrade(0,0)/cn->getLayerDepth(0) << '\n';
    }
  if( flowpathlenofs.good() )
    flowpathlenofs << cn->getFlowPathLength() << '\n';
  tauofs << cn->getTau() << '\n';
  if( qsofs.good() ) qsofs << cn->getQs() << '\n';
}


/*
** The following is designed to allow for compiling under the Borland-style
** template instantiation used by the Linux/GNU and Solaris versions of GCC
*/
#include "../Template_model.h"
#ifdef CHILD_TEMPLATE_IN_HEADER
# include "tOutput.cpp"
#endif

#endif
