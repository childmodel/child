/*************************************************************************/
/**
 **  @file tOutput.cpp
 **  @brief Functions for output classes tOutput and tLOutput
 **
 **  (see tOutput.h for a description of these classes)
 **
 **    Modifications:
 **     - 6/01 GT added optional output of channel widths to file
 **       *.chanwid. Activated if Parker-Paola width model used.
 **       If so, channel depths are also output.
 **     - 4/03 AD added canonical output
 **     - 7/03 AD added tOutputBase and tTSOutputImp
 **     - 8/03: AD Random number generator handling
 **
 **  $Id: tOutput.cpp,v 1.85 2003-08-06 12:38:37 childcvs Exp $
 */
/*************************************************************************/

#include <math.h>    // For fmod function
#include <stdlib.h> // For qsort
#include <string.h>
#include <assert.h>

#if !defined(HAVE_NO_NAMESPACE)
# include <iostream>
using namespace std;
#else
# include <iostream.h>
#endif
#include "tOutput.h"
#include "../errors/errors.h"
#include "../tMeshList/tMeshList.h"
#include "../tStreamNet/tStreamNet.h" // For k2DKinematicWave and kHydrographPeakMethod


/**************************************************************************/
/**
 ** @class tTSOutputImp
 **
 ** Handles application-specific time-series data for the CHILD model.
 **
 ** Modifications:
 ** - 07/03 moved from tLOutput (AD)
 **
 */
/**************************************************************************/
template< class tSubNode >
class tTSOutputImp : public tOutputBase<tSubNode>
{
  tTSOutputImp(const tTSOutputImp&);
  tTSOutputImp& operator=(const tTSOutputImp&);
public:
  tTSOutputImp( tMesh<tSubNode> * meshPtr, tInputFile &infile );
  void WriteTSOutput();
private:
  ofstream volsofs;    // catchment volume
  ofstream dvolsofs;
  ofstream tareaofs;   // total voronoi area of catchment
  ofstream vegcovofs;  // Catchment vegetation cover %
  double mdLastVolume;
};

/*************************************************************************\
 **
 **  Constructor
 **
 **  The constructor takes two arguments, a pointer to the mesh and
 **  a reference to an open input file. It reads the base name for the
 **  output files from the input file.
 **
 **  Input: meshPtr -- pointer to a tMesh object (or descendant), assumed
 **                    valid
 **         infile -- reference to an open input file, assumed valid
 **
\*************************************************************************/
template< class tSubNode >
tOutputBase<tSubNode>::tOutputBase( tMesh<tSubNode> * meshPtr, tInputFile &infile ) :
  m(meshPtr)
{
  assert( meshPtr != 0 );

  infile.ReadItem( baseName, sizeof(baseName), "OUTFILENAME" );
}

/*************************************************************************\
 **
 **  tOutputBase::CreateAndOpenFile
 **
 **  Opens the output file stream pointed to by theOFStream, giving it the
 **  name <baseName><extension>, and checks to make sure that the ofstream
 **  is valid.
 **
 **  Input:  theOFStream -- ptr to an ofstream object
 **          extension -- file name extension (e.g., ".nodes")
 **  Output: theOFStream is initialized to create an open output file
 **  Assumes: extension is a null-terminated string, and the length of
 **           baseName plus extension doesn't exceed kMaxNameSize+6
 **           (ie, the extension is expected to be <= 6 characters)
 **
\*************************************************************************/
template< class tSubNode >
void tOutputBase<tSubNode>::CreateAndOpenFile( ofstream *theOFStream,
					       const char *extension ) const
{
  char fullName[kMaxNameSize+6];  // name of file to be created

  // workaround for an obscure bug re. gcc/valgrind
#ifndef __GNUC__
  assert(strlen(baseName)+strlen(extension) < sizeof(fullName));
#endif

  strcpy( fullName, baseName );
  strcat( fullName, extension );
  theOFStream->open( fullName );

  if( !theOFStream->good() )
    ReportFatalError(
		     "I can't create files for output. Storage space may be exhausted.");
  theOFStream->precision( 12 );
}

/*************************************************************************\
 **
 **  tBaseOutput::WriteTimeNumberElements
 **
 **  write time and number of elements
\*************************************************************************/
template< class tSubNode >
void tOutputBase<tSubNode>::WriteTimeNumberElements( ofstream &fs,
						     double time, int n )
{
  fs << ' ' << time << '\n' << n << '\n';
}

/*************************************************************************\
 **
 **  Constructor
 **
 **  The constructor takes two arguments, a pointer to the mesh and
 **  a reference to an open input file. It opens and initializes the
 **  output files.
 **
 **  Input: meshPtr -- pointer to a tMesh object (or descendant), assumed
 **                    valid
 **         infile -- reference to an open input file, assumed valid
 **
\*************************************************************************/
template< class tSubNode >
tOutput<tSubNode>::tOutput( tMesh<tSubNode> * meshPtr, tInputFile &infile ) :
  tOutputBase<tSubNode>( meshPtr, infile ),  // call base-class constructor
  CanonicalNumbering(true)
{
  CreateAndOpenFile( &nodeofs, SNODES );
  CreateAndOpenFile( &edgofs, SEDGES );
  CreateAndOpenFile( &triofs, STRI );
  CreateAndOpenFile( &zofs, SZ );
  CreateAndOpenFile( &vaofs, SVAREA );
}


/*************************************************************************\
 **
 **  tOutput::WriteOutput
 **
 **  This function writes information about the mesh to four files called
 **  name.nodes, name.edges, name.tri, and name.z, where "name" is a
 **  name that the user has specified in the input file and which is
 **  stored in the data member baseName.
 **
 **  Input: time -- time of the current output time-slice
 **  Output: the node, edge, and triangle ID numbers are modified so that
 **          they are numbered according to their position on the list
 **  Assumes: the four file ofstreams have been opened by the constructor
 **           and are valid
 **
 **  TODO: deal with option for once-only printing of mesh when mesh not
 **        deforming
\*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::WriteOutput( double time )
{
  tMeshListIter<tSubNode> niter( m->getNodeList() ); // node list iterator
  tMeshListIter<tEdge> eiter( m->getEdgeList() );    // edge list iterator
  tListIter<tTriangle> titer( m->getTriList() );     // tri list iterator
  const int nnodes = m->getNodeList()->getSize();  // # of nodes on list
  const int nedges = m->getEdgeList()->getSize();  // "    edges "
  const int ntri = m->getTriList()->getSize();     // "    triangles "

  if (1)//DEBUG
    cout << "tOutput::WriteOutput() time=" << time << endl;

  // Renumber IDs in order by position on list
  if (!CanonicalNumbering)
    RenumberIDInListOrder();
  else
    RenumberIDCanonically();

  // Write node file, z file, and varea file
  WriteTimeNumberElements( nodeofs, time, nnodes);
  WriteTimeNumberElements( zofs, time, nnodes);
  WriteTimeNumberElements( vaofs, time, nnodes);
  if (!CanonicalNumbering) {
    for( tNode *cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() )
      WriteNodeRecord( cn );
  } else {
    // write nodes in ID order
    tIdArray< tSubNode > RNode(*(m->getNodeList()));
    for( int i=0; i<nnodes; ++i )
      WriteNodeRecord( RNode[i] );
  }

  // Write edge file
  WriteTimeNumberElements( edgofs, time, nedges);
  if (!CanonicalNumbering) {
    for( tEdge *ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() )
      WriteEdgeRecord( ce );
  } else {
    // write edges in ID order
    tIdArray< tEdge > REdge(*(m->getEdgeList()));
    for( int i=0; i<nedges; ++i )
      WriteEdgeRecord( REdge[i] );
  }

  // Write triangle file
  WriteTimeNumberElements( triofs, time, ntri);
  if (!CanonicalNumbering) {
    for( tTriangle *ct=titer.FirstP(); !(titer.AtEnd()); ct=titer.NextP() )
      WriteTriangleRecord( ct );
  } else {
    // write triangles in ID order
    tIdArray< tTriangle > RTri(*(m->getTriList()));
    for( int i=0; i<ntri; ++i ) {
      assert( RTri[i]->isIndexIDOrdered() );
      WriteTriangleRecord( RTri[i] );
    }
  }

  nodeofs << flush;
  zofs << flush;
  vaofs << flush;
  edgofs << flush;
  triofs << flush;

  // Call virtual function to write any additional data
  WriteNodeData( time );

  if (0)//DEBUG
    cout << "tOutput::WriteOutput() Output done" << endl;
}

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

/*************************************************************************\
 **
 **  tOutput::RenumberID
 **
 **  Set IDs in list order
 **
 **  AD, April 2003
\*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::RenumberIDInListOrder()
{
  m->ResetNodeID();
  m->ResetEdgeID();
  m->ResetTriangleID();
}

/*************************************************************************\
 **
 **  tOutput::RenumberIDCanonically
 **
 **  Set IDs in a canonical ordering independent of the list ordering
 **  As well, set tNode.edg to the spoke with the lowest destination node ID
 **
 **  AD, April-May 2003
\*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::RenumberIDCanonically()
{
  tMeshListIter<tSubNode> niter( m->getNodeList() ); // node list iterator
  tMeshListIter<tEdge> eiter( m->getEdgeList() );    // edge list iterator
  tListIter<tTriangle> titer( m->getTriList() );     // tri list iterator
  const int nnodes = m->getNodeList()->getSize();  // # of nodes on list
  const int nedges = m->getEdgeList()->getSize();  // "    edges "
  const int ntri = m->getTriList()->getSize();     // "    triangles "

  {
    // First we set the Nodes Id in the order defined below
    // b1 <= b2 then x1 <= x2 then y1<=y2
    tArray< tNode* > RNode(nnodes);
    int i;
    tNode *cn;
    for( cn=niter.FirstP(), i=0; i<nnodes; cn=niter.NextP(), ++i )
      RNode[i] = cn;

    qsort(RNode.getArrayPtr(), RNode.getSize(), sizeof(RNode[0]),
	  orderRNode
	  );
    const int s_ = RNode.getSize();
    for(i=0; i<s_; ++i)
      RNode[i]->setID(i);
    m->SetmiNextNodeID( RNode.getSize() );
  }
  {
    // Set tNode.edg to the spoke that links to the destination node with the
    // lowest ID
    tNode *cn;
    for( cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() ) {
      tSpkIter sI( cn );
      tEdge *thece = cn->getEdg();
      tEdge *ce;
      for( ce = sI.FirstP(); !( sI.AtEnd() ); ce = sI.NextP() ) {
	if (ce->getDestinationPtr()->getID() <
	    thece->getDestinationPtr()->getID())
	  thece = ce;
      }
      cn->setEdg(thece);
    }
  }
  // Set Edge ID with respect to Node ID
  // #1 The pair edge-complement is ordered so that
  // (IDorig IDdest) (IDdest IDorig) with IDorig < IDdest
  {
    tEdge *ce;
    for( ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() ) {
      tListNode< tEdge >* enodePtr1 = eiter.NodePtr();
      eiter.Next();
      tListNode< tEdge >* enodePtr2 = eiter.NodePtr();
      if (ce->getOriginPtr()->getID() > ce->getDestinationPtr()->getID()) {
	// swap edges
	m->getEdgeList()->moveToAfter(enodePtr1, enodePtr2);
	eiter.Next();
      }
    }
  }
  // #2 Then pairs are ordered with IDorig1 < IDorig2 and
  // if IDorig1 == IDorig2 IDdest1 < IDdest2
  {
    tArray< tEdge* > REdge2(nedges/2);
    tEdge *ce;
    int i;
    for( ce=eiter.FirstP(), i=0; !(eiter.AtEnd()); ce=eiter.NextP(), ++i ) {
      REdge2[i] = ce;
      eiter.Next();
    }
    qsort(REdge2.getArrayPtr(), REdge2.getSize(), sizeof(REdge2[0]),
	  orderREdge
	  );
    const int s_ = REdge2.getSize();
    for(i=0; i<s_; ++i) {
      assert (REdge2[i]->getOriginPtr()->getID() <
	      REdge2[i]->getDestinationPtr()->getID() );
      REdge2[i]->setID(2*i);
      REdge2[i]->getComplementEdge()->setID(2*i+1);
    }
    m->SetmiNextEdgID( 2*REdge2.getSize() );
  }
  {
    // Set Triangle Id so that the vertexes are ordered
    tArray< tTriangle* > RTri(ntri);
    int i;
    tTriangle *ct;
    for( ct=titer.FirstP(), i=0; i<ntri; ct=titer.NextP(), ++i ){
      ct->SetIndexIDOrdered(); // set ct->index_ in ID order
      RTri[i] = ct;
    }
    qsort(RTri.getArrayPtr(), RTri.getSize(), sizeof(RTri[0]),
	  orderRTriangle
	  );
    const int s_ = RTri.getSize();
    for(i=0; i<s_; ++i)
      RTri[i]->setID(i);
    m->SetmiNextTriID( RTri.getSize() );
  }
}

// qsort comparison function for canonical nodes ordering
template< class tSubNode >
int tOutput<tSubNode>::orderRNode( const void *a_, const void *b_ )
{
  const tNode *N1 = *static_cast<tNode const *const *>(a_);
  const tNode *N2 = *static_cast<tNode const *const *>(b_);

  const int N1B = N1->getBoundaryFlag();
  const int N2B = N2->getBoundaryFlag();
  if (N1B < N2B)
    return -1;
  if (N1B > N2B)
    return 1;

  const double N1X = N1->getX();
  const double N2X = N2->getX();
  if (N1X < N2X)
    return -1;
  if (N1X > N2X)
    return 1;
  const double N1Y = N1->getY();
  const double N2Y = N2->getY();
  if (N1Y < N2Y)
    return -1;
  if (N1Y > N2Y)
    return 1;
  return 0;
}

// qsort comparison function for canonical edges ordering
template< class tSubNode >
int tOutput<tSubNode>::orderREdge( const void *a_, const void *b_ )
{
  const tEdge *E1 = *static_cast<tEdge const *const *>(a_);
  const tEdge *E2 = *static_cast<tEdge const *const *>(b_);

  const int o1 = E1->getOriginPtr()->getID();
  const int o2 = E2->getOriginPtr()->getID();
  if (o1 < o2) return -1;
  if (o1 > o2) return 1;
  const int d1 = E1->getDestinationPtr()->getID();
  const int d2 = E2->getDestinationPtr()->getID();
  if (d1 < d2) return -1;
  if (d1 > d2) return 1;
  assert(0); /*NOTREACHED*/
  abort();
}

// qsort comparison function for canonical triangles ordering
template< class tSubNode >
int tOutput<tSubNode>::orderRTriangle( const void *a_, const void *b_ )
{
  tTriangle *T1 = *static_cast<tTriangle *const *>(a_);
  tTriangle *T2 = *static_cast<tTriangle *const *>(b_);

  assert(T1 && T2);
  // The comparison itself related to vertexes ID
  for(int i=0;i<3;++i) {
    const int i1 = T1->pPtr(T1->index()[i])->getID();
    const int i2 = T2->pPtr(T2->index()[i])->getID();
    if (i1 < i2) return -1;
    if (i1 > i2) return 1;
  }
  assert(0); /*NOTREACHED*/
  abort();
}

/*************************************************************************\
 **
 **  tOutput::WriteNodeData
 **
 **  This is a virtual function which can be overridden to write any
 **  additional node data. The base class version does nothing.
 **
\*************************************************************************/
template< class tSubNode >
void tOutput<tSubNode>::WriteNodeData( double /* time */ )
{}


/*************************************************************************\
 **
 **  tLOutput constructor
 **
 **  Creates and opens a series of files for drainage areas, slopes, etc.
 **
 **  Modifications:
 **    - 1/00 added "opOpt" and creation of veg output file (GT)
 **    - added flow depth output file (GT 1/00)
 **    - added
\*************************************************************************/
template< class tSubNode >
tLOutput<tSubNode>::tLOutput( tMesh<tSubNode> *meshPtr, tInputFile &infile,
			      tRand *rand_) :
  tOutput<tSubNode>( meshPtr, infile ),  // call base-class constructor
  TSOutput(0),
  rand(rand_),
  counter(0)
{
  int opOpt;  // Optional modules: only output stuff when needed
  CreateAndOpenFile( &randomofs, SRANDOM );
  CreateAndOpenFile( &drareaofs, ".area" );
  CreateAndOpenFile( &netofs, ".net" );
  CreateAndOpenFile( &slpofs, ".slp" );
  CreateAndOpenFile( &qofs, ".q" );
  CreateAndOpenFile( &texofs, ".tx" );
  CreateAndOpenFile( &tauofs, ".tau" );

  // Vegetation cover: if dynamic vegetation option selected
  if( (opOpt = infile.ReadItem( opOpt, "OPTVEG" ) ) != 0)
    CreateAndOpenFile( &vegofs, ".veg" );

  // Flow depth: if kinematic wave option used OR if channel geometry
  // model other than "regime" used
  if( ((opOpt = infile.ReadItem( opOpt, "FLOWGEN" )) == k2DKinematicWave )
      || (opOpt = infile.ReadItem( opOpt, "CHAN_GEOM_MODEL"))>1 )
    CreateAndOpenFile( &flowdepofs, ".dep" );

  // Time-series output: if requested
  if( (opOpt = infile.ReadItem( opOpt, "OPTTSOUTPUT" ) ) != 0) {
    TSOutput = new tTSOutputImp< tSubNode >(meshPtr, infile);
  }

  // Channel width output: if the channel geometry model is other
  // than 1 (code for empirical regime channels)
  if( (opOpt = infile.ReadItem( opOpt, "CHAN_GEOM_MODEL" ) ) > 1 )
    CreateAndOpenFile( &chanwidthofs, ".chanwid" );

  // Flow path length output: if using hydrograph peak method for
  // computing discharge
  if( (opOpt = infile.ReadItem( opOpt, "FLOWGEN" )) == kHydrographPeakMethod )
    CreateAndOpenFile( &flowpathlenofs, ".fplen" );

  // Sediment flux: if not using detachment-limited option
  if( (opOpt = infile.ReadItem( opOpt, "OPTDETACHLIM" ) ) == 0)
    CreateAndOpenFile( &qsofs, ".qs" );
}

/*************************************************************************\
 **
 **  tLOutput constructor
 **
 **  Creates and opens a series of files for drainage areas, slopes, etc.
 **
 **  Modifications:
 **    - 1/00 added "opOpt" and creation of veg output file (GT)
 **    - added flow depth output file (GT 1/00)
 **    - added
\*************************************************************************/
template< class tSubNode >
tLOutput<tSubNode>::~tLOutput()
{
  delete TSOutput;
}

/*************************************************************************\
 **
 **  tLOutput::WriteNodeData
 **
 **  This overridden virtual function writes output for tLNodes, including
 **  drainage areas, flow pathways, slopes, discharges, layer info, etc.
 **
 **  Modifications:
 **    - 1/00 added output to veg output file (GT)
 **    - added output of flow depth; made slope output for all nodes (GT 1/00)
 **    - 6/00 layer info for each time step written to a different file (NG)
 **    - 9/01 added output of flow path length (GT)
\*************************************************************************/
//TODO: should output boundary points as well so they'll map up with nodes
// for plotting. Means changing getSlope so it returns zero if flowedg
// undefined
template< class tSubNode >
void tLOutput<tSubNode>::WriteNodeData( double time )
{
  //for writing out layer info to different files at each time
  const char* const nums("0123456789");

  tMeshListIter<tSubNode> ni( m->getNodeList() ); // node list iterator
  const int nActiveNodes = m->getNodeList()->getActiveSize(); // # active nodes
  const int nnodes = m->getNodeList()->getSize(); // total # nodes

  //taking care of layer file, since new one each time step
  char ext[7];
  strcpy( ext, ".lay");
  if(counter<10)
    strncat( ext, &nums[counter], 1);
  else if(counter>=10){
    strncat(ext, &nums[counter/10], 1);
    strncat(ext, &nums[static_cast<int>( fmod(static_cast<double>(counter),10.0) )], 1);
  }
  CreateAndOpenFile( &layofs, ext );
  counter++;

  // Write current time in each file
  WriteTimeNumberElements( randomofs, time, rand->numberRecords());
  WriteTimeNumberElements( drareaofs, time, nActiveNodes);
  WriteTimeNumberElements( netofs, time, nActiveNodes);
  WriteTimeNumberElements( slpofs, time, nnodes);
  WriteTimeNumberElements( qofs, time, nnodes);
  WriteTimeNumberElements( layofs, time, nActiveNodes);
  WriteTimeNumberElements( texofs, time, nnodes);
  WriteTimeNumberElements( tauofs, time, nnodes);
  if( vegofs.good() )
    WriteTimeNumberElements( vegofs, time, nnodes);
  if( flowdepofs.good() )
    WriteTimeNumberElements( flowdepofs, time, nnodes);
  if( chanwidthofs.good() )
    WriteTimeNumberElements( chanwidthofs, time, nnodes);
  if( flowpathlenofs.good() )
    WriteTimeNumberElements( flowpathlenofs, time, nnodes);
  if( qsofs.good() )
    WriteTimeNumberElements( qsofs, time, nnodes);

  // Write Random number generator state
  rand->dumpToFile( randomofs );
  // Write data
  if (!CanonicalNumbering) {
    tSubNode *cn;   // current node
    for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      WriteActiveNodeData( cn );
    for( cn = ni.FirstP(); !(ni.AtEnd()); cn = ni.NextP() )
      WriteAllNodeData( cn );
  } else {
    // write in node ID order
    tIdArray< tSubNode > RNode(*(m->getNodeList()));
    int i;
    for( i=0; i<nActiveNodes; ++i )
      WriteActiveNodeData( RNode[i] );
    for( i=0; i<nnodes; ++i )
      WriteAllNodeData( RNode[i] );
  }

  randomofs << flush;
  drareaofs << flush;
  netofs << flush;
  slpofs << flush;
  qofs << flush;
  texofs << flush;
  tauofs << flush;
  if( vegofs.good() ) vegofs << flush;
  if( flowdepofs.good() ) flowdepofs << flush;
  if( chanwidthofs.good() ) chanwidthofs << flush;
  if( flowpathlenofs.good() ) flowpathlenofs << flush;
  if( qsofs.good() ) qsofs << flush;

  layofs.close();
}

// Write data, including layer info
template< class tSubNode >
inline void tLOutput<tSubNode>::WriteActiveNodeData( tSubNode *cn )
{
  int i, j;      // counters

  assert( cn!=0 );
  drareaofs << cn->getDrArea() << '\n';
  if( cn->getDownstrmNbr() )
    netofs << cn->getDownstrmNbr()->getID() << '\n';
  layofs << ' ' << cn->getNumLayer() << '\n';
  i=0;
  while(i<cn->getNumLayer()){
    layofs << cn->getLayerCtime(i) << ' ' << cn->getLayerRtime(i) << ' '
	   << cn->getLayerEtime(i) << '\n'
	   << cn->getLayerDepth(i) << ' ' << cn->getLayerErody(i) << ' '
	   << cn->getLayerSed(i) << '\n';
    j=0;
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
  slpofs << (!cn->getBoundaryFlag() ? cn->getSlope():0.) << '\n';
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

/*************************************************************************\
 **
 **  tLOutput::WriteTSOutput
 **  This function writes the total volume of the DEM above the datum to
 **  a file called name.vols, where "name" is a name that the user has
 **  specified in the input file and which is stored in the data member
 **  baseName.
 **
\*************************************************************************/
template< class tSubNode >
void tLOutput<tSubNode>::WriteTSOutput()
{
  if (TSOutput) TSOutput->WriteTSOutput();
}


template< class tSubNode >
bool tLOutput<tSubNode>::OptTSOutput() const { return BOOL(TSOutput!=0); }

/*************************************************************************\
 **
 **  tTSOutputImp constructor
 **
 **  Creates and opens a series of files for time series
 **
 **  Modifications:
\*************************************************************************/
template< class tSubNode >
tTSOutputImp<tSubNode>::tTSOutputImp( tMesh<tSubNode> *meshPtr, tInputFile &infile ) :
  tOutputBase<tSubNode>( meshPtr, infile ),  // call base-class constructor
  mdLastVolume(0.)
{
  int opOpt;  // Optional modules: only output stuff when needed

  CreateAndOpenFile( &volsofs, ".vols" );
  CreateAndOpenFile( &dvolsofs, ".dvols" );
  if( (opOpt = infile.ReadItem( opOpt, "OPTVEG" ) ) != 0)
    CreateAndOpenFile( &vegcovofs, ".vcov" );
  CreateAndOpenFile( &tareaofs, ".tarea" );
}

/*************************************************************************\
 **
 **  tTSOutputImp::WriteTSOutput
 **  This function writes the total volume of the DEM above the datum to
 **  a file called name.vols, where "name" is a name that the user has
 **  specified in the input file and which is stored in the data member
 **  baseName.
 **
\*************************************************************************/
template< class tSubNode >
void tTSOutputImp<tSubNode>::WriteTSOutput()
{
  tMeshListIter<tSubNode> niter( m->getNodeList() ); // node list iterator

  tSubNode * cn;       // current node

  double volume = 0.,
    area = 0.,
    cover = 0.;

  if (0)//DEBUG
    cout << "tTSOutputImp::WriteTSOutput()" << endl;

  for( cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() ) {
    volume += cn->getZ()*cn->getVArea();
    area += cn->getVArea();
  }

  volsofs << volume << endl;
  if( mdLastVolume > 0.0 )
    dvolsofs << volume - mdLastVolume << endl;
  mdLastVolume = volume;
  //tareaofs << area << endl;

  if( vegcovofs.good() ) {
    for( cn = niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() )
      cover += cn->getVegCover().getVeg()*cn->getVArea();
    vegcovofs << cover/area << endl;
  }
}
