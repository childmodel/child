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
 **     - 6/03 QC added output for 20 stratigraphic section positions
 **       this now result in a large number of output files created
 **     - 7/03 AD added tOutputBase and tTSOutputImp
 **     - 8/03: AD Random number generator handling
 **
 **  $Id: tOutput.cpp,v 1.105 2008/07/07 16:18:58 childcvs Exp $
 */
/*************************************************************************/

#include <math.h>    // For fmod function
#include <stdlib.h> // For qsort
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include <iostream>
#include "tOutput.h"
#include "../errors/errors.h"
#include "../tMeshList/tMeshList.h"
#include "../tStreamNet/tStreamNet.h" // For k2DKinematicWave and kHydrographPeakMethod
#include "../tStratGrid/tStratGrid.h"
#include "../tFloodplain/tFloodplain.h"


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
  tTSOutputImp();
public:
  tTSOutputImp( tMesh<tSubNode> * meshPtr, const tInputFile &infile );
  void WriteTSOutput();
private:
  std::ofstream volsofs;    // catchment volume
  std::ofstream dvolsofs;
  std::ofstream tareaofs;   // total voronoi area of catchment
  std::ofstream vegcovofs;  // Catchment vegetation cover %
  double mdLastVolume;
};

/**************************************************************************/
/**
** @class tStratOutputImp
**
** Class tLOutput handles application-specific data for the CHILD model.
** The constructor creates additional output files, and the overloaded
** WriteNodeData function writes the data to files.
**
** Modifications:
**  - 6/03 added 3 x 10 files for stratigraphic sections (z,texture,facies)
**    (QC)
*/
/**************************************************************************/
template< class tSubNode >
class tStratOutputImp : public tOutputBase<tSubNode>
{
  tStratOutputImp();
public:
  tStratOutputImp(tMesh<tSubNode> * meshPtr, const tInputFile &infile );
  void WriteNodeData( double time, int );
  void SetStratGrid(tStratGrid *, tStreamNet *);

private:
  typedef enum{ DIRI, DIRJ } direction_t;
  void WriteStratGridSections( double time );
  void WriteGravelBodies( double time, int );
  void WriteCompleteStratigraphy( double time, int );
  void WritePreservationPotential( double time, int );
  void WriteSingleSection( double time,
			   int section,
			   std::ofstream &xyzofs, std::ofstream &lay3ofs, std::ofstream &lay4ofs,
			   direction_t );
  void Flush();

  // stratigraphic sections 1 to 10
  enum{ nIsections = 5, nSections = 10 };
  std::ofstream xyzofs[nSections];
  std::ofstream layt3ofs[nSections];
  std::ofstream layt4ofs[nSections];
  direction_t direction[nSections];

  // Subsurface gravel bodies
  std::ofstream grxyzofs,gravofs;

  // Preservation potential files
  std::ofstream psurfofs,pssurfofs,pssurf2ofs;

  tStratGrid *stratGrid; // pointer to stratigraphy grid
  tStreamNet *netPtr;
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
tOutputBase<tSubNode>::tOutputBase( tMesh<tSubNode> * meshPtr,
				    const tInputFile &infile ) :
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
void tOutputBase<tSubNode>::CreateAndOpenFile( std::ofstream *theOFStream,
					       const char *extension ) const
{
  char fullName[kMaxNameSize+20];  // name of file to be created

  if (strlen(baseName)+strlen(extension) >= sizeof(fullName)) {
    std::cout << "While opening " << baseName << extension << std::endl;
    ReportFatalError("tOutputCreateAndOpenFile(): buffer too short.");
  }

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
void tOutputBase<tSubNode>::WriteTimeNumberElements( std::ofstream &fs,
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
tOutput<tSubNode>::tOutput( tMesh<tSubNode> * meshPtr,
			    const tInputFile &infile ) :
  tOutputBase<tSubNode>( meshPtr, infile ),  // call base-class constructor
  CanonicalNumbering(true)
{
  this->CreateAndOpenFile( &nodeofs, SNODES );
  this->CreateAndOpenFile( &edgofs, SEDGES );
  this->CreateAndOpenFile( &triofs, STRI );
  this->CreateAndOpenFile( &zofs, SZ );
  this->CreateAndOpenFile( &vaofs, SVAREA );
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
  typename tMesh< tSubNode >::nodeListIter_t niter( this->m->getNodeList() ); // node list iterator
  typename tMesh< tSubNode >::edgeListIter_t eiter( this->m->getEdgeList() ); // edge list iterator
  typename tMesh< tSubNode >::triListIter_t titer( this->m->getTriList() );   // tri list iterator
  const int nnodes = this->m->getNodeList()->getSize();  // # of nodes on list
  const int nedges = this->m->getEdgeList()->getSize();  // "    edges "
  const int ntri = this->m->getTriList()->getSize();     // "    triangles "

  if (1)//DEBUG
    std::cout << "tOutput::WriteOutput() time=" << time << std::endl;

  // Renumber IDs in order by position on list
  if (!CanonicalNumbering)
    RenumberIDInListOrder();
  else
    this->m->RenumberIDCanonically();

  // Write node file, z file, and varea file
  this->WriteTimeNumberElements( nodeofs, time, nnodes);
  this->WriteTimeNumberElements( zofs, time, nnodes);
  this->WriteTimeNumberElements( vaofs, time, nnodes);
  if (!CanonicalNumbering) {
    for( tNode *cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() )
      WriteNodeRecord( cn );
  } else {
    // write nodes in ID order
    typename tMesh< tSubNode >::tIdArrayNode_t  RNode(*(this->m->getNodeList()));
    for( int i=0; i<nnodes; ++i )
      WriteNodeRecord( RNode[i] );
  }

  // Write edge file
  this->WriteTimeNumberElements( edgofs, time, nedges);
  if (!CanonicalNumbering) {
    for( tEdge *ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() )
      WriteEdgeRecord( ce );
  } else {
    // write edges in ID order
    typename tMesh< tSubNode >::tIdArrayEdge_t REdge(*(this->m->getEdgeList()));
    for( int i=0; i<nedges; ++i )
      WriteEdgeRecord( REdge[i] );
  }

  // Write triangle file
  this->WriteTimeNumberElements( triofs, time, ntri);
  if (!CanonicalNumbering) {
    for( tTriangle *ct=titer.FirstP(); !(titer.AtEnd()); ct=titer.NextP() )
      WriteTriangleRecord( ct );
  } else {
    // write triangles in ID order
    typename tMesh< tSubNode >::tIdArrayTri_t RTri(*(this->m->getTriList()));
    for( int i=0; i<ntri; ++i ) {
      assert( RTri[i]->isIndexIDOrdered() );
      WriteTriangleRecord( RTri[i] );
    }
  }

  nodeofs << std::flush;
  zofs << std::flush;
  vaofs << std::flush;
  edgofs << std::flush;
  triofs << std::flush;

  // Call virtual function to write any additional data
  WriteNodeData( time );

  if (0)//DEBUG
    std::cout << "tOutput::WriteOutput() Output done" << std::endl;
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
  this->m->ResetNodeID();
  this->m->ResetEdgeID();
  this->m->ResetTriangleID();
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
tLOutput<tSubNode>::tLOutput( tMesh<tSubNode> *meshPtr,
			      const tInputFile &infile, tRand *rand_) :
  tOutput<tSubNode>( meshPtr, infile ),  // call base-class constructor
  TSOutput(0),
  stratOutput(0),
  rand(rand_),
  counter(0),
  Surfer(false)
  
{
  int opOpt;  // Optional modules: only output stuff when needed
  this->CreateAndOpenFile( &randomofs, SRANDOM );
  this->CreateAndOpenFile( &drareaofs, ".area" );
  this->CreateAndOpenFile( &netofs, ".net" );
  this->CreateAndOpenFile( &slpofs, ".slp" );
  this->CreateAndOpenFile( &qofs, ".q" );
  this->CreateAndOpenFile( &texofs, ".tx" );
  this->CreateAndOpenFile( &tauofs, ".tau" );
  this->CreateAndOpenFile( &permIDofs, ".id" );

  //Layer output: only write layer information if user selects to write it
  OptLayOutput = infile.ReadBool( "OPTLAYEROUTPUT" );

  // Vegetation cover: if dynamic vegetation option selected
  if( (opOpt = infile.ReadItem( opOpt, "OPTVEG" ) ) != 0)
    this->CreateAndOpenFile( &vegofs, SVEG );

  // Flow depth: if kinematic wave option used OR if channel geometry
  // model other than "regime" used
  if( (static_cast<tStreamNet::kFlowGen_t>(opOpt = infile.ReadItem( opOpt, "FLOWGEN" ))
       == tStreamNet::k2DKinematicWave )
      || (opOpt = infile.ReadItem( opOpt, "CHAN_GEOM_MODEL"))>1 )
    this->CreateAndOpenFile( &flowdepofs, ".dep" );

  // Time-series output: if requested
  if( infile.ReadBool( "OPTTSOUTPUT" ) ) {
    TSOutput = new tTSOutputImp< tSubNode >(meshPtr, infile);
  }

  // Channel width output: if the channel geometry model is other
  // than 1 (code for empirical regime channels)
  if( (opOpt = infile.ReadItem( opOpt, "CHAN_GEOM_MODEL" ) ) > 1 )
    this->CreateAndOpenFile( &chanwidthofs, ".chanwid" );

  // Flow path length output: if using hydrograph peak method for
  // computing discharge
  if( static_cast<tStreamNet::kFlowGen_t>(opOpt = infile.ReadItem( opOpt, "FLOWGEN" ))
      == tStreamNet::kHydrographPeakMethod )
    this->CreateAndOpenFile( &flowpathlenofs, ".fplen" );

  // Sediment flux: if not using detachment-limited option
  if( (opOpt = infile.ReadItem( opOpt, "OPTDETACHLIM" ) ) == 0){
     this->CreateAndOpenFile( &qsofs, ".qs" );
     this->CreateAndOpenFile( &qsinofs, ".qsin" );
     this->CreateAndOpenFile( &qsdinofs, ".qsdin" );
     this->CreateAndOpenFile( &dzdtofs, ".dzdt" );
  }  

  this->CreateAndOpenFile( &upofs, ".up" );

  // If Rectangular Stratigraphy Grid, open several files
  // for writing the stratigraphy at fixed positions
  int optStratGrid;
  if( (optStratGrid = infile.ReadItem(optStratGrid,"OPTSTRATGRID")) !=0){
    stratOutput = new tStratOutputImp< tSubNode >(meshPtr, infile);
  }
  {
    int opOpt;
    Surfer = (opOpt = infile.ReadItem( opOpt, "SURFER", false)) != 0;
  }
  //XSurfer = infile.ReadBool( "SURFER", false);
}

/*************************************************************************\
 **
 **  tLOutput destructor
 **
 **  QC/AD 7/2003
\*************************************************************************/
template< class tSubNode >
tLOutput<tSubNode>::~tLOutput()
{
  delete TSOutput;
  delete stratOutput;
  rand = 0;
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
 **    - 5/03 added output in simple x,y,z style for visualisation in Surfer(QC)
 **    - 6/03 added call to function writing stratigraphic sections (QC)
\*************************************************************************/
//TODO: should output boundary points as well so they'll map up with nodes
// for plotting. Means changing getSlope so it returns zero if flowedg
// undefined
template< class tSubNode >
void tLOutput<tSubNode>::WriteNodeData( double time )
{
  //for writing out layer info to different files at each time
  const char* const nums("0123456789");

  typename tMesh< tSubNode >::nodeListIter_t ni( this->m->getNodeList() ); // node list iterator
  const int nActiveNodes = this->m->getNodeList()->getActiveSize(); // # active nodes
  const int nnodes = this->m->getNodeList()->getSize(); // total # nodes

  if(OptLayOutput){
    //taking care of layer and x,y,z file, since new one each time step
    char ext[7];
    strcpy( ext, ".lay");
    if(counter<10)
      strncat( ext, &nums[counter], 1);
    else {
      strncat(ext, &nums[counter/10], 1);
      strncat(ext, &nums[static_cast<int>( fmod(static_cast<double>(counter),10.0) )], 1);
    }
    this->CreateAndOpenFile( &layofs, ext );
  }

#define MY_EXT ".surf"
  char extt[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( extt, "%s%d", MY_EXT, counter );
#undef MY_EXT

  if(Surfer)
    this->CreateAndOpenFile( &surfofs, extt );

  // *Counter that counts the number of write timesteps* 
  counter++;

  // Write current time in each file
  this->WriteTimeNumberElements( randomofs, time, rand->numberRecords());
  this->WriteTimeNumberElements( drareaofs, time, nActiveNodes);
  this->WriteTimeNumberElements( netofs, time, nActiveNodes);
  this->WriteTimeNumberElements( slpofs, time, nnodes);
  this->WriteTimeNumberElements( qofs, time, nnodes);
  if(OptLayOutput)
    this->WriteTimeNumberElements( layofs, time, nActiveNodes);
  this->WriteTimeNumberElements( texofs, time, nnodes);
  if( surfofs.good() )
    this->WriteTimeNumberElements( surfofs, time, nActiveNodes);

  this->WriteTimeNumberElements( tauofs, time, nnodes);
  if( vegofs.good() )
    this->WriteTimeNumberElements( vegofs, time, nnodes);
  if( flowdepofs.good() )
    this->WriteTimeNumberElements( flowdepofs, time, nnodes);
  if( chanwidthofs.good() )
    this->WriteTimeNumberElements( chanwidthofs, time, nnodes);
  if( flowpathlenofs.good() )
    this->WriteTimeNumberElements( flowpathlenofs, time, nnodes);
  if( qsofs.good() )
    this->WriteTimeNumberElements( qsofs, time, nnodes);
  if( qsinofs.good() )
    this->WriteTimeNumberElements( qsinofs, time, nnodes);
  if( qsdinofs.good() )
    this->WriteTimeNumberElements( qsdinofs, time, nnodes);
  if( dzdtofs.good() )
    this->WriteTimeNumberElements( dzdtofs, time, nnodes);
  if( upofs.good() )
    this->WriteTimeNumberElements( upofs, time, nnodes);
  if( permIDofs.good() )
    this->WriteTimeNumberElements( permIDofs, time, nnodes );

  // Write Random number generator state
  rand->dumpToFile( randomofs );
  // Write data
  if (!this->CanonicalNumbering) {
    tSubNode *cn;   // current node
    for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
      WriteActiveNodeData( cn );
    for( cn = ni.FirstP(); !(ni.AtEnd()); cn = ni.NextP() )
      WriteAllNodeData( cn );
  } else {
    // write in node ID order
    typename tMesh< tSubNode >::tIdArrayNode_t RNode(*(this->m->getNodeList()));
    int i;
    for( i=0; i<nActiveNodes; ++i )
      WriteActiveNodeData( RNode[i] );
    for( i=0; i<nnodes; ++i )
      WriteAllNodeData( RNode[i] );
  }

  // Write data specific for the stratGrid class
  // sections, gravel bodies and preservation potential
  if(time > 0 && stratOutput != 0){
    stratOutput->WriteNodeData( time, counter );
  }

  randomofs << std::flush;
  drareaofs << std::flush;
  netofs << std::flush;
  slpofs << std::flush;
  qofs << std::flush;
  texofs << std::flush;
  tauofs << std::flush;
  if( vegofs.good() ) vegofs << std::flush;
  if( flowdepofs.good() ) flowdepofs << std::flush;
  if( chanwidthofs.good() ) chanwidthofs << std::flush;
  if( flowpathlenofs.good() ) flowpathlenofs << std::flush;
  if( qsofs.good() ) qsofs << std::flush;
  if( qsinofs.good() ) qsinofs << std::flush;
  if( qsdinofs.good() ) qsdinofs << std::flush;
  if( upofs.good() ) upofs << std::flush;
  if( dzdtofs.good() ) dzdtofs << std::flush;

  if(OptLayOutput) layofs.close();
  if( surfofs.good() )
    surfofs.close();
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

/***********************************************************************\
**
** SetStratGrid: set stratGrid
**
** Created 6/2003 (QC)
\***********************************************************************/
template< class tSubNode >
void tLOutput<tSubNode>::SetStratGrid(tStratGrid *s_, tStreamNet *netPtr) {
  assert(stratOutput);
  stratOutput->SetStratGrid(s_, netPtr);
}

/*************************************************************************\
 **
 **  tTSOutputImp constructor
 **
 **  Creates and opens a series of files for time series
 **
 **  Modifications:
\*************************************************************************/
template< class tSubNode >
tTSOutputImp<tSubNode>::tTSOutputImp( tMesh<tSubNode> *meshPtr,
				      const tInputFile &infile ) :
  tOutputBase<tSubNode>( meshPtr, infile ),  // call base-class constructor
  mdLastVolume(0.)
{
  int opOpt;  // Optional modules: only output stuff when needed

  this->CreateAndOpenFile( &volsofs, ".vols" );
  this->CreateAndOpenFile( &dvolsofs, ".dvols" );
  if( (opOpt = infile.ReadItem( opOpt, "OPTVEG" ) ) != 0)
    this->CreateAndOpenFile( &vegcovofs, ".vcov" );
  this->CreateAndOpenFile( &tareaofs, ".tarea" );
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
  typename tMesh< tSubNode >::nodeListIter_t niter( this->m->getNodeList() ); // node list iterator

  tSubNode * cn;       // current node

  double volume = 0.,
    area = 0.,
    cover = 0.;

  if (0)//DEBUG
    std::cout << "tTSOutputImp::WriteTSOutput()" << std::endl;

  for( cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() ) {
    volume += cn->getZ()*cn->getVArea();
    area += cn->getVArea();
  }

  volsofs << volume << std::endl;
  if( mdLastVolume > 0.0 )
    dvolsofs << volume - mdLastVolume << std::endl;
  mdLastVolume = volume;
  //tareaofs << area << std::endl;

  if( vegcovofs.good() ) {
    for( cn = niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() )
      cover += cn->getVegCover().getVeg()*cn->getVArea();
    vegcovofs << cover/area << std::endl;
  }
}

/***********************************************************************\
**
** tStratOutputImp constructor
**
** Created 7/2003 (QC)
\***********************************************************************/
template< class tSubNode >
tStratOutputImp<tSubNode>::tStratOutputImp(tMesh<tSubNode> * meshPtr,
					   const tInputFile &infile ) :
  tOutputBase<tSubNode>( meshPtr, infile ),  // call base-class constructor
  stratGrid(0), netPtr(0)
{
  // Set direction of section
  {
    int s;
    for(s=0; s<nIsections;++s){
      direction[s] = DIRI;
    }
    for(s=nIsections; s<nSections;++s){
      direction[s] = DIRJ;
    }
  }

  // Z-elevations along the sections
  {
    for(int s=0; s<nSections;++s){
      char prefix[10];
      sprintf(prefix, ".xyz%d", s+1);
      this->CreateAndOpenFile( &xyzofs[s], prefix );
    }
  }

  // Texture sections 10 x
  {
    for(int s=0; s<nSections;++s){
      char prefix[10];
      sprintf(prefix, ".layt3%d", s+1);
      this->CreateAndOpenFile( &layt3ofs[s], prefix );
    }
  }

  // Facies sections 10 x
  {
    for(int s=0; s<nSections;++s){
      char prefix[10];
      sprintf(prefix, ".layt4%d", s+1);
      this->CreateAndOpenFile( &layt4ofs[s], prefix );
    }
  }

  // Subsurface gravel bodies xyz and layers
  this->CreateAndOpenFile( &grxyzofs, ".xyz");
  this->CreateAndOpenFile( &gravofs, ".grav");

  //Preservation potential files
  this->CreateAndOpenFile( &psurfofs,".presSurface");
  this->CreateAndOpenFile( &pssurfofs,".presSubsurface");
  this->CreateAndOpenFile( &pssurf2ofs,".presSubsurface2");
}

/***********************************************************************\
**
** WriteStratGridSections: Writes z elevation files and stratigraphy
** 			   files at user-defined locations, now every
**			   20 indices
**
** Created 6/2003 (QC)
\***********************************************************************/
template< class tSubNode >
void tStratOutputImp<tSubNode>::WriteStratGridSections(double time)
{
  // In the beginning write the section info to the files
  for(int s=0; s<nSections;++s){
    WriteSingleSection(time, s, xyzofs[s], layt3ofs[s], layt4ofs[s],
		       direction[s]);
  }
}

/***********************************************************************\
**
** WriteSingleSection: Writes a section
**
** Created 7/2003 (QC)
\***********************************************************************/
template< class tSubNode >
void tStratOutputImp<tSubNode>::WriteSingleSection(double time,
						   int section,
						   std::ofstream &xyzofs,
						   std::ofstream &layt3ofs,
						   std::ofstream &layt4ofs,
						   direction_t dir)
{
  const int vmax = dir == DIRI ? stratGrid->getJmax() : stratGrid->getImax();

  //--------------------------------
  // Write time and number of nodes
  // in every file, the Matlab routine
  // needs to know this
  //-------------------------------
  this->WriteTimeNumberElements(xyzofs, time, vmax);

  // Texture sections
  this->WriteTimeNumberElements(layt3ofs, time, vmax);

  // Facies sections
  this->WriteTimeNumberElements(layt4ofs, time, vmax);

  //---------------------------------------------------------------------
  // Section 'section'
  //---------------------------------------------------------------------
  const int vcst=stratGrid->getSectionLocation(section);
  tMatrix<tStratNode> const *StratNodeMatrix = stratGrid->getStratNodeMatrix();
  int v;
  for(v=0;v<vmax;v++){
    tStratNode const &sn = dir == DIRI ?
      (*StratNodeMatrix)(vcst,v) :
      (*StratNodeMatrix)(v,vcst);

    xyzofs << sn.getX() << ' ' << sn.getY() << ' ' << sn.getZ() << '\n';

    // Get the number of layers present
    const int numlayers = sn.getNumLayer();

    layt3ofs << numlayers-1 << '\n';
    layt4ofs << numlayers-1 << '\n';

    int l=1;
    while(l<numlayers){
      const double thickness = sn.getLayerDepth(l);
      const double texture =
	(thickness !=0.0) ?
	sn.getLayerDgrade(l,0)/thickness:
	0.;
      const double lasttime  = sn.getLayerRtime(l);
      const double depotime  = sn.getLayerCtime(l);
      //facies  = sn.getLayerFacies(l);
      const double paleocurrent = sn.getPaleoCurrent(l);

      layt3ofs << thickness << ' ' << texture<< ' '
	       << lasttime << ' ' << depotime << '\n';
      layt4ofs << thickness << ' ' << paleocurrent << ' '
	       << lasttime << ' ' << depotime << '\n';


      l++;   // increment, loop down the layer list.
    }
  } // end of loop over v

  std::cout << "Output::Finished writing section line " << section
       << "..." << '\n';
}

/************************************************************\
 WriteCompleteStratigraphy
 Functions which writes the complete stratigraphic subsurface
 matrix to a file, needed in Geo-Archeaology simulations
 
 Quintijn Clevis, may 2004, @CU Boulder
\*************************************************************/

template< class tSubNode >
void tStratOutputImp<tSubNode>::WriteCompleteStratigraphy(double time, int counter_)
{
	
	// Open the file for writing map-view data
  std::ofstream stratxyzofs;
  std::ofstream stratofs;
  std::ofstream lithoofs;
  std::ofstream channelofs;
  
	// File 1, containing the seperate x,y,z locations
#define MY_EXT ".stratxyz" 
  char ext[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( ext, "%s%d", MY_EXT, counter_ );
#undef MY_EXT
  this->CreateAndOpenFile( &stratxyzofs, ext );
  
  //File 2, containing the full stratigraphic matrix containing archaeology
#define MY_EXT ".strat" 
  char extb[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( extb, "%s%d", MY_EXT, counter_ );
#undef MY_EXT
  this->CreateAndOpenFile( &stratofs, extb );
  
  
  //File 3,containing the full stratigraphic matrix containing lithology
#define MY_EXT ".litho" 
  char extc[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( extc, "%s%d", MY_EXT, counter_ );
#undef MY_EXT
  this->CreateAndOpenFile( &lithoofs, extc );  
  
  
  //File 4, containing the channel locations
#define MY_EXT ".channelmap" 
  char extd[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( extd, "%s%d", MY_EXT, counter_ );
#undef MY_EXT
  this->CreateAndOpenFile( &channelofs, extd );
	
  const int imax = stratGrid->getImax();
  const int jmax = stratGrid->getJmax();
  const int numCells = ((imax-2)*(jmax-2));
  tMatrix<tStratNode> const *StratNodeMatrix = stratGrid->getStratNodeMatrix();
  tMatrix<tTriangle*> const *StratConnect = stratGrid->getStratConnect();
  
  stratxyzofs<< time << '\n' << numCells << '\n'<< imax-1 << '\n' <<jmax-1 << '\n';
  stratofs   << time << '\n' << numCells << '\n'<< imax-1 << '\n' <<jmax-1 << '\n';
  lithoofs   << time << '\n' << numCells << '\n'<< imax-1 << '\n' <<jmax-1 << '\n';
  channelofs << time << '\n' << numCells << '\n'<< imax-1 << '\n' <<jmax-1 << '\n';
  
  // Loop over the inner part of the StratGrid
  int i,j;
  for(j=1;j<jmax-1;j++){
    for(i=1;i<imax-1;i++){

      tStratNode const &sn = (*StratNodeMatrix)(i,j);

      // File 1; Write the location of the node to the xyz file
      stratxyzofs << sn.getX() << ' ' << sn.getY() << ' ' << sn.getZ() << '\n';
      
      // File 2; Write loaction and discharge info
      tTriangle *ct =   (*StratConnect)(i,j);									// fetch  Triangle
      const double sx = sn.getX();                						// i,j's  X-value
      const double sy = sn.getY();	        									// i,j's  Y-value
      
      tLNode *lnds[3];					// put the nodes in an array
      
       if(ct != NULL){				                    // If Triangle present
       	int p;
       	int numnodes=0;
       	for(p=0; p<=2; p++){
	      lnds[numnodes] = static_cast<tLNode *>(ct->pPtr(p));
	      numnodes++;
        }
      }
	     // and interpolate the discharge value
    	const tArray<double> tri_drainage(lnds[0]->getDrArea(),lnds[1]->getDrArea(),lnds[2]->getDrArea() );
      const double drainage = PlaneFit(sx, sy, lnds[0]->get2DCoords(), lnds[1]->get2DCoords(),
				     lnds[2]->get2DCoords(), tri_drainage );
				     
		  // Write File 2
		  channelofs << sn.getX() << ' ' << sn.getY() << ' ' << sn.getZ() << ' '<< drainage/1000.0 <<'\n';		     
				     
      
      // Get the number of layers present at this stratnode location
      const int numlayers = sn.getNumLayer();
      stratofs << numlayers-1 << '\n';
      lithoofs << numlayers-1 << '\n';
      
      int l=1;
      while(l<numlayers){
        const double thickness = sn.getLayerDepth(l);
        const double texture =(thickness !=0.0) ? sn.getLayerDgrade(l,0)/thickness: 0.;
        const double lasttime  = sn.getLayerRtime(l);
        const double depotime  = sn.getLayerCtime(l);
        //facies  = sn.getLayerFacies(l);
        //const double paleocurrent = sn.getPaleoCurrent(l);

        stratofs << thickness << ' ' << texture<< ' '
	       << lasttime << ' ' << depotime << '\n';
	       
	      lithoofs << thickness << ' ' << texture<< ' '
	       << lasttime << ' ' << depotime << '\n';
      


        l++;   // increment, loop down the layer list.
      }
      
    } //-i
  } //-j
  
  //stratxyzofs << flush;
  //stratofs << flush;
  
  stratxyzofs.close();
  stratofs.close();
  lithoofs.close();
  channelofs.close();
  
 
  std::cout << "Output::Finished writing complete subsurface stratigraphy file...\n";

} // end writing complete stratigraphy


/*************************************************************\
 tLOutput::WriteGravelBodies  Function that writes all layers
 that contain gravel into a file. This can be used for visual
 isation using Matlab routines

 Created 6/2003 (QC)
\**************************************************************/

template< class tSubNode >
void tStratOutputImp<tSubNode>::WriteGravelBodies(double time, int counter_)
{

  // Open the file for writing map-view data
  std::ofstream gravelmapofs;

#define MY_EXT ".gravelmap"
  char ext[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( ext, "%s%d", MY_EXT, counter_ );
#undef MY_EXT
  this->CreateAndOpenFile( &gravelmapofs, ext );

  const int imax = stratGrid->getImax();
  const int jmax = stratGrid->getJmax();
  const int numCells = ((imax-2)*(jmax-2));
  tMatrix<tStratNode> const *StratNodeMatrix = stratGrid->getStratNodeMatrix();

  grxyzofs<< time << '\n' << numCells << '\n';
  gravofs << time << '\n' << numCells << '\n';

  // Loop over the inner part of the StratGrid
  int i,j;
  for(j=1;j<jmax-1;j++){
    for(i=1;i<imax-1;i++){

      tStratNode const &sn = (*StratNodeMatrix)(i,j);

      // Write the location of the node
      grxyzofs << sn.getX() << ' ' << sn.getY() << ' ' << sn.getZ() << '\n';

      //--------------------Step 1------------------------------
      //Pre-count the number of suited layers at this position
      const int numlayers = sn.getNumLayer();
      int l;
      l=1;
      double overburden = -1.0;
      double totalgravelthickness = 0.0;
      double totaldepth = 0.0;
      int gravelcount=0;
      while (l<numlayers){
	const double thickness = sn.getLayerDepth(l);
	totaldepth+=thickness;
	const double texture = (thickness !=0.0) ?
	  sn.getLayerDgrade(l,0)/thickness:
	  0.0;

	const double ctime=sn.getLayerCtime(l);

	if((texture >= 0.9 && texture <= 1.0) && (thickness < 10. && ctime >0.)  ){
	  gravelcount++;
	  totalgravelthickness+=thickness;
	  if(overburden == -1){
	    overburden = sn.getZ() - (totaldepth-thickness);
	  }

	}

	l++;
      } //while over layers

      double axisdist = sn.getX() - ( netPtr->getInletNodePtr() )->getY();

      //Something in between......write the file containing the distribution of
      gravelmapofs<<sn.getX()<<' '<<sn.getY()<<' '<<sn.getZ()<<' '<<axisdist<<' '<<totalgravelthickness<<' '<<overburden<< '\n';



      //--------------------Step 2------------------------------
      // Matlab needs to now how many gravel layers are present
      gravofs << gravelcount << '\n';
      double ztop = sn.getZ();

      //--------------------Step 3----------------------------
      // Loop again but now write the layers to the file
      l=1;
      while (l<numlayers){
	const double thickness = sn.getLayerDepth(l);
	const double texture = (thickness !=0.0) ?
	  sn.getLayerDgrade(l,0)/thickness:
	  0.0;

	const double ctime=sn.getLayerCtime(l);

	const double zbase=ztop - thickness;
	if( (texture >= 0.9 && texture <= 1.0)   && (thickness < 10. && ctime > 0.)  ){
	  gravofs << texture << ' ' << ztop << ' ' << zbase << '\n';
	}

	ztop=zbase;
	l++;
      } //while over layers
    } //i
  } //j


  std::cout << "Output::Finished writing gravel bodies subsurface distribution...\n";
} // End of function WriteGravelBody

//-------------------------- a new function below --------------------------------------


/*************************************************************\
** WriteTimeBodies
**
** Function that writes the voxels of sediment with certain
** age bounds to a file. This is used later for visualisation
** using matlab
**
** This routine is called once every 500 yrs, so gives 20
** databatches
\**************************************************************/

// ROUTINE CODING IN PROGRESS NOT FINISHED ! Quintijn 6-11-2003
/*
template< class tSubNode >
void tStratOutputImp<tSubNode>::WriteTimeBodies(double time)
{

 tMatrix<tStratNode> const *StratNodeMatrix = stratGrid->getStratNodeMatrix();

 const int imax = stratGrid->getImax();
 const int jmax = stratGrid->getJmax();

 tArray<double>TimeSliceTop(20);
 tArray<double>TimeSliceBase(20);
 int t;
 nslices = 19;
 for(t=0; t<=19; t++){
  TimeSliceTop[t] =(t+1)*500;
  TimeSliceBase[t]=t*500;
 }



  // Loop over the StratGrid
  int i,j;
  for(i=1;i<imax-1;i++){
    for(j=1;j<jmax-1;j++){

      tStratNode const &sn = (*StratNodeMatrix)(i,j);
      const int numlayers = sn.getNumLayer();
      previous_depth =0.0;
      depth = 0.0;

      int l=1;
      while (l<numlayers){                  // loop downward in the stratigraphy

	for(t=0;t<=nslices;t++){      // loop 'upward' in time, from the first to the current write step

	  const double timeslicebase = TimeSliceBase(t);
	  const double timeslicetop  = TimeSliceTop(t);

	  tArray<double>ztop(20);   // contains for every timeslice the top elevation
	  ztop[t] = 1e32;           // initialize for the current slice step to 1e32, some rediculous value

	  previous_depth = depth;
	  depth = depth +  sn.getLayerDepth(l);                          // maximum real depth of the current layer under the grass...

	  //There are 3 options here, as stratigraphic layers are not always
	  //identical to the width of the timeslices

	  // 1) Both the layer's creation time and Recent activity time fall with the timeslice
	  if( (sn.getLayerCtime(l) >= timeslicebase) && (sn.getLayerCtime(l) <= timeslicetop)
	      && (sn.getLayerRtime(l) >=timeslicebase) && (sn.getLayerRtime(l) <= timeslicetop) )
	    {

	      TimeSliceThickness[t]  += sn.getLayerDepth(l);                        // Write Subsurface
	      if(ztop[t] == 1e32){
	         ztop[t] = previousdepth;	 // we did not touch this timeslice before, so previous depth is the real depth for this unit
	      }

	    }
	  // 2) layer is partly in this timeslice, but also partly in the timeslice above
	  // find the thickness of the part in this tmeslice using linear interpolation of
	  // the sedimentation rate
	  else if(sn.getLayerCtime(l) >= timeslicebase && sn.getLayerCtime(l) <= timeslicetop
	          && sn.getLayerRtime(l) > timeslicetop)
	    {
	    	if(sn.getLayerRtime(l) - sn.getLayerCtime(l) > 0.0 && sn.getLayerDepth(l)> 0.0){
	    	 double Rsed        =  sn.getLayerDepth(l)/(sn.getLayerRtime(l) - sn.getLayerCtime(l) );
	    	 double TimeInSlice =  (timeslicetop - sn.getLayerCtime(l));       // thickness deposited in this timeslice

	    	 TimeSliceThickness[t]  += TimeInSlice*Rsed;                           // add part
	    	 if(ztop[t] == 1e32){
	            ztop[t] = previousdepth;  // we did not touch this timeslice before, so previous depth is the real depth for this unit
	         }


	    	} // both times are ok
	    }
	 // 3) layer is in this slice but also in the timeslice below
	 else if(sn.getLayerCtime(l) < timeslicebase && sn.getLayerRtime(l) > timeslicebase && sn.getLayerRtime(l) <=timeslicetop){

	        if(sn.getLayerRtime(l) - sn.getLayerCtime(l) > 0.0 && sn.getLayerDepth(l)> 0.0){
	         double Rsed        =  sn.getLayerDepth(l)/(sn.getLayerRtime(l) - sn.getLayerCtime(l) );
	         double TimeInSlice =  sn.getLayerRtime(l)-timeslicebase;          // thickness deposited in this timeslice


	         TimeSliceThickness[t] +=  TimeInSlice*Rsed;   	                  // add part

	         if(ztop[t] == 1e32){
	            ztop[t] = previousdepth;   // we did not touch this timeslice before, so previous depth is the real depth for this unit
	         }


	        }
	   }

       } // loop over timeslices
      // We also want to know how the subsurface-age distribution looks like

       // continue with the next layer.........
       l++;

      } // loop layers

      */
      /********************************************************************************\
      ** OK, we now have two important arrays.....
      ** 1) TimeSliceThickness[t] contains the thickness of sediment with a certan 500 yrs timeslice
      ** 2) ztop[t] contains the top elevation for this timeslice
      **
      ** Loop again over the timeslices and write the information to file
      \*********************************************************************************/
  /*

      for(t=0;t<=nslices;t++){



      }



    } // i    loop over the stratgrid cell locations
  } // j

  }

//--------------------finished new function above -------------------------

*/

/**************************************************************\
 WritePreservationPotential

 This files assumes that the statigraphy can be devided in time
 slices of ~size 'writestep', e.g 100-1000 yrs. The function checks
 how much (%,fraction) of the originally sediment deposited in a
 specific timeslice is preserved in subsequent timeslices.

 -Modification 1-2-2004(QC)
 Preservation potential values are generally high and
 less sensitive to perturbations in the forcing when the preservation
 ratio is calculted on the entire floodplain, an area which is 3,5 times
 wider than a final channel belt, e.g the zone of meandering.
 Ratios are therefore now also calculted on a smaller zone with a with
 close to the channel, meander belt width (~1000, 1500m)

 Created 06/2003 (QC)
 \*************************************************************/

template< class tSubNode >
void tStratOutputImp<tSubNode>::WritePreservationPotential(double time,
							   int counter_)
{
  // Surfer style file containing stratnodes x,y,z, mindist to meander
  // and toplayer's Ct, Rt, and Thickn.
  double depth;
  double mbelt_Ytop;
  double mbelt_Ybase;

  mbelt_Ytop  = 21000.0;   //y-coordinate narrow channel belt/ meander zone
  mbelt_Ybase = 18900.0;   //of width 1200m


  std::ofstream topofs;

#define MY_EXT ".top"
  char ext[sizeof(MY_EXT)+10];  // name of file to be created

  sprintf( ext, "%s%d", MY_EXT, counter_ );
#undef MY_EXT
  this->CreateAndOpenFile( &topofs, ext );

  //File header, time and what's in the column below
  topofs<<time<<'\n';
  topofs<<"X"<<' '<<"Y"<<' '<<"Z"<<' '<<"mind"<<' '<<"axisd"<<' '
	<<"Ct1"<<' '<<"Rt1"<<' '
	<<"D1"<<' '<<"DA1"<<' '<<"DA2"<<' '<<"DA3"<<' '<<"DA4"<<'\n';
  tMatrix<tStratNode> const *StratNodeMatrix = stratGrid->getStratNodeMatrix();

  const int current_ts = int( ROUND( time/stratGrid->getnWrite() ) );
  //std::cout<<"Current-Ts= " << current_ts << '\n';

  tArray<double> surface(current_ts+1);
  tArray<double> subsurface(current_ts+1);		// material in the entire floodplain zone
  tArray<double> subsurface_mbelt(current_ts+1);	// material in the narrow meanderbelt zone
  tArray<double> depth_age(15);
  tLNode *meandernode = NULL;

  const int imax = stratGrid->getImax();
  const int jmax = stratGrid->getJmax();

  // Reset the values for the two preservation arrays
  int ts;
  for(ts=0;ts<=current_ts;ts++){
    surface[ts]   = 0.;
    subsurface[ts]= 0.;
    subsurface_mbelt[ts] = 0.;
  }

  stratGrid->setOutputTime(current_ts, time);

  // Loop over the StratGrid
  for(int i=1;i<imax-1;i++){
    for(int j=1;j<jmax-1;j++){

      //reset the array containg the sediment age at depths, 1,3,4,5 meter.
      for(int d=0;d<=10;d++){
	depth_age[d] = 0.0;
      }

      tStratNode const &sn = (*StratNodeMatrix)(i,j);
      const int numlayers = sn.getNumLayer();
      depth = 0.0;

      // FIXME could be a loop on tStratNode::layerlist
      int l=1;
      while (l<numlayers){                  // loop downward in the stratigraphy

	for(ts=1;ts<=current_ts;ts++){      // loop 'upward' in time, from the first to the current write step

	  const double timeslicebase = stratGrid->getOutputTime(ts-1);
	  const double timeslicetop  = stratGrid->getOutputTime(ts);
	  depth += sn.getLayerDepth(l);     // maximum real depth of the current layer under the grass...

	  //There are 3 options here, as stratigraphic layers are not always
	  //identical to the width of the timeslices

	  // 1) Both the layer's creation time and Recent activity time fall with the timeslice
	  const double thisLayerCtime = sn.getLayerCtime(l);
	  const double thisLayerRtime = sn.getLayerRtime(l);
	  const double thisLayerDepth = sn.getLayerDepth(l);

	  if( (thisLayerCtime >= timeslicebase) && (thisLayerCtime <= timeslicetop)
	      && (thisLayerRtime >=timeslicebase) && (thisLayerRtime <= timeslicetop) )
	    {
	      subsurface[ts]      += thisLayerDepth;                        // Write Subsurface, entire floodplain
	      if(l==1) surface[ts]+= thisLayerDepth;

	      if(sn.getY() >= mbelt_Ybase && sn.getY() <= mbelt_Ytop){
	      	subsurface_mbelt[ts] += thisLayerDepth;
	      }

	      //DEBUG, trying to track origin of negative values
	      if(subsurface[ts] < 0.0){
		std::cout<<"Error in StratOutput::preservation potential, type 1 "<<std::endl;
		std::cout<<"at ts ="<<ts<< " layer "<< l << "subsurface value = "<<subsurface[ts]<<std::endl;
		std::cout<<"timeslicetop= "<<timeslicetop<<std::endl;
		std::cout<<"timeslicebase= "<<timeslicebase<<std::endl;
		std::cout<<"layerdepth ="<<	thisLayerDepth;
		exit(1);
	      }


	    }
	  // 2) layer is partly in this timeslice, but also partly in the timeslice above
	  // find the thickness of the part in this tmeslice using linear interpolation of
	  // the sedimentation rate
	  else if(thisLayerCtime >= timeslicebase && thisLayerCtime <= timeslicetop
	          && thisLayerRtime > timeslicetop)
	    {
	      if(thisLayerRtime - thisLayerCtime > 0.0 && thisLayerDepth> 0.0){
		const double Rsed        =  thisLayerDepth/(thisLayerRtime - thisLayerCtime );
		const double TimeInSlice =  (timeslicetop - thisLayerCtime);       // thickness deposited in this timeslice
		subsurface[ts]      += TimeInSlice*Rsed;                     // add part
		if(l==1) surface[ts] +=  TimeInSlice*Rsed;                   // surface layer

		if(sn.getY() >= mbelt_Ybase && sn.getY() <= mbelt_Ytop){
		  subsurface_mbelt[ts] += TimeInSlice*Rsed;
		}

		if(subsurface[ts] < 0.0  ){     // DEBUG
		  std::cout<<"Error in StratOutput::write preservation potential, type 2 "<<std::endl;
		  std::cout<<"at ts ="<<ts<< " layer "<< l << "subsurface value = "<<subsurface[ts]<<std::endl;
		  std::cout<<"timeslicetop= "<<timeslicetop<<" timeslicebase= "<<timeslicebase<<std::endl;
		  std::cout<<"layerR time= "<<thisLayerRtime<<" layerC time= "<<thisLayerCtime<<std::endl;
		  std::cout<<"TimeInSlice= "<<TimeInSlice<<std::endl;
		  std::cout<<"layerdepth ="<<	thisLayerDepth;
		  std::cout<<"subsurface[ts]= "<<subsurface[ts] <<std::endl;
		  exit(1);
		}

	      } // both times are ok
	    }
	  // 3) layer is in this slice but also in the timeslice below
	  else if(thisLayerCtime < timeslicebase && thisLayerRtime > timeslicebase && thisLayerRtime <=timeslicetop){

	    if(thisLayerRtime - thisLayerCtime > 0.0 && thisLayerDepth> 0.0){
	      const double Rsed        =  thisLayerDepth/(thisLayerRtime - thisLayerCtime );
	      const double TimeInSlice =  thisLayerRtime-timeslicebase;         // thickness deposited in this timeslice
	      subsurface[ts]      +=  TimeInSlice*Rsed;                   // add part
	      if(l==1) surface[ts] +=  TimeInSlice*Rsed;                  // surface layer

	      if(sn.getY() >= mbelt_Ybase && sn.getY() <= mbelt_Ytop){
		subsurface_mbelt[ts] += TimeInSlice*Rsed;
	      }

	      if(subsurface[ts] < 0.0 ){     // DEBUG
		std::cout<<"Error in StratOutput::write preservation potential, type 3 "<<std::endl;
		std::cout<<"at ts ="<<ts<< " layer "<< l << "subsurface value = "<<subsurface[ts]<<std::endl;
		std::cout<<"timeslicetop= "<<timeslicetop<<" timeslicebase= "<<timeslicebase<<std::endl;
		std::cout<<"layerR time= "<<thisLayerRtime<<" layerC time= "<<thisLayerCtime<<std::endl;
		std::cout<<"TimeInSlice= "<<TimeInSlice<<std::endl;
		std::cout<<"layerdepth ="<<	thisLayerDepth;
		std::cout<<"subsurface[ts]= "<<subsurface[ts] <<std::endl;
		exit(1);
	      }
	    }
	  }

	} // loop over timeslices
	// We also want to know how the subsurface-age distribution looks like

	// continue with the next layer.........
	l++;
      } // loop layers

      /*********************************************************************************************\
       // Write the properties of the top layer to separate file, we want to know them as a funtion
       // of the distance to the main meandering channel, therefore a costly linear search, sorry....
      \*********************************************************************************************/

      int counter=0;
      double dist    = 1000000.0;
      double mindist = 1000000.0;
      assert( netPtr != 0);
      meandernode = netPtr->getInletNodePtrNC();
      while(meandernode != NULL && counter < 300){                 // nb wil not work if num meander nodes > 300 !
	const double XX = sn.getX()- meandernode->getX();
	const double YY = sn.getY()- meandernode->getY();
	dist = sqrt( XX*XX + YY*YY );
	if(dist < mindist){
	  mindist=dist;
	}

	meandernode=meandernode->getDownstrmNbr();
	counter++;
	if(meandernode->getBoundaryFlag() != kNonBoundary) break;
      }
      /*****************************************************\
       // And,...what is the distance to the channel axis ?
       // in case this channel axis runs parralel to the x-axis
       // of the StratGrid
       // axisdist = sn->getY() - getInletNodePtr()->getY
       // if axisdist = negative, soutthern floodplain half
       // if axisdist = positve, southern floodplain half
       //
       // Oeps, be carefull
      \******************************************************/
      double axisdist = sn.getX() - ( netPtr->getInletNodePtr() ->getY());

      depth_age[1] = sn.getAgeAtDepth( 0.5 );
      depth_age[2] = sn.getAgeAtDepth( 1.0 );
      depth_age[3] = sn.getAgeAtDepth( 2.0 );
      depth_age[4] = sn.getAgeAtDepth( 3.0 );
      depth_age[5] = sn.getAgeAtDepth( 4.0 );
      depth_age[6] = sn.getAgeAtDepth( 5.0 );
      depth_age[7] = sn.getAgeAtDepth( 6.0 );
      depth_age[8] = sn.getAgeAtDepth( 7.0 );
      depth_age[9] = sn.getAgeAtDepth( 8.0 );
      depth_age[10] = sn.getAgeAtDepth( 9.0 );

      topofs <<sn.getX()<< ' '<<sn.getY()<< ' '<<sn.getZ()<< ' '
	     << mindist << ' ' << axisdist<<' '
	     << sn.getLayerCtime(1) << ' ' <<sn.getLayerRtime(1) << ' '<<sn.getLayerDepth(1)<<' '
	     <<depth_age[1]<<' '<<depth_age[2]<<' '<<depth_age[3]<<' '
	     <<depth_age[4]<<' '<<depth_age[5]<<' '<<depth_age[6]<<' '
	     <<depth_age[7]<<' '<<depth_age[8]<<' '<<depth_age[9]<<' '
	     <<depth_age[10]<<'\n';



    } //     loop over the stratgrid cell locations
  }

  // Set the base values for the current timestep
  stratGrid->setSurface(current_ts, surface[current_ts]);
  stratGrid->setSubsurface(current_ts, subsurface[current_ts]);
  stratGrid->setSubsurface_mbelt(current_ts, subsurface_mbelt[current_ts]);

  /**********************************************************************\
    Complement the two files, for the surface and subsurface
    preservation potential
    File structure is: time  p_of_unit1  p_of_unit2 ...
                       time2 p_of_unit1  p_of_unit2 p_of_unit3...
                       time3 p_of_unit1  p_of_unit2 p_of_unit3 p_of_unit4

 \***********************************************************************/
  psurfofs << time;                                    // time1..2..3 etc.
  pssurfofs << time;
  for(ts=1;ts<=current_ts;ts++){
    //---------------------------------------
    const double surface_preservation =
      (stratGrid->getSurface(ts) > 0.)?
      surface[ts]/stratGrid->getSurface(ts):
      0.;
    //---------------------------------------
    const double subsurface_preservation =
      (stratGrid->getSubsurface(ts) > 0.)?
      subsurface[ts]/stratGrid->getSubsurface(ts):
      0.;

    psurfofs  << ' ' << surface_preservation;
    pssurfofs << ' ' << subsurface_preservation;
  }

  psurfofs << '\n';                                         // end of line in file
  pssurfofs<< '\n';                                         // end of line in file

  /************************************************************************\
   // Write a file with the format  -Timeslice1  -preservation -totalvolume
   //				   -Timeslice2  -preservation -totalvolume
   //				   -Timeslice3  -preservation .....
  \************************************************************************/
  pssurf2ofs.close();					    // close the existing pssurf2ofs file
  this->CreateAndOpenFile( &pssurf2ofs,".presSubsurface2");       // overwrites the existing file (emty & rewrite)

  pssurf2ofs<<time<<'\n';
  pssurf2ofs<<"T_unit"<<' '<<"p_fldpl"<<' '<<"vol_fldpl"<<' '<<"p_mbelt"<<' '<<"vol_mbelt"<<'\n';       // file header

  for(ts=1;ts<=current_ts;ts++){
    const double subsurface_preservation = (stratGrid->getSubsurface(ts) > 0.)?
      subsurface[ts]/stratGrid->getSubsurface(ts):0.;

    const double subsurface_preservation_mbelt = (stratGrid->getSubsurface_mbelt(ts) > 0.)?
      subsurface_mbelt[ts]/stratGrid->getSubsurface_mbelt(ts):0.;

    double theTime  = stratGrid->getOutputTime(ts);

    //assuming stratgrid cell dx is 50.0m, for the volume calulation
    pssurf2ofs<<theTime<<' '<<subsurface_preservation<<' '<<subsurface[ts]*50.0*50.0<<' '<<
      subsurface_preservation_mbelt<<' '<<subsurface_mbelt[ts]*50.0*50.0<<'\n';

  }
  if (1) //DEBUG
    std::cout << "Output::Finished writing Preservation potential of fluvial units..." << '\n';
} //Presevation potential

/***********************************************************************\
**
** SetStratGrid: set stratGrid
**
** Created 6/2003 (QC)
\***********************************************************************/
template< class tSubNode >
void tStratOutputImp<tSubNode>::SetStratGrid(tStratGrid *s_, tStreamNet *netPtr_) {
  stratGrid = s_;
  netPtr = netPtr_;
}

/***********************************************************************\
**
** Flush
**
** Created 6/2003 (QC)
\***********************************************************************/
template< class tSubNode >
void tStratOutputImp<tSubNode>::Flush() {
  for(int s=0; s<nSections;++s){
    xyzofs[s] << std::flush;
    layt3ofs[s] << std::flush;
    layt4ofs[s] << std::flush;
  }
  grxyzofs << std::flush;
  gravofs << std::flush;
  psurfofs << std::flush;
  pssurfofs << std::flush;
  pssurf2ofs << std::flush;
}

/*************************************************************************\
**
**  tStratOutputImp::WriteNodeData
**
\*************************************************************************/
template< class tSubNode >
void tStratOutputImp<tSubNode>::WriteNodeData( double time, int counter )
{
  WriteStratGridSections(time);
  WriteGravelBodies(time, counter);
  WriteCompleteStratigraphy(time, counter);
  WritePreservationPotential(time, counter);
  Flush();
}
