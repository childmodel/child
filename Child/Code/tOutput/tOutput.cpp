/*************************************************************************\
**
**  tOutput.cpp: Functions for output objects
**
**  
**
**  $Id: tOutput.cpp,v 1.10 1998-06-17 23:37:33 nmgaspar Exp $
\*************************************************************************/

#include "tOutput.h"

/*************************************************************************\
**
**  Constructor
**
**  The constructor takes two arguments, a pointer to the grid mesh and
**  a reference to an open input file. It reads the base name for the
**  output files from the input file, and opens and initializes these.
**
**  Input: gridPtr -- pointer to a tGrid object (or descendant), assumed
**                    valid
**         infile -- reference to an open input file, assumed valid
**
\*************************************************************************/
template< class tSubNode >
tOutput<tSubNode>::tOutput( tGrid<tSubNode> * gridPtr, tInputFile &infile )
{
   assert( gridPtr > 0 );
   g = gridPtr;

   infile.ReadItem( baseName, "OUTFILENAME" );

   CreateAndOpenFile( &nodeofs, ".nodes" );
   CreateAndOpenFile( &edgofs, ".edges" );
   CreateAndOpenFile( &triofs, ".tri" );
   CreateAndOpenFile( &zofs, ".z" );
   
   
}


/*************************************************************************\
**
**  CreateAndOpenFile
**
**  Opens the output file stream pointed by theOFStream, giving it the
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
void tOutput<tSubNode>::CreateAndOpenFile( ofstream *theOFStream,
                                           char *extension )
{
   char fullName[kMaxNameSize+6];
   
   strcpy( fullName, baseName );
   strcat( fullName, extension );
   theOFStream->open( fullName );

   if( !theOFStream->good() )
       ReportFatalError(
           "I can't create files for output. Memory may be exhausted." );
        
}




/*************************************************************************\
**
**  WriteOutput
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
   tGridListIter<tSubNode> niter( g->getNodeList() );
   tGridListIter<tEdge> eiter( g->getEdgeList() );
   tListIter<tTriangle> titer( g->getTriList() );
   tNode * cn;
   tEdge * ce;
   tTriangle * ct;
   int id;
   int nnodes = g->getNodeList()->getSize();
   int nedges = g->getEdgeList()->getSize();
   int ntri = g->getTriList()->getSize();

   cout << "tOutput::WriteOutput()\n";
   
   // Renumber IDs in order by position on list
   for( cn=niter.FirstP(), id=0; id<nnodes; cn=niter.NextP(), id++ )
       cn->setID( id );
   for( ce=eiter.FirstP(), id=0; id<nedges; ce=eiter.NextP(), id++ )
       ce->setID( id );
   for( ct=titer.FirstP(), id=0; id<ntri; ct=titer.NextP(), id++ )
       ct->setID( id );

   // Write node file and z file
   nodeofs << " " << time << endl << nnodes << endl;
   zofs << " " << time << endl << nnodes << endl;
   for( cn=niter.FirstP(); !(niter.AtEnd()); cn=niter.NextP() )
   {
      nodeofs << cn->getX() << " " << cn->getY() << " "
              << cn->getEdg()->getID() << " " << cn->getBoundaryFlag() << endl;
      zofs << cn->getZ() << endl;
   }
   
   // Write edge file
   edgofs << " " << time << endl << nedges << endl;
   for( ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() )
      edgofs << ce->getOriginPtrNC()->getID() << " "
             << ce->getDestinationPtrNC()->getID() << " "
             << ce->getCCWEdg()->getID() << endl;
   
   // Write triangle file
   int i;
   triofs << " " << time << endl << ntri << endl;
   for( ct=titer.FirstP(); !(titer.AtEnd()); ct=titer.NextP() )
   {
      for( i=0; i<=2; i++ )
          triofs << ct->pPtr(i)->getID() << " ";
      for( i=0; i<=2; i++ )
      {
          if( ct->tPtr(i) ) triofs << ct->tPtr(i)->getID() << " ";
          else triofs << "-1 ";
      }
      triofs << ct->ePtr(0)->getID() << " " 
             << ct->ePtr(1)->getID() << " " 
             << ct->ePtr(2)->getID() << endl;
   }
   
   WriteNodeData( time );
   
}

template< class tSubNode >
void tOutput<tSubNode>::WriteNodeData( double time ) 
{}


template< class tSubNode >
tLOutput<tSubNode>::tLOutput( tGrid<tSubNode> *g, tInputFile &infile ) 
        : tOutput<tSubNode>( g, infile )  // call base-class constructor
{
   CreateAndOpenFile( &drareaofs, ".area" );
   CreateAndOpenFile( &netofs, ".net" );
   CreateAndOpenFile( &slpofs, ".slp" );
   CreateAndOpenFile( &qofs, ".q" );
   CreateAndOpenFile( &layofs, ".lay" );
   CreateAndOpenFile( &texofs, ".tex" );
   
}


//TODO: should output boundary points as well so they'll map up with nodes
// for plotting. Means changing getSlope so it returns zero if flowedg
// undefined
template< class tSubNode >
void tLOutput<tSubNode>::WriteNodeData( double time )
{
   tGridListIter<tSubNode> ni( g->getNodeList() );
   tSubNode *cn;
   int nActiveNodes = g->getNodeList()->getActiveSize();
   int nnodes = g->getNodeList()->getSize();
   
   drareaofs << " " << time << "\n " << nActiveNodes << endl;
   netofs << " " << time << "\n " << nActiveNodes << endl;
   slpofs << " " << time << "\n " << nActiveNodes << endl;
   qofs << " " << time << "\n " << nnodes << endl;
   layofs << " " << time << "\n" << nnodes << endl;
   texofs << " " << time << "\n" << nnodes << endl;
   
   for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
   {
      assert( cn>0 );
      drareaofs << cn->getDrArea() << endl;
      if( cn->getDownstrmNbr() )
          netofs << cn->getDownstrmNbr()->getID() << endl;
      slpofs << cn->getSlope() << endl;
   }

   int i, j;
   for( cn = ni.FirstP(); !(ni.AtEnd()); cn = ni.NextP() ){
      qofs << cn->getQ() << endl;
      texofs << cn->getLayerDgrade(0,0)/cn->getLayerDepth(0) << endl;      
      i=0;
      while(i<cn->getNumLayer()){
         layofs << i+1 << " " << cn->getLayerCtime(i) << " " << cn->getLayerRtime(i) << endl;
         layofs << cn->getLayerDepth(i) << " " << cn->getLayerErody(i) << " " << cn->getLayerSed(i) << endl;
         if(cn->getLayerSed(i)>0){
            j=0;
            while(j<cn->getNumg()){
               layofs << cn->getLayerDgrade(i,j) << " ";
               j++;
            }
            layofs << endl;
         }
         i++;
      }
      
   }
   
}

