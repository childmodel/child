/*************************************************************************\
**
**  tOutput.cpp: Functions for output objects
**
**  
**
**  $Id: tOutput.cpp,v 1.1 1998-01-21 01:25:43 gtucker Exp $
\*************************************************************************/

#include "tOutput.h"
#include <string.h>

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
tOutput::tOutput( tGrid<tNode> * gridPtr, tInputFile &infile )
{
   char fullName[kMaxNameSize+6];
   
   assert( gridPtr > 0 );
   g = gridPtr;

   infile.ReadItem( baseName, "OUTFILENAME" );
   strcpy( fullName, baseName );
   strcat( fullName, ".nodes" );
   nodeofs.open( fullName );
   strcpy( fullName, baseName );
   strcat( fullName, ".edges" );
   edgofs.open( fullName );
   strcpy( fullName, baseName );
   strcat( fullName, ".tri" );
   triofs.open( fullName );
   strcpy( fullName, baseName );
   strcat( fullName, ".z" );
   zofs.open( fullName );

   if( !nodeofs.good() || !edgofs.good() || !triofs.good() || !zofs.good() )
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
void tOutput::WriteOutput( float time )
{
/*   tGridListIter<tNode> niter( g->GetNodeList() );
   tGridListIter<tEdge> eiter( g->GetEdgeList() );
   tListIter<tTriangle> titer( g->GetTriList() );*/
   tNode * cn;
   tEdge * ce;
   tTriangle * ct;
   int id;
   int nnodes = g->GetNodeList()->getSize();
/*   int nedges = g->GetEdgeList()->getSize();
   int ntri = g->GetTriList()->getSize();*/
   
   // Renumber IDs in order by position on list
/*   for( cn=niter.FirstP(), id=0; id<nnodes; cn=niter.NextP(), id++ )
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
              << cn->GetEdg()->getID() << " " << cn->getBoundaryFlag() << endl;
      zofs << cn->getZ() << endl;
   }
   
   // Write edge file
   edgofs << " " << time << endl << nedges << endl;
   for( ce=eiter.FirstP(); !(eiter.AtEnd()); ce=eiter.NextP() )
      edgofs << ce->getOriginPtrNC()->getID() << " "
             << ce->getDestinationPtrNC()->getID() << " "
             << ce->GetCCWEdg()->getID() << endl;
   
   // Write triangle file
   triofs << " " << time << endl << ntri << endl;
   for( ct=titer.FirstP(); !(titer.AtEnd()); ct=titer.NextP() )
       triofs << ct->pPtr(0)->getID() << " " 
              << ct->pPtr(1)->getID() << " " 
              << ct->pPtr(2)->getID() << " " 
              << ct->tPtr(0)->getID() << " " 
              << ct->tPtr(1)->getID() << " " 
              << ct->tPtr(2)->getID() << " " 
              << ct->ePtr(0)->getID() << " " 
              << ct->ePtr(1)->getID() << " " 
              << ct->ePtr(2)->getID() << endl;*/
   
   
}
