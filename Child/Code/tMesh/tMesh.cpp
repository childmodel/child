/***************************************************************************\
**
**  tMesh.cpp: Functions for class tMesh (see tMesh.h) plus global
**             functions used by tMesh methods (formerly tGrid)
**
**  Summary of recent changes:
**    - modified LocateTriangle to implement triangle search
**      starting from a given location; modified constructors to set
**      initial value of mSearchOriginTriPtr, and modified ExtricateTri...
**      to avoid dangling ptr. GT, 1/2000
**    - added initial densification functionality, GT Sept 2000
**
**  $Id: tMesh.cpp,v 1.107 2002-06-24 14:03:16 arnaud Exp $
\***************************************************************************/

#ifndef __GNUC__
#include "tMesh.h"
#endif

#include <stdlib.h>

#include "ParamMesh_t.h"

/***************************************************************************\
**  Templated global functions used by tMesh here
\***************************************************************************/

/***************************************************************************\
**
**  Next3Delaunay
**
**  global function; determines whether nbr node currently pointed to
**  by iterator and the next two in the nbr list form a Delaunay triangle.
**
**  Inputs:  nbrList -- list of pointers to nodes
**           nbrIter -- iterator for this list
**  Returns: 1 if the next 3 nodes on the list are Delaunay, 0 otherwise
**  Called by: tMesh::RepairMesh
**
\***************************************************************************/
template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &nbrList,
                   tPtrListIter< tSubNode > &nbrIter )
{
   static int ncalls = 0;
   ncalls++;
   tSubNode *cn, *nbrnd;
   
   assert( (&nbrList != 0) && (&nbrIter != 0) );

   nbrnd = nbrIter.DatPtr();
     //cout << "Next3Delaunay? no calls =  " << ncalls << endl;
     //global function prototypes
     //int PointsCCW( tArray< double > &, tArray< double > &,
     //           tArray< double > & );
     //int TriPasses( tArray< double > &, tArray< double > &,
     //           tArray< double > &, tArray< double > & );
     //cout << "N3D: pt a\n";
   tPtrListIter< tSubNode > nbrIterCopy( nbrList );
     //cout << "N3D: pt b\n";
   int i = nbrIter.Where();
     //cout << "N3D: pt c\n";
   nbrIterCopy.Get(i);
     //cout << "N3D: point d\n";
   // tPtrListNode< tSubNode > *tempptrnode = nbrIter.NodePtr();
   tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   tArray< double > p1( nbrIterCopy.NextP()->get2DCoords() );
   tArray< double > p2( nbrIterCopy.NextP()->get2DCoords() );
     //cout << "N3D: point B\n";

   // If points aren't counter-clockwise, we know it's not Delaunay
   if( !PointsCCW( p0, p1, p2 ) ) return 0;

   // Otherwise, compare it to each of the other potential triangles
   // p0-p1-ptest (?) where ptest is one of the other points in the
   // ring
   tArray< double > ptest;
   cn = nbrIterCopy.NextP();  // Move to next point in the ring
   while( cn != nbrnd )       // Keep testing 'til we're back to p0
   {
      ptest = cn->get2DCoords();
      if( !TriPasses( ptest, p0, p1, p2 ) )
      {
           //cout << "Next3Delaunay? No" << endl;
         return 0;
      }
        //else cout << "Next3Del? this tri passed..\n";
      cn = nbrIterCopy.NextP();  // Next point in ring
   }
   //cout << "Next3Delaunay? Yes" << endl;
   return 1;
}


/***************************************************************************\
**
**  PointAndNext2Delaunay
**
**  Global function that determines whether nbr node currently pointed to
**  by iterator and the next two in the nbr list form a Delaunay triangle.
**  Similar to Next3Delaunay but p2 is an arbitrary node (testNode) rather
**  than one of the neighbor list nodes.
**
**  Inputs:  testNode -- a node to check (this is "p2")
**           nbrList -- list of pointers to nodes
**           nbrIter -- iterator for this list
**  Returns: 1 if they are Delaunay, 0 otherwise
**
**  IS THIS EVER CALLED? OBSOLETE? TODO -- delete if not
**
\***************************************************************************/
template< class tSubNode >
int PointAndNext2Delaunay( tSubNode &testNode, tPtrList< tSubNode > &nbrList,
                           tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) && (&testNode != 0) );
   //cout << "PointAndNext2Delaunay?" << endl;
     //global function prototypes
     //int PointsCCW( tArray< double > &, tArray< double > &,
     //           tArray< double > & );
     //int TriPasses( tArray< double > &, tArray< double > &,
     //           tArray< double > &, tArray< double > & );
   tPtrListIter< tSubNode > nbrIterCopy( nbrList );
   int i = nbrIter.Where();
     //cout << "Where: " << i << endl;
   nbrIterCopy.Get( i );
     //assert( nbrIterCopy.Get( i ) );
   assert( nbrIterCopy.DatPtr() == nbrIter.DatPtr() );
   // tPtrListNode< tSubNode > *tempptrnode = nbrIter.NodePtr();
   tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   assert( nbrIterCopy.Next() );
   tArray< double > p1( nbrIterCopy.DatPtr()->get2DCoords() );
   tArray< double > p2( testNode.get2DCoords() );

   // If the points aren't CCW then we know it's not Delaunay
   if( !PointsCCW( p0, p1, p2 ) ) return 0;

   // Otherwise, call TriPasses to compare
   tArray< double > ptest;
   assert( nbrIterCopy.Next() );
   while( nbrIterCopy.DatPtr() != nbrIter.DatPtr() )
   {
      ptest = nbrIterCopy.DatPtr()->get2DCoords();
      if( !TriPasses( ptest, p0, p1, p2 ) )
      {
         //cout << "No" << endl;
         return 0;
      }
      assert( nbrIterCopy.Next() );
   }
   //cout << "Yes" << endl;
   return 1;
}


/**************************************************************************\
**  FUNCTIONS FOR CLASS tMesh
\**************************************************************************/

//default constructor
template< class tSubNode >     //tMesh
tMesh< tSubNode >::
tMesh() 
  :
  nnodes(0),
  nedges(0),
  ntri(0),
  seed(0),
  layerflag(FALSE),
  miNextNodeID(0),
  miNextEdgID(0),
  miNextTriID(0),   
  mSearchOriginTriPtr(0)
{
   cout<<"tMesh()"<<endl;
}

//copy constructor (created 11/99, GT)
//WARNING: this constructor relies on the behavior of assignment operations
//in tMeshList, tList, and tPtrList, which may or not give the desired
//results! Caveat emptor! -GT
template< class tSubNode >
tMesh<tSubNode>::tMesh( tMesh *originalMesh )
  :
  nnodes(originalMesh->nnodes),
  nedges(originalMesh->nedges),
  ntri(originalMesh->ntri),
  nodeList(originalMesh->nodeList),
  edgeList(originalMesh->edgeList),
  triList(originalMesh->triList),
  miNextNodeID(originalMesh->miNextNodeID),
  miNextEdgID(originalMesh->miEdgNodeID),
  miNextTriID(originalMesh->miNextTriID),   
  seed(originalMesh->seed),
  layerflag(originalMesh->layerflag),
  mSearchOriginTriPtr(0)
{}


/**************************************************************************\
**
**   tMesh( infile ): Reads from infile whether it is to reconstruct a mesh
**                    from input, construct a mesh from a list of (x,y,z,b)
**                    points (b=boundary code), or construct a mesh from
**                    scratch, given parameters in infile.
**
**   Created: 2/11/98, SL
**   Calls: tInputFile::ReadItem,
**        MakeMeshFromInputData( infile ), MakeMeshFromScratch( infile ),
**        of MakeMeshFromPoints( infile )
**   Parameters: infile -- main input file containing option for input
**                         reading and any other needed parameters (but
**                         not mesh point coords and connectivity data;
**                         if needed, these are in separate files
**   Notes: needs to find 0, 1 or 2 under the heading of "OPTREADINPUT"
**               in infile.
**   Modifications: 
**    - 4/98 GT added MakeMeshFromPoints function
**    - Added 7/98 - will read in layering information from a Child output
**         file if OPTREADLAYER is set to 1 (NMG).
**    - new variable mSearchOriginTriPtr initialized, first to zero,
**      then to center of domain. GT 1/2000  
**
\**************************************************************************/
template< class tSubNode >
tMesh< tSubNode >::
tMesh( tInputFile &infile )
  :
  nnodes(0),
  nedges(0),
  ntri(0),
  nodeList(),
  miNextNodeID(0),
  miNextEdgID(0),
  miNextTriID(0),   
  seed(0),
  layerflag(FALSE),
  mSearchOriginTriPtr(0)
{
   // mSearchOriginTriPtr:
   // initially set search origin (tTriangle*) to zero:
   // in initialisation list

   // As "layerflag" is used in this constructor, we compute it now.
   {
     int help;
     help = infile.ReadItem( help, "OPTINTERPLAYER" );
     if(help>0) layerflag=TRUE;
     else layerflag=FALSE;
   }
   // option for reading/generating initial mesh
   int read;
   read = infile.ReadItem( read, "OPTREADINPUT" );
   switch (read){
   case 0:
     MakeMeshFromScratch( infile ); //create new mesh with parameters
     break;
   case 1:
     {
       int lay;  // option for reading layer info
       MakeMeshFromInputData( infile ); //create mesh by reading data files
       lay = infile.ReadItem( lay, "OPTREADLAYER" );
       if( lay == 1 )
	 MakeLayersFromInputData( infile );
     }
     break;
   case 2:
     MakeMeshFromPoints( infile );  //create new mesh from list of points
     break;
   case 3:
     MakeRandomPointsFromArcGrid( infile ); //create mesh from regular grid
     break;
   case 4:
     MakeHexMeshFromArcGrid( infile );
     break;
   default:
     cerr << "Valid options for reading mesh input are:\n"
	  << "  0 -- create rectangular offset mesh\n"
	  << "  1 -- read mesh from input data files\n"
	  << "  2 -- create mesh from a list of (x,y,z,b) points\n"
	  << "  3 -- create random mesh from ArcGrid ascii output\n"
	  << "  4 -- create hexagonal mesh from ArcGrid ascii output\n";
     ReportFatalError( "Invalid mesh input option requested." );
     break;
   }

   // find geometric center of domain:
   double cx = 0.0;
   double cy = 0.0;
   double sumarea = 0.0;
   double carea;
   tMeshListIter< tSubNode > nI( getNodeList() );
   tNode* cn;
   for( cn = nI.FirstP(); !nI.AtEnd(); cn = nI.NextP() )
   {
      carea = cn->getVArea();
      cx += cn->getX() * carea;
      cy += cn->getY() * carea;
      sumarea += carea;
   }
   assert( sumarea>0.0 );
   cx /= sumarea;
   cy /= sumarea;
   // find triangle in which these coordinates lie and designate it the
   // search origin:
   mSearchOriginTriPtr = LocateTriangle( cx, cy );

}

//destructor
template< class tSubNode >
tMesh< tSubNode >::
~tMesh() {cout << "    ~tMesh()" << endl;}                    //tMesh


/************************************************************************\
**
**  tMesh::MakeLayersFromInputData( infile ):
**
**  This function reads in layer information from a CHILD output file and
**  copies it to the nodes which have already been created.
**
**  TODO: is there a way to handle this outside of tMesh, so tMesh
**  doesn't have to know anything about layering?
**
\************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakeLayersFromInputData( tInputFile &infile )
{
   int i, item, numl;
   int righttime;
   double time, intime;
   double ditem;
   char thestring[80], inname[80];
   char headerLine[kMaxNameLength];
   ifstream layerinfile;
   infile.ReadItem( thestring, "INPUTDATAFILE" );

   //cout<<"in MakeLayersFromInputData..."<<endl;
   
   strcpy( inname, thestring );
   strcat( inname, ".lay" );
   layerinfile.open(inname); /* Layer input file pointer */
   assert( layerinfile.good() );

   intime = infile.ReadItem( intime, "INPUTTIME" );
   //find specified input times in input data files and read no. items.
   //nodes:
   righttime = 0;
   while( !( layerinfile.eof() ) && !righttime )
   {
      layerinfile.getline( headerLine, kMaxNameLength );
      if( headerLine[0] == kTimeLineMark )
      {
         layerinfile.seekg( -layerinfile.gcount(), ios::cur );
         layerinfile >> time;
         //cout << "from file, time = " << time << endl;
         if( time == intime ) righttime = 1;
      }
   }
   if( !( layerinfile.eof() ) ) layerinfile >> nnodes;
   else
   {
      cerr << "Couldn't find specified input time in layer file" << endl;
      ReportFatalError( "Input error" );
   }

   tLayer layhelp;
   int numg;
   numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   layhelp.setDgradesize(numg);

   int g;
   tLNode * cn;
   //int nActNodes = getNodeList()->getActiveSize();
   //int NNodes = getNodeList()->get
   tMeshListIter<tLNode> ni( getNodeList() );

   for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
   {
      layerinfile >> numl;
      for(i = 1; i<=numl; i++){
         layerinfile >> ditem;
         layhelp.setCtime(ditem);
         layerinfile >> ditem;
         layhelp.setRtime(ditem);
         layerinfile >> ditem;
         layhelp.setEtime(ditem);
         layerinfile >> ditem;
         layhelp.setDepth(ditem);
         layerinfile >> ditem;
         layhelp.setErody(ditem);
         layerinfile >> item;
         layhelp.setSed(item);
         for(g=0; g<numg; g++){
            layerinfile >> ditem;
            layhelp.setDgrade(g, ditem);
         }
         cn->InsertLayerBack( layhelp );
      }
          
   }

   tArray<double> dgradebrhelp( numg );   
   double sumbr = 0;
   i=0;
   char add='1';
   char name[20];
   double help;
   
   while ( i<numg ){
      // Reading in proportions for intital regolith and bedrock
      strcpy( name, "BRPROPORTION");
      strcat( name, &add ); 
      help = infile.ReadItem( help, name);
      dgradebrhelp[i]=help;
      sumbr += help;
      i++;
      add++;
   }

   assert(sumbr>0.999 & sumbr<1.001);
   
   layhelp.setCtime(0);
   layhelp.setRtime(0);
   //layhelp.setFlag(0);
   layhelp.setErody(0);
   layhelp.setSed(0);
   ditem=layhelp.getDepth();
   for(g=0; g<numg; g++){
      layhelp.setDgrade(g, ditem*dgradebrhelp[g]);
   }   
   
   for( cn = cn; !(ni.AtEnd()); cn=ni.NextP() ){
      cn->InsertLayerBack( layhelp );
   }
   
}


/**************************************************************************\
**
**   tMesh::MakeMeshFromInputData
**
**   Formerly tMesh( tListInputData &, tlnodflag ). Constructs
**   tListInputData obj., makes mesh from data in that object.
**                    
**
**   Created: 2/11/98, SL
**   Calls: tListInputData( infile ), UpdateMesh(), CheckMeshConsistency()
**   Inputs: infile -- main input file from which various items are read
**                     (such as the names of the files containing mesh
**                     data)
**
**   Modifications:
**    - 2nd edge iterator used in CCW-setup loop to enhance speed. GT 8/98
**    - 2/02 Fixed bug in which edges connecting 2 closed boundaries were
**      not being correctly flagged as "no flux" edges (GT)
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromInputData( tInputFile &infile )
{
   int i;
   tListInputData< tSubNode > input( infile );
   seed = 0;
   // set the number of nodes, edges, and triangles in the mesh
   //assert( lnodflag );
   nnodes = input.x.getSize();
   nedges = input.orgid.getSize();
   ntri = input.p0.getSize();
   //cout << "nnodes, nedges, ntri: " << nnodes << " " << nedges << " " << ntri << endl << flush;
   assert( nnodes > 0 );
   assert( nedges > 0 );
   assert( ntri > 0 );

   // Create the node list by creating a temporary node and then iteratively
   // (1) assigning it values from the input data and (2) inserting it onto
   // the back of the node list.
   cout << "Creating node list..." << flush;
   tSubNode tempnode( infile );
   int bound;
   for( i = 0; i< nnodes; i++ )
   {
      tempnode.set3DCoords( input.x[i], input.y[i], input.z[i] );
      tempnode.setID( i );
      bound = input.boundflag[i];
      assert( bound >= 0 && bound <= 2 );
      tempnode.setBoundaryFlag( bound );
      if( bound == kNonBoundary )
          nodeList.insertAtActiveBack( tempnode );
      else if( bound == kOpenBoundary )
          nodeList.insertAtBoundFront( tempnode );
      else
          nodeList.insertAtBack( tempnode );       //kClosedBoundary
      //cout << input.x[i] << input.y[i] << input.z[i]
      //     << input.boundflag[i] << endl;
      //cout << tempnode.getBoundaryFlag() << " ";
      //cout << nodeList.getLast()->getDataPtr()->getBoundaryFlag() << endl;
   }
   cout << "done.\n";
   
   // Create and initialize the edge list by creating two temporary edges
   // (which are complementary, ie share the same endpoints) and then
   // iteratively assigning values to the pair and inserting them onto the
   // back of the edgeList
   cout << "Creating edge list..." << flush;
   tMeshListIter< tSubNode > nodIter( nodeList );
   tEdge tempedge1, tempedge2;
   int obnd, dbnd;
   for( miNextEdgID = 0; miNextEdgID < nedges-1; miNextEdgID+=2 )
   {
      // Assign values: ID, origin and destination pointers
      tempedge1.setID( miNextEdgID );
      tempedge2.setID( miNextEdgID + 1 );
      //cout << input.orgid[i] << " " << input.destid[i] << endl;
      //cout << nodIter.Get( input.orgid[i] ) << " ";
      //cout << nodIter.Get( input.destid[i] ) << endl;
      nodIter.Get( input.orgid[miNextEdgID] );
          //{
      tempedge1.setOriginPtr( &(nodIter.DatRef()) );
      tempedge2.setDestinationPtr( &(nodIter.DatRef()) );
      obnd = nodIter.DatRef().getBoundaryFlag();
         //cout << nodIter.DatRef().getID() << "->";
         //}
      nodIter.Get( input.destid[miNextEdgID] );
          //{
      tempedge1.setDestinationPtr( &(nodIter.DatRef()) );
      tempedge2.setOriginPtr( &(nodIter.DatRef()) );
      dbnd = nodIter.DatRef().getBoundaryFlag();
         //cout << nodIter.DatRef().getID() << endl;
         //}

	 // set the "flowallowed" status (FALSE if either endpoint is a
         // closed boundary, or both are open boundaries) 
	 // and insert edge pair onto the list --- active
         // part of list if flow is allowed, inactive if not
         //cout << "BND: " << obnd << " " << dbnd << " " << kClosedBoundary
         //     << endl;
      if( obnd == kClosedBoundary || dbnd == kClosedBoundary
	  || (obnd==kOpenBoundary && dbnd==kOpenBoundary) )
      {
         /*cout << "setting edges " << tempedge1.getID() << " and "
              << tempedge2.getID() << " as no-flux" << endl;*/
         tempedge1.setFlowAllowed( 0 );
         tempedge2.setFlowAllowed( 0 );
         edgeList.insertAtBack( tempedge1 );
         edgeList.insertAtBack( tempedge2 );
      }
      else
      {
         /*cout << "setting edges " << tempedge1.getID() << " and "
              << tempedge2.getID() << " as OPEN" << endl;*/
         tempedge1.setFlowAllowed( 1 );
         tempedge2.setFlowAllowed( 1 );
         edgeList.insertAtActiveBack( tempedge1 );
         edgeList.insertAtActiveBack( tempedge2 );
         //cout << "EDGFA " << tempedge2.FlowAllowed() << endl;
      }
   }
   cout << "done.\n";

   //DEBUG
   cout << "JUST ADDED EDGES:\n";
   tMeshListIter< tEdge > ei( edgeList );
   tEdge * ce;

   /*Xfor( ce=ei.FirstP(); !(ei.AtEnd()); ce=ei.NextP() )
     {
       ce->TellCoords();
       cout << ce->FlowAllowed() << endl;
     }*/

   // set up the lists of edges (spokes) connected to each node
   // (GT added code to also assign the 1st edge to "edg" as an alternative
   // to spokelist implementation)
   cout << "setting up spoke lists..." << flush;
   int e1;
   int ne;
   tMeshListIter< tEdge >
       edgIter( edgeList );
   tSubNode * curnode;
   assert( nodIter.First() );
   i = 0;
   do                                        //for( i=0; i<nnodes; i++ )
   {
      curnode = nodIter.DatPtr();
      e1 = input.edgid[curnode->getID()];  //fix of above error
      //cout << "spokes for Node " << curnode->getID() << endl;
      if( edgIter.Get( e1 ) )
      {
          curnode->insertBackSpokeList( &(edgIter.DatRef()) );
          //nodIter.DatRef().insertBackSpokeList( &(edgIter.DatRef()) );
          curnode->setEdg( edgIter.DatPtr() );
          /*cout << "Node " << curnode->getID() << " has edg "
               << (curnode->getEdg())->getID() << endl;*/
      }
      for( ne = input.nextid[e1]; ne != e1; ne = input.nextid[ne] )
      {
         if( ne>=nedges )
         {
             cerr << "Warning: edge " << e1 << " has non-existant ccw edge "
                  << ne << endl;
             cerr << "This is likely to be a problem in the edge input file"
                  << endl;
         }
         if( edgIter.Get( ne ) )
              curnode->insertBackSpokeList( &(edgIter.DatRef()) );
          //nodIter.DatRef().insertBackSpokeList( &(edgIter.DatRef()) );
      }
      i++;
   }
   while( nodIter.Next() );
   cout << "done.\n";

   // Assign ccwedg connectivity (that is, tell each edge about its neighbor
   // immediately counterclockwise) (added by gt Dec 97 to re-implement
   // previous data structure) TODO: can this be done using makeCCWEdges?
   cout << "Setting up CCW edges..." << flush;
   tEdge * curedg, * ccwedg;
   int ccwedgid;
   tMeshListIter<tEdge> ccwIter( edgeList ); // 2nd iter for performance
   //X for( i=0; i<nedges; i++ )
   for( i=0, curedg=edgIter.FirstP(); i<nedges; i++, curedg=edgIter.NextP() )
   {
      //X curedg = edgIter.GetP( i );
      ccwedgid = input.nextid[i];
      //X ccwedg = edgIter.GetP( ccwedgid );
      ccwedg = ccwIter.GetP( ccwedgid ); //test
      curedg->setCCWEdg( ccwedg );
      /*cout << "Edg " << ccwedgid << " (" << ccwedg->getOriginPtr()->getID()
           << " " << ccwedg->getDestinationPtr()->getID() << ") is ccw from "
           << curedg->getID() << " ("
           << curedg->getOriginPtr()->getID()
           << " " << curedg->getDestinationPtr()->getID() << ") " << endl;*/
      
   }
   cout << "done.\n";
   
/*   cout << "doing something else w/ spokes\n" << flush;
   tPtrListIter< tEdge > spokIter;
   assert( nodIter.First() );
   do                                        //for( i=0; i<nnodes; i++ )
   {
      spokIter.Reset( nodIter.DatRef().getSpokeListNC() );
      //cout << " node " << nodIter.DatRef().getID() << " with spoke edges";
      i = 0;
      do
      {
         if( i > 0 ) spokIter.Next();
         //cout << " " << spokIter.DatPtr()->getID();
         i++;
      }
      while( spokIter.NextIsNotFirst() );
      //cout << endl;
   }
   while( nodIter.Next() );*/

   cout << "setting up triangle connectivity..." << flush;
   for ( i=0; i<ntri; i++ )
   {
      //cout << "TRI " << i << endl << flush;
      tTriangle newtri;
      newtri.setID( i );
      if( nodIter.Get( input.p0[i] ) )
          newtri.setPPtr( 0, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p1[i] ) )
          newtri.setPPtr( 1, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p2[i] ) )
          newtri.setPPtr( 2, &(nodIter.DatRef()) );
      if( edgIter.Get( input.e0[i] ) )
          newtri.setEPtr( 0, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e1[i] ) )
          newtri.setEPtr( 1, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e2[i] ) )
          newtri.setEPtr( 2, &(edgIter.DatRef()) );
      triList.insertAtBack( newtri );
   }
   
   tListIter< tTriangle >
       triIter( triList ), triIter2( triList );
   tTriangle * ct, * nbrtri;
   for( i=0, ct=triIter.FirstP(); i<ntri; ct=triIter.NextP(), i++ )
   {
      nbrtri = ( input.t0[i]>=0 ) ? triIter2.GetP( input.t0[i] ) : 0;
      ct->setTPtr( 0, nbrtri );
      nbrtri = ( input.t1[i]>=0 ) ? triIter2.GetP( input.t1[i] ) : 0;
      ct->setTPtr( 1, nbrtri );
      nbrtri = ( input.t2[i]>=0 ) ? triIter2.GetP( input.t2[i] ) : 0;
      ct->setTPtr( 2, nbrtri );
   }
   
   cout<<"done.\n";

   // The user may wish to densify the starting mesh uniformly by adding
   // a new node at the circumcenter of every triangle. That option is
   // implemented here. We simply iterate through the list of triangles
   // and add a node at the circumcenter of each. In doing so we take
   // advantage of the fact that the circumcenter is also a Voronoi vertex,
   // and is pointed to by each of the clockwise-directed edges in the
   // triangle. The z value of each new node is obtained by linear (plane)
   // interpolation from the 3 triangle vertices.
   //   "initMeshDensLevel" serves as both a flag indicating whether the
   // user wants densification, and as an indicator of the number of passes
   // (the "level") to make -- ie, the number of times we sweep through
   // adding a node in each triangle.
   //   Added Sept. 2000, GT
   int initMeshDensLevel;
   initMeshDensLevel = infile.ReadItem( initMeshDensLevel, "OPTINITMESHDENS" );
   if( initMeshDensLevel)
   {
     int j;  // Level counter
     int nnewpoints;  // No. of new points added in a given pass
     tArray<double> newx, newy, newz;   // Lists of new coords
     tArray<double> zvals(3);   // z values of a triangle's 3 nodes
     tempnode.setBoundaryFlag( 0 );  // assumed all interior points
     for( j=1; j<=initMeshDensLevel; j++ )
     {
       // Set up for this pass
       cout << "Densifying initial mesh (level " << j << ")\n";
       UpdateMesh();
       nnewpoints = ntri = triList.getSize();  // no. of triangles in the list
       newx.setSize( nnewpoints );
       newy.setSize( nnewpoints );
       newz.setSize( nnewpoints );

       // Compute and store the x,y,z coordinates of the points to be added
       ct = triIter.FirstP();     // start with the first triangle
       for( i=0; i<ntri; i++ )    // loop through the triangles
       {
	 assert( ct!=0 );
	 tArray<double> xy = ct->ePtr(0)->getRVtx();  // get the coords
	 newx[i] = xy[0];
	 newy[i] = xy[1];

	 // Now find the z coordinate using interpolation
	 zvals[0] = ct->pPtr(0)->getZ();
	 zvals[1] = ct->pPtr(1)->getZ();
	 zvals[2] = ct->pPtr(2)->getZ();
	 newz[i] = PlaneFit( xy[0], xy[1], ct->pPtr(0)->get2DCoords(), 
			       ct->pPtr(1)->get2DCoords(), 
			       ct->pPtr(2)->get2DCoords(), zvals );
	 ct = triIter.NextP();
       }

       // Now loop through and add the nodes
       for( i=0; i<nnewpoints; i++ )
       {
	 tempnode.set3DCoords( newx[i], newy[i], newz[i] );  // assign them
	 tempnode.setID( nnodes+i );
	 AddNode( tempnode );        // Add the new node
       }
     }  // end of current densification level
   } // end of optional mesh densification  

   /*tMeshListIter< tEdge > ei( edgeList );
   tEdge * ce;*/
   cout << "JUST BEFORE UPDATEMESH\n";
   for( ce=ei.FirstP(); !(ei.AtEnd()); ce=ei.NextP() )
     {
       ce->TellCoords();
       cout << ce->getVEdgLen() << " " << ce->getBoundaryFlag() << endl;
     }



   UpdateMesh();
   CheckMeshConsistency();

   cout << "end of tMesh( input )" << endl << flush;
}


/**************************************************************************\
**
**   Note: This function works but is not used because it turned out to be
**     slower than using AddNode. It needs an effective region limitation
**     on checked points.
**
**   BatchAddNodes():
**      A nicer way (who wants to code it? oo-oo! stephen does, stephen does!):
**        1) make boundary nodes and edges between them;
**        2) add all other nodes at once;
**        3) do triangulation once by making a temporary list of boundary
**           edges and, for each edge, finding the node which would make
**           the biggest angle with the endpoints of the edge; remove edge(s)
**           from temp. list, make new edge(s) and triangle; add new edge(s)
**           to temp. list. (from Du)
**
**  Algorithm in more detail (modified from Du, C., 1996. An algorithm for
**  automatic Delaunay triangulation of arbitrary planar domains, _Adv. in
**  Eng. Software, v. 2, p. 21-26.):
**     (1) add existing edges to boundary list in CCW order, and ONLY the edges
**         pointing from node to next CCW node in boundary (complement not
**         added);
**     (2) add all nodes to temporary list for triangulation;
**     (3) at each iteration, remove first edge from boundary list;
**     (4) find 1st candidate node for triangulation with endpoints of
**         latter boundary edge such that candidate node is (a) not one of
**         other two, (b) CCW with other two, and (c) edges connecting it
**         to other two would not intersect any existing edges (only perform
**         this last check if node passed TriPasses(...) test in (5) because
**         the intersection check is most costly);
**     (5) check 1st candidate against others that meet the above conditions
**         to find node that forms the largest possible angle between the
**         boundary edge endpoints (in TriPasses(...));
**     (6) add new edges to node found in (5) if necessary both to mesh
**         and to boundary edge list (for latter, only those directed in
**         CCW direction around current working boundary);
**     (7) remove any already existing edge to node found in (5) from
**         boundary list;
**     (8) make new triangle from boundary edge endpoints and node found
**         in (5);
**     (9) remove nodes of triangle formed in (8) that are not endpoints
**           of any boundary edge;
**    (10) when boundary edge list is empty, we're done!
**                    
**
**   Created: 10/98, SL
**
**   Modification, SL, 10/98: Instead of calling
**   tMesh::IntersectsAnyEdge(), I made a (global) function that
**   takes pointers to the edge in question and the working boundary list
**   as arguments and checks only for intersection of the edge in question
**   with edges in the working boundary list. Du [1996] indicates that
**   such a limited check is sufficient as long as the list of possible
**   connecting nodes is correct, i.e., on the working boundary or contained
**   within the space defined by the working boundary.            
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
BatchAddNodes()
{
   int i;
   double mincosine, testcosine;
   tPtrList< tSubNode > tmpnodList, tmpnbrList;
   tPtrList< tEdge > tmpbndList;
   tMeshListIter< tSubNode > nI( nodeList );
   tMeshListIter< tEdge > eI( edgeList );
   tListIter< tTriangle > tI( triList );
   tPtrListIter< tEdge > bI( tmpbndList );
   tPtrListIter< tSubNode > tnI( tmpnodList ), nbrI;
   tEdge* ce;
   tEdge* be = 0;
   tEdge tempedge;
   tSubNode* n0, * n1, * n2, * ntest, * cn;
   tTriangle* ct;
   tArray< double > p0, p1, p2, ptest;

   tempedge.setID( -4 );

   // use fact that AddEdge( node0, node1, node2 ), where nodes are listed
   // in ccw order, adds two edges, and the first of the pair is always
   // directed from node0 to node1; thus, if we started with a hull
   // of nodes and added edges around the hull in ccw order, the 0th, 2nd,
   // 4th, ..., edges will be directed such that their left sides are the
   // inside.
   // construct tmpbndList by putting only the 0th, 2nd, etc., edges in
   // list:
   for( ce = eI.FirstP(); !( eI.AtEnd() ); eI.Next(), ce = eI.NextP() )
       tmpbndList.insertAtBack( ce );
   tmpbndList.makeCircular();
   //cout << tmpbndList.getSize() << " edges in initial boundary list\n";
   // DEBUG:
   for( ce = bI.FirstP(); !( bI.AtEnd() ); ce = bI.NextP() )
   {
      // assert boundary edges are arranged head to tail:
      assert( ce->getDestinationPtr() == bI.ReportNextP()->getOriginPtr() &&
              ce->getOriginPtr() == bI.ReportPrevP()->getDestinationPtr() );
      //cout << "edge " << ce->getID() << " from "
      //     << ce->getOriginPtr()->getID() << " to "
      //     << ce->getDestinationPtr()->getID() << endl;
   }
   // END DEBUG
   // put all nodes in tmpnodList:
   for( cn = nI.FirstP(); !( nI.AtEnd() ); cn = nI.NextP() )
       tmpnodList.insertAtBack( cn );
   while( tmpbndList.removeFromFront( ce ) ) // current bndy edge removed from list
   {      
      n0 = static_cast< tSubNode* >(ce->getOriginPtrNC());
      n1 = static_cast< tSubNode* >(ce->getDestinationPtrNC());
      //cout << "bndy edge " << ce->getID() << ", endpts at "
      //     << n0->getX() << ", " << n0->getY() << ", and "
      //     << n1->getX() << ", " << n1->getY() << endl;
      if( !( tmpnodList.isEmpty() ) )
      {
         for( n2 = tnI.FirstP(); !( tnI.AtEnd() ); n2 = tnI.NextP() )
         {
            if( n2 != n0 && n2 != n1 ) // if n2 isn't one of other two
            {
               if( PointsCCW( n0->get2DCoords(), n1->get2DCoords(),
                              n2->get2DCoords() ) ) // points CCW? if so...
               {
                  if( n0->EdgToNod( n2 ) == 0 ) // already edge to test node? if not...
                  {
                     tempedge.setOriginPtr( n0 );
                     tempedge.setDestinationPtr( n2 );
                     // would new edge intersect others?
                     if( IntersectsAnyEdgeInList( &tempedge, tmpbndList ) == 0 ) 
                     {                                         
                        // if not, already other edge to test node?
                        if( n1->EdgToNod( n2 ) == 0 ) 
                        {
                           //  if not, if new edge would not intersect any others
                           tempedge.setOriginPtr( n1 );
                           tempedge.setDestinationPtr( n2 );
                           if( IntersectsAnyEdgeInList( &tempedge, tmpbndList ) == 0 )
                               break; // found n2 cand.
                        }
                        else break;
                     }
                  }
                  // already other edge to test node?
                  else if( n1->EdgToNod( n2 ) == 0 ) 
                  {
                     // if not, if new edge would not intersect any others
                     tempedge.setOriginPtr( n1 );
                     tempedge.setDestinationPtr( n2 );
                     if( IntersectsAnyEdgeInList( &tempedge, tmpbndList ) == 0 )
                         break; // found n2 cand.
                  }
                  else break; // if both edges already exist, this must be the one!
               }
            }
         }
      }
      else break; // done!
      mincosine = FindCosineAngle0_2_1( n0->get2DCoords(),
                                        n1->get2DCoords(),
                                        n2->get2DCoords() );
      // test all nodes to find "correct" triangle with edge endpts:
      for( ntest = tnI.NextP(); !( tnI.AtEnd() ); ntest = tnI.NextP() )
      {
         // check that test node isn't one of other three:
         if( ntest != n0 && ntest != n1 && ntest != n2 ) 
         {
            // check that edge endpoint nodes and test node are CCW:
            if( PointsCCW( n0->get2DCoords(),
                           n1->get2DCoords(),
                           ntest->get2DCoords() ) )
            {
               // nodes are CCW;
               // now check whether test node, ntest, is better than n2:
               testcosine = FindCosineAngle0_2_1( n0->get2DCoords(),
                                                  n1->get2DCoords(),
                                                  ntest->get2DCoords() );
               //if( !TriPasses( ntest->get2DCoords(),
               //            n0->get2DCoords(), 
               //            n1->get2DCoords(),
               //            n2->get2DCoords() ) ) 
               if( testcosine < mincosine )
               {
                  // Node ntest passed the angle test; now make sure new edges won't
                  // intersect any old edges: if an edge doesn't already exist, put
                  // a dummy edge there and see if it intersects anything; if no
                  // new edges would intersect any old ones, set n2 = ntest
                  // (i.e., choose ntest):
                  // already edge to test node? if not...
                  if( n0->EdgToNod( ntest ) == 0 ) 
                  {
                     tempedge.setOriginPtr( n0 );
                     tempedge.setDestinationPtr( ntest );
                     // would new edge intersect others?
                     if( IntersectsAnyEdgeInList( &tempedge, tmpbndList ) == 0 ) 
                     {                                         // if not...
                        // already other edge to test node? if not...
                        if( n1->EdgToNod( ntest ) == 0 ) 
                        {
                           tempedge.setOriginPtr( n1 );
                           tempedge.setDestinationPtr( ntest );
                           // and no intersection, set n2 to ntest:
                           if( IntersectsAnyEdgeInList( &tempedge, tmpbndList ) == 0 )
                           {
                              n2 = ntest;
                              mincosine = testcosine;
                           }
                        }
                        else
                        {
                           n2 = ntest; // if so...set n2 to ntest
                           mincosine = testcosine;
                        }
                     }
                  }
                  // already other edge to test node? if not...
                  else if( n1->EdgToNod( ntest ) == 0 ) 
                  {
                     tempedge.setOriginPtr( n1 );
                     tempedge.setDestinationPtr( ntest );
                     // and new edge wouldn't intersect, set n2 to ntest:
                     if( IntersectsAnyEdgeInList( &tempedge, tmpbndList ) == 0 )
                     {
                         n2 = ntest;
                         mincosine = testcosine;
                     }
                  }
                  else
                  {
                     n2 = ntest; // set n2 to ntest
                     mincosine = testcosine;
                  }
               }
            }
         }
      }
      // at end of all this, we should have the "correct" node, n2
      // 4 DEBUG
      //cout << "found node " << n2->getID() << " at "
      //     << n2->getX() << ", " << n2->getY() << ", to make tri with nodes "
      //     << n0->getID() << " and " << n1->getID() << endl;
      // END DEBUG
      tmpnbrList.Flush();
      // check for edges to "new" node
      if( n2->EdgToNod( n0 ) == 0 )
      {
         AddEdge( n2, n0, n1 ); // if edge does not exist, add it to mesh
         tmpbndList.insertAtBack( n0->EdgToNod( n2 ) ); // and to boundary list
      }
      else
      {
         if( bI.Get( n2->EdgToNod( n0 )->getID() ) ) // remove edge from boundary list
         {
            tmpbndList.moveToBack( bI.NodePtr() );
            tmpbndList.removeFromBack( be );
         }
         else cerr << "n2-n0 edge "
                   << n2->EdgToNod( n0 )->getID() << " not in temp bound list\n";
      }
      if( n1->EdgToNod( n2 ) == 0 )
      {
         AddEdge( n1, n2, n0 ); // if edge does not exist, add it to mesh
         tmpbndList.insertAtBack( n2->EdgToNod( n1 ) ); // and to boundary list
      }
      else
      {
         if( bI.Get( n1->EdgToNod( n2 )->getID() ) ) // remove edge from boundary list
         {
            tmpbndList.moveToBack( bI.NodePtr() );
            tmpbndList.removeFromBack( be );
         }
         else cerr << "n1-n2 edge "
                   << n2->EdgToNod( n0 )->getID() << " not in temp bound list\n";
      }
      // make new triangle:
      tmpnbrList.insertAtBack( n0 );
      tmpnbrList.insertAtBack( n1 );
      tmpnbrList.insertAtBack( n2 );
      tmpnbrList.makeCircular();
      nbrI.Reset( tmpnbrList );
      MakeTriangle( tmpnbrList, nbrI );
      // check nodes of new triangle (at end of triList); remove nodes not connected
      // to boundary edge from tmpnodList:
      ct = tI.LastP();
      for( i=0; i<3; i++ )
      {
         cn = static_cast<tSubNode *>(ct->pPtr( i ));
         assert( cn != 0 );
         be = 0;
         for( ce = bI.FirstP(); !( bI.AtEnd() ); ce = bI.NextP() )
             if( ce->getOriginPtr() == cn )
                 be = ce;
         if( be == 0 )
         {
            if( tnI.Get( cn->getID() ) )
            {
               tmpnodList.moveToBack( tnI.NodePtr() );
               tmpnodList.removeFromBack( cn );
            }
            else cerr << "node " << cn->getID() << " was not in temp list\n";
         }
      }  
      cout << tmpbndList.getSize() << " edges in current boundary list\n";
      cout << tmpnodList.getSize() << " nodes in current temp list\n";
      ce = bI.FirstP(); // not sure why (or whether still) this line is necessary,
                        // but it fixed a bug
   }
}


/**************************************************************************\
**
**   tMesh::MakeMeshFromScratch( infile )
**
**   Formerly tMesh( infile ). Makes
**   mesh from scratch; reads parameters
**   from input file to get mesh size, spacing, method of point
**   placement.
**
**      Could probably be done more gracefully, but here's how it does it:
**        1) makes boundary nodes and edges between them;
**        2) makes triangulation with only boundary nodes;
**        3) adds the rest of the nodes one at a time, i.e., triangulation
**           redone for each node added.
**
**      A nicer way (who wants to code it? oo-oo! stephen does, stephen does!):
**        1) make boundary nodes and edges between them;
**        2) add all other nodes at once;
**        3) do triangulation once by making a temporary list of boundary
**           edges and, for each edge, finding the node which would make
**           the biggest angle with the endpoints of the edge; remove edge(s)
**           from temp. list, make new edge(s) and triangle; add new edge(s)
**           to temp. list. (from Du)
**
**   Created: 2/11/98, SL
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: xGrid, yGrid, boundType, mElev, ptPlace, delGrid, numPts,
**               upperZ, xout, yout
**   Modified: 3/13/98--now makes rows offset so pattern is "hexagonal"
**      rather than square
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakePointBoundary( const ParamMMFS_t &Param, tInputFile &infile,
		   tPtrList< tSubNode > &bndList)
{
   int i,                        // counters 
       n;                        // no. of nodes along a side
   double dist;                  // current distance along boundary
   tSubNode tempnode( infile );  // temporary node used to create node list
   tMeshListIter< tSubNode > nodIter( nodeList );

   //MAKE BOUNDARY
   if( Param.boundType == kCornerOutlet )
   {
      miNextNodeID = 0;
      tempnode.setBoundaryFlag( kOpenBoundary );
      tempnode.set3DCoords( 0, 0, 0 );
      tempnode.setID( miNextNodeID );
      n = ROUND( Param.xGrid / Param.delGrid );
      tempnode.setBoundaryFlag( kOpenBoundary );
      nodeList.insertAtBack( tempnode );
      bndList.insertAtBack( nodIter.LastP() );
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=1, miNextNodeID++; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( Param.boundType == kOpenSide )
   {
      cout << "OPEN SIDE boundary\n";
      n = ROUND( Param.xGrid / Param.delGrid );
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=1, miNextNodeID=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>=0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   if( Param.boundType == kOppositeSidesOpen )
   {
      n = ROUND( Param.xGrid / Param.delGrid );
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=1, miNextNodeID=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<=n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kOpenBoundary );
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n-1; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, Param.upperZ );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBoundFront( tempnode );
         bndList.insertAtBack( nodIter.FirstBoundaryP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>=0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( Param.boundType == kAllSidesOpen )
   {
      miNextNodeID = 0;
      n = ROUND( Param.xGrid / Param.delGrid );
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( Param.boundType == kSpecifyOutlet )
   {
      // Create nodes for bottom (Y=0) boundary and place them on list
      n = ROUND( Param.xGrid / Param.delGrid );
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=0, miNextNodeID=0; i<n; i++, miNextNodeID++ )
      {
         // Assign node coords to tempnode and add tempnode to list
         dist = i * Param.delGrid + 0.01 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );

         // If user wants outlet on this side between this and the next pt,
         // create the outlet now
         if( Param.yout == 0 && Param.xout > dist && Param.xout < dist + Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            miNextNodeID++;
            tempnode.setID( miNextNodeID );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.FirstBoundaryP() );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }

      // Create nodes for right (X=Param.xGrid) boundary and place them on list
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( Param.xout == Param.xGrid && Param.yout > dist && Param.yout < dist + Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            miNextNodeID++;
            tempnode.setID( miNextNodeID );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.FirstBoundaryP() );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }

      // Create nodes for top (Y=Param.yGrid) boundary and place them on list
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( Param.yout == Param.yGrid && Param.xout < dist && Param.xout > dist - Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            miNextNodeID++;
            tempnode.setID( miNextNodeID );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.FirstBoundaryP() );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }

      // Create nodes for left (X=0) boundary and place them on list
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( Param.xout == 0 && Param.yout < dist && Param.yout > dist - Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            miNextNodeID++;
            tempnode.setID( miNextNodeID );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.FirstBoundaryP() );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
   }
}

template< class tSubNode >
void tMesh< tSubNode >::
MakePointInterior( const ParamMMFS_t &Param, tInputFile &infile )
{
   int i, j,                     // counters 
       nx, ny;                   // no. of nodes along a side
   double slope;
   tArray< double > xyz(3);
   tSubNode tempnode( infile );  // temporary node used to create node list

   // Add the interior points. 
   // Variations on the theme are these:
   // 1 - If points are uniform, set up a staggered mesh (alternate rows
   //     are offset by + or - 25% of Param.delGrid, respectively)
   // 2 - If points are "perturbed", do the same but add an extra random offset
   //     of +/- 25% of Param.delGrid
   // 3 - If top and bottom sides (but not left & right) are open boundaries,
   //     modify elevations to set up a slope from top to bottom
   // 4 - If points are "random", simply pick Param.numPts random locations within
   //     the interior
   tempnode.setBoundaryFlag( kNonBoundary );
   if( Param.ptPlace == kUniformMesh || Param.ptPlace == kPerturbedMesh )
   {
      nx = ROUND( Param.xGrid / Param.delGrid );  // no. points in x direction
      ny = ROUND( Param.yGrid / Param.delGrid );  // no. points in y direction
      for( i=1; i<nx; i++ )
      {
         for( j=1; j<ny; j++, miNextNodeID++ )
         {
            //rows are offset such that there should be an
            //edge leading to a corner outlet -- NB: no longer true!
            // Random offset amplified to 25%!
            // TODO: ensure integrity of corner outlet!
            xyz[0] = i * Param.delGrid - 0.25 * Param.delGrid * (j%2)
                + 0.25 * Param.delGrid * ((j+1)%2);
            xyz[1] = j * Param.delGrid;
            if( Param.ptPlace == kPerturbedMesh )
            {
               xyz[0] += 0.5 * Param.delGrid * ( ran3( &seed ) - 0.5 );
               xyz[1] += 0.5 * Param.delGrid * ( ran3( &seed ) - 0.5 );
            }
            xyz[2] = Param.mElev + Param.mElev * ( ran3( &seed ) - 0.5 );
            if( Param.boundType == kOppositeSidesOpen || Param.kSloped)
            {
               slope = Param.upperZ / Param.yGrid;
               xyz[2] += slope * xyz[1] - Param.mElev;
            }
            tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
            tempnode.setID( miNextNodeID );
            AddNode( tempnode );
         }
      }
   }
   else if( Param.ptPlace == kRandomMesh )
   {
      for( i=0; i<Param.numPts; i++ )
      {
         // Randomize x,y, and z coordinates
         xyz[0] = ran3(&seed) * Param.xGrid;
         xyz[1] = ran3(&seed) * Param.yGrid;
         xyz[2] = Param.mElev + Param.mElev * ( ran3( &seed ) - 0.5 );
         if( xyz[0] != 0 && xyz[0] != Param.xGrid && xyz[1] != 0 && xyz[1] != Param.yGrid )
         {
            tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
            tempnode.setID( miNextNodeID );
            AddNode( tempnode );
            miNextNodeID++;
         }
      }
   }

   // If user wants a specified outlet location, place the outlet point
   // now (unless the outlet is right along one of the boundaries, in which
   // case it will already have been created during boundary setup)
   // Added by GT, 5/14/99
   // Note: potential gotcha is that we don't check to see if there's 
   // already another point at the same location. TODO (see tInlet)
   if( Param.boundType==kSpecifyOutlet && Param.xout!=0 && Param.yout!=0 )
   {
      tempnode.setBoundaryFlag( kOpenBoundary );
      tempnode.set3DCoords( Param.xout, Param.yout, 0 );
      tempnode.setID( miNextNodeID );
      miNextNodeID++;
      AddNode( tempnode );
   }
}

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromScratch( tInputFile &infile )
{
   //cout << "In MGFS, calling node constr w/ infile\n";

   tSubNode *node0, *node1, *node2;
   tPtrList< tSubNode > bndList;

   seed = infile.ReadItem( seed, "SEED" );

   // Parameters defined in Input File
   ParamMMFS_t Param(infile);

   // Make Boundary
   MakePointBoundary(Param, infile, bndList);
   bndList.makeCircular();
   cout << "made points; now adding edges\n";
   
   // Add edges
   tPtrListIter< tSubNode > bndIter( bndList );
   for( node0 = bndIter.FirstP(); !( bndIter.AtEnd() ); node0 = bndIter.NextP() )
   {
      node1 = bndIter.ReportNextP();
      node2 = bndIter.ReportPrevP();
      AddEdge( node0, node1, node2 );
   }
   nnodes = nodeList.getSize();
   nedges = edgeList.getSize();
   ntri = 0;

   //DumpEdges:
   /*cout << "edges:" << endl;
   tMeshListIter< tEdge > edgIter( edgeList );
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << endl;
           }*/

   cout << "calling repair mesh for initial boundary\n";
   int meshok = RepairMesh( bndList );
   assert( meshok );

   // Add the interior points.
   cout << "filling in points\n";
   MakePointInterior(Param, infile);

   // Now finalize the initialization by updating mesh properties
   MakeCCWEdges();
   UpdateMesh(); //calls CheckMeshConsistency()  TODO: once bug-free,
   CheckMeshConsistency();                     //remove CMC call from UM
}

/**************************************************************************\
**
**   tMesh::MakeMeshFromPoints
**
**   Constructs a mesh from a given set of (x,y,z,b) points, where
**   b=boundary code. Unlike MakeMeshFromInputData no connectivity
**   information is needed, just the coordinates and boundary codes.
**
**   The format of the file containing points is:
**        NP
**        X0 Y0 Z0 B0
**        X1 Y1 Z1 B1 ...etc.
**   where NP is the number of points in the file, X, Y, and Z are
**   x, y, and z coords, and B is the boundary code.
**
**   Variations: to reduce memory overhead, the routine could be modified
**   to read one (x,y,z,b) set at a time rather than reading and
**   temporarily storing all the values. Also, instead of giving a
**   boundary code for every point, a separate list of open boundary
**   could be read (eg, from a separate file or from the main input file)
**
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: infile -- main parameter input file
**   Assumes: infile is valid and open
**   Created: 4/98 GT
**   Modified: 
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromPoints( tInputFile &infile )
{
   int i;                           // loop counter
   int numpts;                      // no. of points in mesh
   tArray<double> x, y, z;          // arrays of x, y, and z coordinates
   tArray<int> bnd;                 // array of boundary codes 
   char pointFilenm[80];            // name of file containing (x,y,z,b) data
   ifstream pointfile;              // the file (stream) itself
   double minx = 1e12, miny = 1e12, // minimum x and y coords
       maxx = 0, maxy=0,            // maximum x and y coords 
       dx, dy;                      // max width and height of region
   tSubNode tempnode( infile ),     // temporary node used in creating new pts
       *stp1, *stp2, *stp3;         // supertriangle vertices

   cout<<"MakeMeshFromPoints"<<endl<<flush;
   
   // get the name of the file containing (x,y,z,b) data, open it,
   // and read the data into 4 temporary arrays
   infile.ReadItem( pointFilenm, "POINTFILENAME" );
   pointfile.open( pointFilenm );
   if( !pointfile.good() )
   {
      cerr << "Point file name: '" << pointFilenm << "'\n";
      ReportFatalError( "I can't find a file by this name." );
   }
   pointfile >> numpts;
   x.setSize( numpts );
   y.setSize( numpts );
   z.setSize( numpts );
   bnd.setSize( numpts );
   for( i=0; i<numpts; i++ )
   {
      if( pointfile.eof() )
          ReportFatalError( "Reached end-of-file while reading points." );
      pointfile >> x[i] >> y[i] >> z[i] >> bnd[i];
      //if( bnd[i]<0 || bnd[i]>2 )
      //    ReportWarning( "Invalid boundary code." );
      if( x[i]<minx ) minx = x[i];
      if( x[i]>maxx ) maxx = x[i];
      if( y[i]<miny ) miny = y[i];
      if( y[i]>maxy ) maxy = y[i];
      
   }
   pointfile.close();
   cout << "finished reading in points"<< endl;
   dx = maxx - minx;
   dy = maxy - miny;

   // Create the 3 nodes that form the supertriangle and place them on the
   // node list in counter-clockwise order. (Note that the base and height
   // of the supertriangle are 7 times the
   // width and height, respectively, of the rectangle that encloses the
   // points.) Assigning the IDs allows us to retrieve and delete these
   // nodes when we're done creating the mesh.
   cout << "creating supertri: min & max are (" << minx << "," << miny << ") (" << maxx << "," << maxy << ")\n";
   
   tempnode.set3DCoords( minx-3*dx, miny-3*dy, 0.0 );
   tempnode.setBoundaryFlag( kClosedBoundary );
   tempnode.setID( -1 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( maxx+3*dx, miny-3*dy, 0.0 );
   tempnode.setID( -2 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( minx+0.5*dx, maxy+3*dy, 0.0 );
   tempnode.setID( -3 );
   nodeList.insertAtBack( tempnode );

   cout << "Supertri coords: " << minx-3*dx << "," << miny-3*dy << "  " << maxx+3*dx << "," << miny-3*dy << "  " << minx+0.5*dx << "," << maxy+3*dy << endl;

   // set # of nodes, edges, and triangles
   nnodes = 3;
   nedges = ntri = 0;

   // Create the edges that connect the supertriangle vertices and place
   // them on the edge list.
   // (To do this we need to retrieve pointers from the nodeList)
   tMeshListIter<tSubNode> nodIter( nodeList );
   stp1 = nodIter.FirstP();
   stp2 = nodIter.NextP();
   stp3 = nodIter.NextP();
   AddEdge( stp1, stp2, stp3 );  // edges 1->2 and 2->1
   AddEdge( stp2, stp3, stp1 );  // edges 2->3 and 3->2
   AddEdge( stp3, stp1, stp2 );  // edges 3->1 and 1->3

   // set up the triangle itself and place it on the list. To do this, we
   // just set up a list of pointers to the three nodes in the super tri
   // and pass the list (along with an iterator) to MakeTriangle.
   tPtrList<tSubNode> supertriptlist;
   supertriptlist.insertAtBack( stp1 );
   supertriptlist.insertAtBack( stp2 );
   supertriptlist.insertAtBack( stp3 );
   supertriptlist.makeCircular();
   tPtrListIter<tSubNode> stpIter( supertriptlist );
   MakeTriangle( supertriptlist, stpIter );

   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: " << nedges << " NT: " << ntri << endl << flush;
   cout << "c4\n";

   // Now add the points one by one to construct the mesh.
   for( i=0; i<numpts; i++ )
   {
      //cout << "IN MGFP c0, ADDING NODE " << i << endl;
      //Xtempnode.setID( miNextNodeID );
      tempnode.set3DCoords( x[i],y[i],z[i] );
      tempnode.setBoundaryFlag( bnd[i] );
      //if(bnd[i]==kNonBoundary && z[i]<0)
      //  cout<<"problem at x "<<x[i]<<" y "<<y[i]<<endl<<flush;
      AddNode( tempnode );
   }

   cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: " << nedges << " NT: " << ntri << endl << flush;

   // We no longer need the supertriangle, so remove it by deleting its
   // vertices.
   DeleteNode( stp1, kNoRepair );
   DeleteNode( stp2, kNoRepair );
   DeleteNode( stp3, kNoRepair );
   
   cout << "3 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: " << nedges << " NT: " << ntri << endl << flush;
   
   // Update Voronoi areas, edge lengths, etc., and test the consistency
   // of the new mesh.
   UpdateMesh();
   CheckMeshConsistency( );

   // Clean up (probably not strictly necessary bcs destructors do the work)
   supertriptlist.Flush();

   //XPtrListDebug::TellAll();
   
}

/**************************************************************************\
**
**   tMesh::MakeRandomPointsFromArcGrid
**
**		Routine to make irregular mesh from regular Arc grid by randomly 
**		sampling the grid such that the average resolution of the irregular
**    mesh is equal to that of the Arc grid.
**
**    Not quite random: before a node is added, it is checked for its
**    proximity to existing nodes. If it is as close or closer than 1/10th 
**    the nominal discretization to an existing node, it is not added.
**
**		Designed to read from a grid containing points either within an
**		already isolated basin or containing "no data".
**
**
**
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: infile -- main parameter input file
**   Assumes: infile is valid and open
**   Created: 10/98 SL
**   Modified: 
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakeRandomPointsFromArcGrid( tInputFile &infile )
{
   int i, j;                        // loop counter
   //Xn;                            // counter
   //int updatemesh;                  // flag to indicate whether AddNode
                                    // should call UpdateMesh
   int numrows, numcols, numpts;    // no. of rows, cols, points in mesh
   int delgrid,                     // arc grid cell size (side, meters) 
       nodata;                      // value signifying and no. w/ "no data"
   int numActNbrs;
   double x, y;                     // x, y coordinates
   tArray<int> bnd;                 // array of boundary codes 
   char arcgridFilenm[80];          // name of file containing arc grid data
   ifstream gridfile;               // the file (stream) itself
   double minx = 1e12, miny = 1e12, // minimum x and y coords
       maxx = 0, maxy=0,            // maximum x and y coords 
       minz= 1e12, /*maxz = 0,*/    // min. and max. elevs.
       //minzx, minzy,            // x and y coords of min. elevation
       dx, dy,                      // max width and height of region (meters)
       di, dj,                      // width and height of region in pixels
       xgen, ygen,                  // randomly generated x and y val's
       zinterp;                     // interp'd elev.
   double mindist;
   double delx, dely, dist;
   //XtSubNode *closestPtr;
   tSubNode tempnode( infile ),     // temporary node used in creating new pts
       *stp1, *stp2, *stp3;         // supertriangle vertices
   //Xdouble dumval;
   char dumhead[3];
   tSubNode *cn, *minzPtr(0);
   tEdge* ce;
   tMeshListIter< tSubNode > nI( nodeList );
   tPtrList<tSubNode> supertriptlist, deletelist;
   tPtrListIter<tSubNode> stpIter( supertriptlist ), dI( deletelist );
   tPtrListIter< tEdge > sI;
 
   // get the name of the file containing (x,y,z,b) data, open it,
   // and read the data into 4 temporary arrays
   infile.ReadItem( arcgridFilenm, "ARCGRIDFILENAME" );
   gridfile.open( arcgridFilenm );
   if( !gridfile.good() )
   {
      cerr << "Arc grid file name: '" << arcgridFilenm << "'\n";
      ReportFatalError( "I can't find a file by this name." );
   }
   gridfile >> dumhead >> numcols >> dumhead >> numrows >> dumhead 
            >> minx >> dumhead >> miny >> dumhead >> delgrid 
            >> dumhead >> nodata;
   cout << "Arc grid with: " << numcols << " cols; " << numrows
        << " rows; LL x " << minx << "; LL y " << miny
        << "; grid spacing (m) " << delgrid << "; nodata value "
        << nodata << endl;
   
   maxx = minx + ( numcols - 1 ) * delgrid;
   maxy = miny + ( numrows - 1 ) * delgrid;
   dx = maxx - minx;
   dy = maxy - miny;
   di = numcols;
   dj = numrows;
   // create matrix from input file:
   tMatrix< double > elev( numrows, numcols );
   for( j=0; j<numrows; j++ )
   {
      for( i=0; i<numcols; i++ )
      {
         if( gridfile.eof() )
             ReportFatalError( "Reached end-of-file while reading points." );
         gridfile >> elev( j, i );
      }
   }
   gridfile.close();
   cout << "finished reading file," << gridfile << endl;
   // Create the 3 nodes that form the supertriangle and place them on the
   // node list in counter-clockwise order. (Note that the base and height
   // of the supertriangle are 5 times the
   // width and height, respectively, of the rectangle that encloses the
   // points.) Assigning the IDs allows us to retrieve and delete these
   // nodes when we're done creating the mesh.
   tempnode.set3DCoords( minx-2*dx, miny-2*dy, 0.0 );
   tempnode.setBoundaryFlag( kClosedBoundary );
   tempnode.setID( -3 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( maxx+2*dx, miny-2*dy, 0.0 );
   tempnode.setID( -2 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( minx+0.5*dx, maxy+2*dy, 0.0 );
   tempnode.setID( -1 );
   nodeList.insertAtBack( tempnode );

   // set # of nodes, edges, and triangles
   nnodes = 3;
   nedges = ntri = 0;

   // Create the edges that connect the supertriangle vertices and place
   // them on the edge list.
   // (To do this we need to retrieve pointers from the nodeList)
   stp1 = nI.FirstP();
   stp2 = nI.NextP();
   stp3 = nI.NextP();
   AddEdge( stp1, stp2, stp3 );  // edges 1->2 and 2->1
   AddEdge( stp2, stp3, stp1 );  // edges 2->3 and 3->2
   AddEdge( stp3, stp1, stp2 );  // edges 3->1 and 1->3

   // set up the triangle itself and place it on the list. To do this, we
   // just set up a list of pointers to the three nodes in the super tri
   // and pass the list (along with an iterator) to MakeTriangle.
   supertriptlist.insertAtBack( stp1 );
   supertriptlist.insertAtBack( stp2 );
   supertriptlist.insertAtBack( stp3 );
   supertriptlist.makeCircular();
   MakeTriangle( supertriptlist, stpIter );
   cout << "formed supertriangle\n";
   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: " << nedges << " NT: " << ntri << endl;

   cout << "begin interpolation\n";
   // Read and initialize seed for random number generation
   seed = infile.ReadItem( seed, "SEED" );
   srand48( seed );
   numpts = numcols * numrows;
   tempnode.setBoundaryFlag( kNonBoundary );
   //Xn = 0;
   mindist = delgrid / 10.0;
   for( i=0; i<numpts; ++i )
   {
      xgen = drand48() * (di - 1.0);
      ygen = drand48() * (dj - 1.0);

      zinterp = InterpSquareGrid( xgen, ygen, elev, nodata );
      
      // reset to "real" coords and add node:
      if( zinterp != nodata )
          tempnode.setBoundaryFlag( kNonBoundary );
      else tempnode.setBoundaryFlag( kClosedBoundary );
      x = xgen * delgrid + minx;
      y = ygen * delgrid + miny;
      tempnode.set3DCoords( x, y, zinterp );
      // check whether node is closer than delgrid/10 to another node:
      cn = nI.FirstP();
      dist = delgrid;
      while( dist > mindist && !( nI.AtEnd() ) )
      {
         delx = cn->getX() - tempnode.getX();
         dely = cn->getY() - tempnode.getY();
         dist = sqrt( delx * delx + dely * dely );
         cn = nI.NextP();
      }
      // if not too close, add it:
      if( dist > mindist )
      {
         cn = AddNode( tempnode, /*updatemesh =*/ 0 );
         if( zinterp != nodata && zinterp < minz )
         {
            minz = zinterp;
            minzPtr = cn;
         }
      }
      else --i; // otherwise, decrement i and try again
   }      
   cout << "finished interpolation:";
   cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: "
        << nedges << " NT: " << ntri << endl;

   // make lowest node outlet (open boundary)
   assert( minzPtr != 0 );
   minzPtr->setBoundaryFlag( kOpenBoundary );
   nI.Get( minzPtr->getID() );
   nodeList.moveToBoundFront( nI.NodePtr() );   
   cout << "created open boundary outlet: " << nI.FirstBoundaryP()->getID() << "\n";
   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: "
        << nedges << " NT: " << ntri << endl;

   // remove closed boundary nodes that do not have non-boundary nbrs:
   for( cn = nI.FirstBoundaryP(); !( nI.AtEnd() ); cn = nI.NextP() )
   {
      numActNbrs = 0;
      sI.Reset( cn->getSpokeListNC() );
      for( ce = sI.FirstP(); !( sI.AtEnd() ); ce = sI.NextP() )
          if( ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary )
              ++numActNbrs;
      if( numActNbrs ) cn->setZ( minz );
      else deletelist.insertAtBack( cn ); 
   }
   for( cn = dI.FirstP(); !( dI.AtEnd() ); cn = dI.FirstP() )
   {
      DeleteNode( cn, kNoRepair );
      deletelist.removeFromFront( cn );
   }
   cout << "deleted superfluous boundary nodes\n";
   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: "
        << nedges << " NT: " << ntri << endl;

   // Update Voronoi areas, edge lengths, etc., and test the consistency
   // of the new mesh.
   UpdateMesh();
}


/**************************************************************************\
**
**   tMesh::MakeHexMeshFromArcGrid
**
**		Routine to make irregular mesh from regular Arc grid by randomly 
**		sampling the grid such that the average resolution of the irregular
**    mesh is equal to that of the Arc grid.
**
**		Designed to read from a grid containing points either within an
**		already isolated basin or containing "no data".
**
**   IMPORTANT: Designed for use with an isolated subbasin--only makes
**     one outlet (open boundary) node
**
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: infile -- main parameter input file
**   Assumes: infile is valid and open
**   Created: 9/98 SL
**   Modified: 
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakeHexMeshFromArcGrid( tInputFile &infile )
{
   /*bool*/int keepgoing;
   int i, j;                        // loop counter
   //Xn;                            // counter
   //int updatemesh;                  // flag to indicate whether AddNode
                                    // should call UpdateMesh
   int numrows, numcols;            // no. of rows, cols, points in mesh
   int delgrid,                     // arc grid cell size (side, meters) 
       nodata;                      // value signifying and no. w/ "no data"
   int numActNbrs;
   double x, y;                     // x, y coordinates
   tArray<int> bnd;                 // array of boundary codes 
   char arcgridFilenm[80];          // name of file containing arc grid data
   ifstream gridfile;               // the file (stream) itself
   double minx = 1e12, miny = 1e12, // minimum x and y coords
       maxx = 0, maxy=0,            // maximum x and y coords 
       minz= 1e12, /*maxz = 0,*/        // min. and max. elevs.
       /*minzx, minzy,*/                // x and y coords of min. elevation
       dx, dy,                      // max width and height of region (meters)
       di, dj,                      // width and height of region in pixels
       xgen, ygen,                  // randomly generated x and y val's
       zinterp;                     // interp'd elev.
   tSubNode tempnode( infile ),     // temporary node used in creating new pts
       *stp1, *stp2, *stp3;         // supertriangle vertices
   //Xdouble dumval;
   char dumhead[3];
   tSubNode *cn, *minzPtr(0);
   tEdge* ce;
   tMeshListIter< tSubNode > nI( nodeList );
   tPtrList<tSubNode> supertriptlist, deletelist;
   tPtrListIter<tSubNode> stpIter( supertriptlist ), dI( deletelist );
   tPtrListIter< tEdge > sI;
 
   // get the name of the file containing (x,y,z,b) data, open it,
   // and read the data into 4 temporary arrays
   infile.ReadItem( arcgridFilenm, "ARCGRIDFILENAME" );
   gridfile.open( arcgridFilenm );
   if( !gridfile.good() )
   {
      cerr << "Arc grid file name: '" << arcgridFilenm << "'\n";
      ReportFatalError( "I can't find a file by this name." );
   }
   gridfile >> dumhead >> numcols >> dumhead >> numrows >> dumhead 
            >> minx >> dumhead >> miny >> dumhead >> delgrid 
            >> dumhead >> nodata;
   cout << "Arc grid with: " << numcols << " cols; " << numrows
        << " rows; LL x " << minx << "; LL y " << miny
        << "; grid spacing (m) " << delgrid << "; nodata value "
        << nodata << endl;
   
   maxx = minx + ( numcols - 1 ) * delgrid;
   maxy = miny + ( numrows - 1 ) * delgrid;
   dx = maxx - minx;
   dy = maxy - miny;
   di = numcols;
   dj = numrows;
   // create matrix from input file:
   tMatrix< double > elev( numrows, numcols );
   for( j=0; j<numrows; j++ )
   {
      for( i=0; i<numcols; i++ )
      {
         if( gridfile.eof() )
             ReportFatalError( "Reached end-of-file while reading points." );
         gridfile >> elev( j, i );
      }
   }
   gridfile.close();
   cout << "finished reading file," << gridfile << endl;
   // Create the 3 nodes that form the supertriangle and place them on the
   // node list in counter-clockwise order. (Note that the base and height
   // of the supertriangle are 5 times the
   // width and height, respectively, of the rectangle that encloses the
   // points.) Assigning the IDs allows us to retrieve and delete these
   // nodes when we're done creating the mesh.
   tempnode.set3DCoords( minx-2*dx, miny-2*dy, 0.0 );
   tempnode.setBoundaryFlag( kClosedBoundary );
   tempnode.setID( -3 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( maxx+2*dx, miny-2*dy, 0.0 );
   tempnode.setID( -2 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( minx+0.5*dx, maxy+2*dy, 0.0 );
   tempnode.setID( -1 );
   nodeList.insertAtBack( tempnode );

   // set # of nodes, edges, and triangles
   nnodes = 3;
   nedges = ntri = 0;

   // Create the edges that connect the supertriangle vertices and place
   // them on the edge list.
   // (To do this we need to retrieve pointers from the nodeList)
   stp1 = nI.FirstP();
   stp2 = nI.NextP();
   stp3 = nI.NextP();
   AddEdge( stp1, stp2, stp3 );  // edges 1->2 and 2->1
   AddEdge( stp2, stp3, stp1 );  // edges 2->3 and 3->2
   AddEdge( stp3, stp1, stp2 );  // edges 3->1 and 1->3

   // set up the triangle itself and place it on the list. To do this, we
   // just set up a list of pointers to the three nodes in the super tri
   // and pass the list (along with an iterator) to MakeTriangle.
   supertriptlist.insertAtBack( stp1 );
   supertriptlist.insertAtBack( stp2 );
   supertriptlist.insertAtBack( stp3 );
   supertriptlist.makeCircular();
   MakeTriangle( supertriptlist, stpIter );
   cout << "formed supertriangle\n";
   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: " << nedges << " NT: " << ntri << endl;

   // place points in hexagonal mesh, i.e., equilateral triangles:
   xgen = ygen = 0.0;
   j = 0;
   miNextNodeID = 0;
   keepgoing = 1;
   while( keepgoing )
   {
      zinterp = InterpSquareGrid( xgen, ygen, elev, nodata );
      
      // reset to "real" coords and add node:
      if( zinterp != nodata )
          tempnode.setBoundaryFlag( kNonBoundary );
      else tempnode.setBoundaryFlag( kClosedBoundary );
      x = xgen * delgrid + minx;
      y = ygen * delgrid + miny;
      tempnode.setID( miNextNodeID );
      miNextNodeID++;
      tempnode.set3DCoords( x, y, zinterp );
      cn = AddNode( tempnode, /*updatemesh =*/ 0 );
      if( zinterp != nodata && zinterp < minz )
      {
         minz = zinterp;
         minzPtr = cn;
      }
      if( ygen <= dj - 1.732 )
      {
         // along x-dir:
         if( ( j%2 == 0 && xgen < di - 1.0 ) || ( j%2 == 1 && xgen < di - 2.0 ) )
             ++xgen;
         // at end of row, start odd row:
         else 
         {
            ++j;
            xgen = 0.5 * (j%2);
            ygen += 0.866; // sqrt(3)/2
         }
      }
      else keepgoing = 0;
   }
   cout << "finished interpolation:";
   cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: "
        << nedges << " NT: " << ntri << endl;

   // make lowest node outlet (open boundary)
   assert( minzPtr != 0 );
   minzPtr->setBoundaryFlag( kOpenBoundary );
   nI.Get( minzPtr->getID() );
   nodeList.moveToBoundFront( nI.NodePtr() );   
   cout << "created open boundary outlet: " << nI.FirstBoundaryP()->getID() << "\n";
   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: "
        << nedges << " NT: " << ntri << endl;

   // remove closed boundary nodes that do not have non-boundary nbrs:
   for( cn = nI.FirstBoundaryP(); !( nI.AtEnd() ); cn = nI.NextP() )
   {
      numActNbrs = 0;
      sI.Reset( cn->getSpokeListNC() );
      for( ce = sI.FirstP(); !( sI.AtEnd() ); ce = sI.NextP() )
          if( ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary )
              ++numActNbrs;
      if( numActNbrs ) cn->setZ( minz );
      else deletelist.insertAtBack( cn ); 
   }
   for( cn = dI.FirstP(); !( dI.AtEnd() ); cn = dI.FirstP() )
   {
      DeleteNode( cn, kNoRepair );
      deletelist.removeFromFront( cn );
   }
   cout << "deleted superfluous boundary nodes\n";
   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: "
        << nedges << " NT: " << ntri << endl;

   // Update Voronoi areas, edge lengths, etc., and test the consistency
   // of the new mesh.
   UpdateMesh();
}


/*****************************************************************************\
**
**  tMesh::CheckMeshConsistency
**
**  Performs a series of tests to make sure the mesh connectivity is correct.
**  Should be called immediately after reading in a user-defined mesh (but
**  of course can also be used for debugging).
**
**  The consistency checks include the following:
**
**  1) Each edge:
**     - Has valid origin and destination pointers
**     - Has a valid counter-clockwise edge, which shares the same origin but
**       not the same destination
**     - Is paired with its complement in the list
**
**  2) Each node:
**     - Points to a valid edge which has the node as its origin
**     - If the node is not a boundary, it has at least one neighbor that
**       is not a closed boundary (unless boundaryCheckFlag is FALSE).
**     - Has a consistent spoke list (ie, you can go around the spokes and
**       get back to where you started)
**
**  3) Each triangle:
**     - Has 3 valid points and edges
**     - Each edge Ei has Pi as its origin and P((i+2)%3) as its
**       destination
**     - If an opposite triangle Ti exists, points P((i+1)%3) and
**       P((i+2)%3) are the same as points PO((n+2)%3) and PO((n+1)%3) in
**       the opposite triangle, where PO denotes a point in the opposite
**       triangle and n is the vertex ID (0, 1, or 2) of the point in the
**       opposite triangle that is opposite from the shared face.
**     - If an opposite triange Ti does not exist, points P((i+1)%3) and
**       and P((i+2)%3) should both be boundary points.
**
**      Parameters:  boundaryCheckFlag -- defaults to TRUE; if FALSE,
**                                        node connection to open node or
**                                        open boundary isn't tested
**      Data members updated: none
**      Called by: 
**      Calls:
**      Notes: does not check whether ID's are within the valid range;
**             that's assumed to be taken care of by the input routines
**        
**      Created: GT 1/98
**      Modifications:
**        - 4/98 GT added boundaryCheckFlag with default TRUE, so that
**               boundary checks can be disabled when the routine is called
**               in the middle of mesh creation operation as a debug/test
**               helper.
**
\*****************************************************************************/
#define kMaxSpokes 100
template<class tSubNode>
void tMesh< tSubNode >::
CheckMeshConsistency( int boundaryCheckFlag ) /* default: TRUE */
{
   tMeshListIter<tSubNode> nodIter( nodeList );
   tMeshListIter<tEdge> edgIter( edgeList );
   tListIter<tTriangle> triIter( triList );
   tPtrListIter< tEdge > sIter;
   tNode * cn, * org, * dest;
   tEdge * ce, * cne, * ccwedg;
   tTriangle * ct, * optr;
   int boundary_check_ok, i, nvop;

   // Edges: make sure complementary pairs are together in the list
   // (each pair Ei and Ei+1, for i=0,2,4,...nedges-1, should have the same
   // endpoints but the opposite orientation)
   for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() )
   {
      cne = edgIter.NextP();
      if( ce->getOriginPtrNC() != cne->getDestinationPtrNC()
          || ce->getDestinationPtrNC() != cne->getOriginPtrNC() )
      {
          cerr << "EDGE #" << ce->getID()
               << " must be followed by its complement in the list\n";
          goto error;
      }
      
   }

   // Edges: check for valid origin, destination, and ccwedg
   for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() )
   {
      if( !(org=ce->getOriginPtrNC() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not have a valid origin point\n";
         goto error;
      }
      if( !(dest=ce->getDestinationPtrNC() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not have a valid destination point\n";
         goto error;
      }
      if( !(ccwedg=ce->getCCWEdg() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not point to a valid counter-clockwise edge\n";
         goto error;
      }
      if( ccwedg->getOriginPtrNC()!=org )
      {
         cerr << "EDGE #" << ce->getID()
              << " points to a CCW edge with a different origin\n";
         goto error;
      }
      if( ccwedg->getDestinationPtrNC()==dest )
      {
         cerr << "EDGE #" << ce->getID()
              << " points to a CCW edge with the same destination\n";
         goto error;
      }
      if( org==dest )
      {
         cerr << "EDGE #" << ce->getID()
              << " has the same origin and destination nodes\n";
         goto error;
      }
      
   }
   //cout << "EDGES PASSED\n";

   // Nodes: check for valid edg pointer, spoke connectivity, and connection
   // to at least one non-boundary or open boundary node
   for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP() )
   {
      // edg pointer
      if( !(ce = cn->getEdg()) )
      {
         cerr << "NODE #" << cn->getID()
              << " does not point to a valid edge\n";
         goto error;
      }
      if( ce->getOriginPtrNC()!=cn )
      {
         cerr << "NODE #" << cn->getID()
              << " points to an edge that has a different origin\n";
         goto error;
      }

      // Boundary check and spoke consistency: if node is NOT a boundary,
      // it should be adjacent to at least one non-boundary or open boundary
      // point. Here we also test for an infinite loop in spoke connectivity.
      //   (Note that the boundary test always passes if the boundaryCheckFlag
      // is FALSE, meaning that we're in the middle of an operation that
      // could legitimately add open points without connection to an
      // open node or boundary --- this is added to allow for frequent
      // consistency checks even in the middle of mesh creation operations,
      // for testing/debugging purposes).
      boundary_check_ok = ( cn->getBoundaryFlag()==kNonBoundary &&
                            boundaryCheckFlag ) ? 0 : 1;
      i = 0;
      // Loop around the spokes until we're back at the beginning
      do
      {
         
         if( ce->getDestinationPtrNC()->getBoundaryFlag()!=kClosedBoundary )
             boundary_check_ok = 1;  // OK, there's at least one open nbr
         i++;
         if( i>kMaxSpokes ) // Uh-oh, an infinite loop
         {
            cerr << "NODE #" << cn->getID()
                 << ": infinite loop in spoke connectivity\n";
            goto error;
         }

         // Make sure node is the origin --- and not the destination
         if( ce->getOriginPtrNC()!=cn )
         {
            cerr << "EDGE #" << ce->getID()
                 << " is in the spoke chain of NODE " << cn->getID()
                 << " but does not have the node as an origin\n";
            goto error;
         }
         if( ce->getDestinationPtrNC()==cn )
         {
            cerr << "EDGE #" << ce->getID()
                 << " is in the spoke chain of NODE " << cn->getID()
                 << " but has the node as its destination\n";
            goto error;
         }   
         
      } while( (ce=ce->getCCWEdg())!=cn->getEdg() );
      if( !boundary_check_ok )
      {
         cerr << "NODE #" << cn->getID()
              << " is surrounded by closed boundary nodes\n";
         goto error;
      }
      
      //make sure node coords are consistent with edge endpoint coords:
      sIter.Reset( cn->getSpokeListNC() );
      for( ce = sIter.FirstP(); !(sIter.AtEnd()); ce = sIter.NextP() )
      {
         if( ce->getOriginPtrNC()->getX() != cn->getX() ||
             ce->getOriginPtrNC()->getY() != cn->getY() )
         {
            cerr << "NODE #" << cn->getID()
                 << " coords don't match spoke origin coords\n";
            goto error;
         }
      }
      
   }
     //cout << "NODES PASSED\n";

   // Triangles: check for valid points and connectivity
   for( ct=triIter.FirstP(); !(triIter.AtEnd()); ct=triIter.NextP() )
   {
      for( i=0; i<=2; i++ )
      {
         // Valid point i?
         if( !(cn=ct->pPtr(i)) )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << " has an invalid point " << i << endl;
            goto error;
         }
         // Valid edge i?
         if( !(ce=ct->ePtr(i)) )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << " has an invalid edge " << i << endl;
            goto error;
         }
         // Edge and point consistency
         if( ce->getOriginPtrNC()!=cn )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << ": edge " << i << " does not have point " << i
                 << " as origin\n";
            goto error;
         }
         // changed from (i+1) to (i+2) for "right-hand" format gt 3/98
         if( ce->getDestinationPtrNC()!=ct->pPtr((i+2)%3) )
         {
            cerr << "TRIANGLE #" << ct->getID()
                 << ": edge " << i << " does not have point " << (i+1)%3
                 << " as destination\n";
            goto error;
         }
         // Opposite triangle: if it exists, check common points
         if( (optr = ct->tPtr(i)) )
         {
            nvop = optr->nVOp(ct); // Num (0,1,2) of opposite vertex in optr
            if( nvop < 3 )
            {
               if( ct->pPtr((i+1)%3) != optr->pPtr((nvop+2)%3)
                   || ct->pPtr((i+2)%3) != optr->pPtr((nvop+1)%3) )
               {
                  cerr << "TRIANGLE #" << ct->getID()
                       << ": opposite triangle " << i << " does not share nodes "
                       << (ct->pPtr((i+1)%3))->getID() << " and "
                       << (ct->pPtr((i+2)%3))->getID() << endl;
                  goto error;
               }
            }
            else
            {
               cerr << "TRIANGLE #" << ct->getID()
                    << ": opposite triangle " << i << ", triangle #"
                    << optr->getID() << ", does not have current tri as neighbor\n";
               goto error;
            }
         }
         // If no opposite triangle, make sure it really is a boundary
         else
         {
            if( (ct->pPtr((i+1)%3))->getBoundaryFlag()==kNonBoundary
                || (ct->pPtr((i+2)%3))->getBoundaryFlag()==kNonBoundary )
            {
               cerr << "TRIANGLE #" << ct->getID()
                    << ": there is no neighboring triangle opposite node "
                    << cn->getID() << " but one (or both) of the other nodes "
                    << "is a non-boundary point\n";
               goto error;
            }
         }       
      }
   }
     //cout << "TRIANGLES PASSED\n";
   //cout << "MESH PASSED\n";
   return;
   
  error:
   ReportFatalError( "Error in mesh consistency." );
   
}
#undef kMaxSpokes

template< class tSubNode >
void tMesh< tSubNode >::
Print()                                                  //tMesh
{
   triList.print();
   nodeList.print();
   edgeList.print();
}

template< class tSubNode >
void tMesh< tSubNode >::
MakeCCWEdges()
{
   tMeshListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn;
   for( cn = nodIter.FirstP(); !( nodIter.AtEnd() ); cn = nodIter.NextP() )
   {
      cn->makeCCWEdges();
   }
}


/*****************************************************************************\
**
**  tMesh::setVoronoiVertices
**
**  Each Delaunay triangle is associated with an intersection between
**  three Voronoi cells, called a Voronoi vertex. These Voronoi vertices
**  are used in computing the area of each Voronoi cell. The Voronoi
**  vertex associated with each triangle is the circumcenter of the
**  triangle. This routine finds the Voronoi vertex associated with
**  each triangle by finding the triangle's circumcenter. 
**
**  The vertex coordinates are stored in the three clockwise-
**  oriented tEdge objects associated with each triangle. This is a 
**  space-for-time tradeoff: the coordinates could be stored in the triangles, 
**  saving redundancy (3 copies of each point are stored here), but in
**  that case each tEdge would have to point back to a triangle, and an
**  additional level of indirection would be needed in accessing the Voronoi
**  vertices associated with a particular node.
**
**  Note that computation of the Voronoi vertices can be prone to numerical
**  errors, leading to inconsistent Voronoi polygons with cells that
**  overlap or have loops. See note under tNode::ComputeVoronoiArea().
**
**    Assumes: correct triangulation with valid edge pointers in each tri.
**    Data members modified: none
**    Other objects modified: Voronoi vertices set for each tEdge
**    Modifications:
**     - reverted to earlier triangle-based computation, from an edge-based
**       computation that takes 3x as long because NE = 3NT. In so doing,
**       the definition of the Voronoi vertex stored in a tEdge is changed
**       to "left-hand", meaning the V. vertex associated with the edge's
**       lefthand triangle (the vertex itself may or may not lie to the left
**       of the edge). 1/98 GT
**     - also moved circumcenter computation into a tTriangle mbr fn.
**     - copied function to tMesh member from tStreamNet member, gt 3/98.
**       Other fns now use "right-hand" definition; this fn may have to
**       be changed.
**
\*****************************************************************************/
template <class tSubNode>
void tMesh<tSubNode>::setVoronoiVertices()
{
   //double x, y, x1, y1, x2, y2, dx1, dy1, dx2, dy2, m1, m2;
   //tArray< double > xyo, xyd1, xyd2, xy(2);
   //cout << "setVoronoiVertices()..." << endl;
   tArray< double > xy;
   tListIter< tTriangle > triIter( triList );
   tTriangle * ct;

   // Find the Voronoi vertex associated with each Delaunay triangle
   for( ct = triIter.FirstP(); !(triIter.AtEnd()); ct = triIter.NextP() )
   {
      xy = ct->FindCircumcenter();    
      //cout << "setVoronoiVertices(): " << xy[0] << " " << xy[1];
      // Assign the Voronoi point as the left-hand point of the three edges 
      // associated with the current triangle
      ct->ePtr(0)->setRVtx( xy );
      ct->ePtr(1)->setRVtx( xy );
      ct->ePtr(2)->setRVtx( xy );

      // debug output
      /*cout << "FOR edges: ";
      int i;
      for( i=0; i<=2; i++ )
          cout << ct->ePtr(i)->getID() << " ("
               << ct->ePtr(i)->getOriginPtr()->getID() << ","
               << ct->ePtr(i)->getDestinationPtr()->getID() << ") ";
      cout << ", v verts are:\n";
      xy = ct->ePtr(0)->getRVtx();
      cout << "  setVoronoiVertices(): " << xy[0] << " " << xy[1] << endl;*/
   }
   //cout << "setVoronoiVertices() finished" << endl;
}


/**************************************************************************\
**
**  tMesh::CalcVoronoiEdgeLengths
**
**  Updates the length of the Voronoi cell edge associated with each
**  triangle edge. Because complementary edges are stored pairwise on
**  the edge list, we can save computation time by only computing the
**  vedglen once for the first of the pair, then assigning it to the
**  second. For boundary triangle edges, the corresponding Voronoi edge
**  is infinitely long, so the calculation is only done for interior
**  (active) edges.
**
**  Data mbrs modified:  (none)
**  Calls:  tEdge::CalcVEdgLen, tEdge::setVEdgLen
**  Assumes:  complementary edges are stored pairwise on the edge list;
**				Voronoi vertices are up to date
**  Notes:  replaces tNode::CalcSpokeVEdgLengths()
**
\**************************************************************************/
template <class tSubNode>
void tMesh<tSubNode>::CalcVoronoiEdgeLengths()
{
	tEdge *ce;
	double vedglen;
	tMeshListIter<tEdge> edgIter( edgeList );

	for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
	{
		vedglen = ce->CalcVEdgLen();  // Compute Voronoi edge length
		ce = edgIter.NextP();         // Advance to complement edge and
		ce->setVEdgLen( vedglen );    //   and assign the same edge length.
	}
}


/**************************************************************************\
**
**  tMesh::CalcVAreas
**
**  Computes Voronoi area for each active (non-boundary) node in the
**  mesh (Voronoi area is only defined for interior nodes). Accomplishes
**  this by calling ComputeVoronoiArea for each node. (see meshElements)
**
\**************************************************************************/
template <class tSubNode>
 void tMesh<tSubNode>::CalcVAreas()
{
   //cout << "CalcVAreas()..." << endl << flush;
   //Xdouble area;
   tSubNode* curnode;
   //XtEdge *ce;
   //XtArray< double > xy;
   tMeshListIter< tSubNode > nodIter( nodeList );
   
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
         //area = VoronoiArea( curnode );
      curnode->ComputeVoronoiArea();
      
   }
   //cout << "CalcVAreas() finished" << endl;
}


/**************************************************************************\
**
**  tMesh::DeleteNode( tListNode<tSubNode> *, int =1 )
**    (see DeleteNode( tSubNode *, int =1 ) below)
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
DeleteNode( tListNode< tSubNode > *nodPtr, int repairFlag )
{
   //cout << "DeleteNode: " << nodPtr->getDataPtr()->getID() << endl;
   if( !DeleteNode( nodPtr->getDataPtrNC(), repairFlag ) ) return 0;
   return 1;   
}


/**************************************************************************\
**
**  tMesh::DeleteNode( tSubNode *, int =1 )
**
**  Deletes a node from the mesh. This is done by first calling
**  ExtricateNode to detach the node by removing its edges and their
**  associated triangles, then removing the node from the nodeList.
**  Normally, RepairMesh is then called to retriangulate the "hole" left
**  behind in the mesh. (However, if the node was on the hull of the
**  mesh there's no "hole" to fix --- the caller is assumed to be smart
**  enough to recognize this case and let us know about it by setting
**  repairFlag to kNoRepair. This is the case, for example, when deleting
**  the nodes that form a "supertriangle" as in MakeMeshFromPoints).
**
**  Once the mesh is repaired, the nodes are renumbered and as a safety
**  measure for debugging/testing purposes, UpdateMesh is called.
**
**  Data mbrs modified:  nnodes, nedges, and ntri are updated;
**                       the node is deleted from nodeList; other edges &
**                       triangles are removed and/or modified by
**                       ExtricateNode and RepairMesh (qv)
**  Calls:  tMesh::ExtricateNode, tMesh::RepairMesh, plus utility member
**               functions of tNode, tMeshList, etc. (and temporarily,
**               UpdateMesh)
**  Returns:  error code: 0 if either ExtricateNode or RepairMesh fails,
**            1 otherwise.
**  Assumes:  
**  Notes:
**    - repairFlag defaults to TRUE if not specified
**    - if node is on the hull and repairFlag is TRUE, the result can be
**        an infinite loop in RepairMesh as it tries to mend a hole that
**        doesn't exist. This condition isn't tested for, so be careful.
**  Created: fall, '97 SL
**  Modifications: added repairFlag 4/98 GT
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
DeleteNode( tSubNode *node, int repairFlag )
{
   tPtrList< tSubNode > nbrList;
   tListNode< tSubNode > *nodPtr;
   tMeshListIter< tSubNode > nodIter( nodeList );
   nodIter.Get( node->getID() );
   tSubNode nodeVal;
   
   //cout << "DeleteNode: " << node->getID() << " at " << node->getX() << " "
   //<< node->getY() << " " << node->getZ() << endl;
   //assert( repairFlag || node->getBoundaryFlag()==kClosedBoundary );
   
   nodPtr = nodIter.NodePtr();
   if( !( ExtricateNode( node, nbrList ) ) ) return 0;
   if( node->getBoundaryFlag() )
   {
      nodeList.moveToBack( nodPtr );
      nodeList.removeFromBack( nodeVal );
   }
   else
   {
      nodeList.moveToFront( nodPtr );
      nodeList.removeFromFront( nodeVal );
   }
   
   //cout << "Removed node " << nodeVal.getID() << " at x, y "
   // << nodeVal.getX() << ", " << nodeVal.getY() << "; " << endl;
   nnodes = nodeList.getSize();
   nedges = edgeList.getSize();
   ntri = triList.getSize();
   //cout << "nn " << nnodes << "  ne " << nedges << "  nt " << ntri << endl;
   
   /*tPtrListIter< tSubNode > nbrIter( nbrList );
   cout << "leaving hole defined by " << endl << "   Node  x  y " << endl;
   for( i=0, nbrIter.First(); nbrIter.NextIsNotFirst(); i++ )
   {
      if( i>0 ) nbrIter.Next();
      cout << "   " << nbrIter.DatPtr()->getID() << "     "
           << nbrIter.DatPtr()->getX() << "  "
           << nbrIter.DatPtr()->getY() << endl;
   }*/

   if( repairFlag )
   {
      nbrList.makeCircular();
      if( !RepairMesh( nbrList ) ) return 0;
   }
   
   //reset node id's
   assert( nodIter.First() );
   miNextNodeID = 0;
   do
   {
      nodIter.DatRef().setID( miNextNodeID );
      miNextNodeID++;
   }
   while( nodIter.Next() );
     //cout << "Mesh repaired" << endl;
   /*tSubNode *cn;
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      cout << "node " << cn->getID() << endl;
   }*/

   if( repairFlag ) UpdateMesh();
   return 1;
}


/**************************************************************************\
**
**  tMesh::ExtricateNode
**
**  Detaches a node from the mesh by deleting all of its edges (which in
**  turn removes the affected triangles). Returns a list of the node's
**  former neighbors by modifying the nbrList input parameter. Also
**  returns a code that indicates failure if the node still has a non-empty
**  spoke list after edge deletion.
**
**  Data mbrs modified:  nnodes; edges and triangles are removed from
**                       edgeList and triList
**  Calls:  tMesh::DeleteEdge and utility member functions of tNode,
**               tPtrList, tPtrListIter
**  Output:  list of node's (former) neighbors, in nbrList
**  Returns:  1 if all edges successfully deleted, 0 if not
**  Calls: DeleteEdge
**  Assumes:  
**  Notes:
**  Created: SL fall, '97
**  Modifications: if node is a closed boundary, any of its neighbors that
**             are non-boundaries are switched to closed boundaries, so
**             that nodes along the edge of the domain (including nodes of
**             a "supertriangle" used in MakeMeshFromPoints) may be removed
**             without causing errors, GT 4/98
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
ExtricateNode( tSubNode *node, tPtrList< tSubNode > &nbrList )
{
   //cout << "ExtricateNode: " << node->getID() << endl;
   tPtrListIter< tEdge > spokIter( node->getSpokeListNC() );
   tEdge edgeVal1, edgeVal2, *ce;
   tSubNode *nbrPtr;
   
     //cout << "Removing spokes: " << endl;
     //assert( ExtricateEdge( edgptrIter.DatPtr() ) );
   for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.FirstP() )
   {
      nbrPtr = static_cast< tSubNode * >(ce->getDestinationPtrNC());
      nbrList.insertAtBack( nbrPtr );
      if( node->getBoundaryFlag()                      // If node is a bdy make
          && nbrPtr->getBoundaryFlag()==kNonBoundary )// sure nbrs are also
      {                                                // boundaries.
         nbrPtr->ConvertToClosedBoundary();
         nodeList.moveToBack( nbrPtr );
      }
      if( !DeleteEdge( ce ) ) return 0;
   }  

   //cout<<"nnodes decremented now "<<nnodes<<endl;
   if( node->getSpokeList().isEmpty() ) return 1;
   return 0;
}

/**************************************************************************\
**
**  tMesh::DeleteEdge
**
**  Deletes an edge from the mesh, returning 1 if deletion succeeds and
**  0 if not. Starts by calling ExtricateEdge to detach the edge from
**  the other mesh elements. This function actually deletes two directed
**  edges: edgePtr and its complement.
**
**  Inputs:  edgePtr -- ptr to the edge to be deleted
**  Returns:  1 if successful, 0 if not
**  Calls: ExtricateEdge 
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
DeleteEdge( tEdge * edgePtr )
{
   //cout << "DeleteEdge(...) " << edgePtr->getID() << endl;
   //edgePtr->TellCoords();
   tEdge edgeVal1, edgeVal2;

   // Detach the edge from other mesh elements
   if( !ExtricateEdge( edgePtr ) ) return 0;

   // Remove the edge and its complement from the list
   // Note, extricate edge does not actually remove the edge from
   // the edge list, it only moves the two edges (one 'line' that
   // has two directions) to the end or front of the list.
   // These two edges are then actually removed here.
   if( edgePtr->getBoundaryFlag() )
   {
      if( !( edgeList.removeFromBack( edgeVal1 ) ) ) return 0;
      if( !( edgeList.removeFromBack( edgeVal2 ) ) ) return 0;
   }
   else
   {
      if( !( edgeList.removeFromFront( edgeVal1 ) ) ) return 0;
      if( !( edgeList.removeFromFront( edgeVal2 ) ) ) return 0;
   }
   //    cout << "  edges " << edgeVal1.getID() << " and "
   //         <<  edgeVal2.getID() << " between nodes "
   //         << edgeVal1.getOriginPtr()->getID() << " and "
   //         << edgeVal2.getOriginPtr()->getID() << " removed" << endl;
   
   if( &edgeVal1 == 0 || &edgeVal2 == 0 ) return 0;
   return 1;
}


/**************************************************************************\
**
**  tMesh::ExtricateEdge
**
**  Here we detach an edge and its complement from the surrounding mesh 
**  elements prior to deletion. Adjacent triangle(s) are also deleted
**  via a call to DeleteTriangle. Calls the virtual node function
**  WarnSpokeLeaving to signal the affected nodes to take appropriate
**  action. (Appropriate action might depend on the application; that's
**  why it is a virtual function that can be handled by any descendents
**  of tNode). The two complementary edges are then placed at the back
**  of the edge list, where DeleteEdge can conveniently find them.
**
**  Inputs:  edgePtr -- ptr to the edge to be deleted
**  Returns: 1 if successful, 0 otherwise
**  Calls: DeleteTriangle, <tSubNode>::WarnSpokeLeaving
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
ExtricateEdge( tEdge * edgePtr )
{
   //cout << "ExtricateEdge: " << edgePtr->getID() << endl;
   //edgePtr->TellCoords();
   assert( edgePtr != 0 );
   //temporary objects:
   tEdge *tempedgePtr=0, *ce, *cce, *spk;
   tEdge *ceccw, *cceccw;
   tMeshListIter< tEdge > edgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   tPtrList< tEdge > *spkLPtr;
   tListNode< tEdge > *listnodePtr;
   tArray< tTriangle * > triPtrArr(2);

   //cout << "find edge in list; " << flush;
   ce = edgIter.GetP( edgePtr->getID() );  //NB: why necessary? isn't ce the
   // same as edgePtr??  Yes, why???
   // Puts edgIter at edgePtr's node!
   // Remove the edge from it's origin's spokelist
   //cout << "update origin's spokelist if not done already; " << flush;
   spkLPtr = &( ce->getOriginPtrNC()->getSpokeListNC() );
   spokIter.Reset( *spkLPtr );
   for( spk = spokIter.FirstP(); spk != ce && !( spokIter.AtEnd() ); spk = spokIter.NextP() );
   if( spk == ce )
   {
      spk = spokIter.NextP();
      spkLPtr->removePrev( tempedgePtr, spokIter.NodePtr() );
   }

   // Find the triangle that points to the edge
   //cout << "find triangle; " << flush;
   triPtrArr[0] = TriWithEdgePtr( edgePtr ); 
   // Find the edge's complement
   listnodePtr = edgIter.NodePtr();
   assert( listnodePtr != 0 );
   
   //cout << "find complement; " << flush;
   if( edgePtr->getID()%2 == 0 ) cce = edgIter.NextP();
   else if( edgePtr->getID()%2 == 1 ) cce = edgIter.PrevP();
   else return 0; //NB: why whould this ever occur??

   // Find the triangle that points to the edges complement
   //cout << "find other triangle; " << flush;
   triPtrArr[1] = TriWithEdgePtr( cce );
   //if triangles exist, delete them
   //cout << "conditionally calling deletetri from extricateedge\n";
   if( triPtrArr[0] != 0 )
       if( !DeleteTriangle( triPtrArr[0] ) ) return 0;
   if( triPtrArr[1] != 0 )
       if( !DeleteTriangle( triPtrArr[1] ) ) return 0;
   //update complement's origin's spokelist
   spkLPtr = &(cce->getOriginPtrNC()->getSpokeListNC());
   spokIter.Reset( *spkLPtr );
   for( spk = spokIter.FirstP(); spk != cce && !( spokIter.AtEnd() );
        spk = spokIter.NextP() );
   if( spk == cce )
   {
      spk = spokIter.NextP();
      spkLPtr->removePrev( tempedgePtr, spokIter.NodePtr() );
   }

   //Need to make sure that edg member of node was not pointing
   //to one of the edges that will be removed.  Also, may be implications
   //for some types of subnodes, so take care of that also.
   tSubNode * nodece = static_cast<tSubNode *>(ce->getOriginPtrNC());
   nodece->WarnSpokeLeaving( ce ); 
   tSubNode * nodecce = static_cast<tSubNode *>(cce->getOriginPtrNC());
   nodecce->WarnSpokeLeaving( cce );

   //now, take care of the edges who had as thier ccwedge ce or cce
   ceccw=ce->getCCWEdg();
   tempedgePtr=ceccw;
   do{
      tempedgePtr=tempedgePtr->getCCWEdg();
   }while(tempedgePtr->getCCWEdg() != ce);
   //set tempedgeptrs ccwedge to ceccw
   tempedgePtr->setCCWEdg( ceccw);
   
   cceccw=cce->getCCWEdg();
   tempedgePtr=cceccw;
   do{
      tempedgePtr=tempedgePtr->getCCWEdg();
   }while(tempedgePtr->getCCWEdg() != cce);
   //set tempedgeptrs ccwedge to cceccw
   tempedgePtr->setCCWEdg(cceccw);

   //Since WarnSpokeLeaving can set a node to a boundary node if
   //There is no longer a legit place to flow, we need to check 
   //to see if nodece and nodecce are now boundaries, and take proper action.

   if( ce->getBoundaryFlag() )
   {
      //move edges to back of list
      edgeList.moveToBack( listnodePtr );
      edgeList.moveToBack( edgIter.NodePtr() );
   }
   else
   {
      //move edges to front of list
      edgeList.moveToFront( edgIter.NodePtr() );
      edgeList.moveToFront( listnodePtr );
   }
   nedges-=2;
   return 1;
}


/***************************************************************************\
**
**  tMesh::LocateTriangle
**
**  Locates the triangle in which point (x,y) falls. The algorithm exploits
**  the fact that the 3 triangle points are always in counter-clockwise
**  order, so that the point is contained within a given triangle (p0,p1,p2)
**  if and only if the point lies to the left of vectors p0->p1, p1->p2,
**  and p2->p0. Here's how it works:
**   1 - start with a given triangle (currently first on the list, but a
**       smarter initial guess could be used -- TODO)
**   2 - lv is the number of successful left-hand checks found so far:
**       initialize it to zero
**   3 - check whether (x,y) lies to the left of p(lv)->p((lv+1)%3)
**   4 - if so, increment lv by one (ie, move on to the next vector)
**   5 - if not, (x,y) is to the right of the current face, so move to
**       the triangle that lies opposite that face and reset lv to zero
**   6 - continue steps 3-5 until lv==3, which means that we've found
**       our triangle.
**   7 - so far, a point "on the line", i.e., colinear w/ two of the
**       three points, still passes; that's OK unless that line is on
**       the boundary, so we need to check
**
**  Input: x, y -- coordinates of the point
**  Modifies: (nothing)
**  Returns: a pointer to the triangle that contains (x,y)
**  Assumes: the point is contained within one of the current triangles
**
\***************************************************************************/
template< class tSubNode >
tTriangle * tMesh< tSubNode >::
LocateTriangle( double x, double y )
{
   //cout << "\nLocateTriangle (" << x << "," << y << ")\n";
   int n, lv=0;
   tListIter< tTriangle > triIter( triList );  //lt
   //XtTriangle *lt = &(triIter.DatRef());
   tTriangle *lt = ( mSearchOriginTriPtr > 0 ) ? mSearchOriginTriPtr
       : triIter.FirstP();
   double a, b, c;
   int online = -1;
   tArray< double > xy1, xy2;
   
   /* it starts from the first triangle, 
      searches through the triangles until the point is on
      the same side of all the edges of a triangle.
      "lt" is the current triangle and "lv" is the edge number. */
   for (n=0 ;(lv!=3)&&(lt); n++)
   {
      xy1 = lt->pPtr(lv)->get2DCoords();
      xy2 = lt->pPtr( (lv+1)%3 )->get2DCoords();
      a = (xy1[1] - y) * (xy2[0] - x);
      b = (xy1[0] - x) * (xy2[1] - y);
      c = a - b;
      //cout << "find tri for point w/ x, y, " << x << ", " << y
      //  << "; no. tri's " << ntri << "; now at tri " << lt->getID() << endl;
      // lt->TellAll();
        cout << flush;
      
      if ( c > 0.0 )
      {
         //cout << "    Moving on...\n";
         lt=lt->tPtr( (lv+2)%3 );
         lv=0;
         online = -1;
      }
      else
      {
         //cout << "    So far so good...\n";
         if( c == 0.0 ) online = lv;
         lv++;
      }
      
      /*if( n >= ntri + 20 )
        {
        DumpTriangles();
        DumpNodes();
        }*/
      //cout << "NTRI: " << ntri << flush;
      assert( n < 3*ntri );
   }
   //cout << "FOUND point in:\n";
   //if( lt != 0 ) lt->TellAll(); //careful with this! TellAll() will crash
   //else ReportFatalError( "point out of bounds" );
   //if lt == 0, i.e., point is out of bounds,
   //and we don't want that;
   //calling code is built to deal with lt == 0.
   if( online != -1 )
       if( lt->pPtr(online)->getBoundaryFlag() != kNonBoundary &&
           lt->pPtr( (online+1)%3 )->getBoundaryFlag() != kNonBoundary ) //point on bndy
           return 0;
   //else cout << "location out of bounds\n";
   return(lt);
}


/**************************************************************************\
**
**  tMesh::LocateNewTriangle
**
**  HOW IS THIS DIFFERENT FROM LOCATETRI?
**
**  Called by: AddNodeAt
**
\**************************************************************************/
template< class tSubNode >
tTriangle * tMesh< tSubNode >::
LocateNewTriangle( double x, double y )
{
   //cout << "LocateTriangle" << endl;
   int n, lv=0;
   tListIter< tTriangle > triIter( triList );  //lt
   tTriangle *lt = triIter.FirstP();
   tSubNode *p1, *p2;
   
   tArray< double > xy1, xy2;
   /* it starts from the first triangle, 
      searches through the triangles until the point is on
      the same side of all the edges of a triangle.
      "lt" is the current triangle and "lv" is the edge number.
   */
   for (n=0 ;(lv!=3)&&(lt); n++)
   {
      p1 = static_cast<tSubNode *>(lt->pPtr(lv));
      if( p1->Meanders() ) xy1 = p1->getNew2DCoords();
      else xy1 = p1->get2DCoords();
      p2 = static_cast<tSubNode *>(lt->pPtr( (lv+1)%3 ));
      if( p2->Meanders() ) xy1 = p2->getNew2DCoords();
      else xy2 = p2->get2DCoords();
      if ( ( (xy1[1] - y) * (xy2[0] - x) ) > ( (xy1[0] - x) * (xy2[1] - y)) )
      {
         lt=lt->tPtr( (lv+2)%3 );
         lv=0;
      }
      else {lv++;}
      /*if( !(n < ntri) )
        cout << "tri not found for point w/ x, y, " << x << ", " << y
        << "; no. tri's " << ntri << "; now at tri " << lt->getID() << endl;*/
      //assert( n < ntri + 20 );
   }
   return(lt);
}


/**************************************************************************\
**
**  tMesh::TriWithEdgePtr
**
**  Finds and returns the triangle that points to edgPtr as one of its
**  clockwise-oriented edges.
**
\**************************************************************************/
template< class tSubNode >
tTriangle *tMesh< tSubNode >::
TriWithEdgePtr( tEdge *edgPtr )
{
   assert( edgPtr != 0 );
   tTriangle *ct;
   //cout << "TriWithEdgePtr " << edgPtr->getID();
   tListIter< tTriangle > triIter( triList );

   for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
       if( ct != 0 ) //TODO: is this test nec? why wd it be zero?
           if( ct->ePtr(0) == edgPtr ||
               ct->ePtr(1) == edgPtr ||
               ct->ePtr(2) == edgPtr ) return ct;
   return 0;
}


/**************************************************************************\
**
**  tMesh::DeleteTriangle
**
**  Deletes a triangle from the mesh. Starts off with a call to 
**  ExtricateTriangle to detach the triangle from other mesh elements,
**  after which the triangle is at the front of the triangle list,
**  from whence it is then deleted.
**
**  Inputs:  triPtr -- ptr to the triangle to be deleted
**  Returns:  1 if successful, 0 if not
**  Calls: ExtricateTriangle
**  Called by: DeleteEdge, AddNode, AddNodeAt
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
DeleteTriangle( tTriangle * triPtr )
{
   //cout << "DeleteTriangle(...) " << triPtr->getID() << endl;
   //triPtr->TellAll();
   tTriangle triVal;

   if( !ExtricateTriangle( triPtr ) ) return 0;
   //if( !( triList.removeFromFront( triVal ) ) ) return 0;
   if( !( triList.removeFromFront(triVal) ) )
   {
      cerr << "DeleteTriangle(): triList.removeFromFront( triPtr ) failed\n";
      return 0;
   }
   if( &triVal == 0 ) return 0;
   return 1;
}

/**************************************************************************\
**
**  tMesh::ExtricateTriangle
**
**  Detaches a triangle from surrounding mesh elements and places it at
**  the head of the triangle list, where it can be easily deleted by
**  DeleteTriangle.
**
**  Inputs: triPtr -- ptr to the triangle to be extricated
**  Returns: 1 if successful, 0 if not
**  Called by: DeleteTriangle
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
ExtricateTriangle( tTriangle *triPtr )
{
   //cout << "ExtricateTriangle" << endl;
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct;

   // Find the triangle on the list
   for( ct = triIter.FirstP(); ct != triPtr && !( triIter.AtEnd() );
        ct = triIter.NextP() );
   if( ( triIter.AtEnd() ) ) return 0;

   // Tell your neighboring triangles that you're about to disappear by
   // setting their corresponding tPtr's to zero
   int i, j;
   for( i=0; i<3; i++ ) for( j=0; j<3; j++ )
       if( triPtr->tPtr(i) != 0 )
           if( triPtr->tPtr(i)->tPtr(j) == triPtr ) 
               triPtr->tPtr(i)->setTPtr( j, 0 );

   // Move the triangle to the head of the list where it can be deleted
   triList.moveToFront( triIter.NodePtr() );

   ntri--;

   // Make sure mSearchOriginTriPtr doesn't point to the triangle being
   // extricated. If it does, reassign it to a neighbor or zero..
   if( triPtr == mSearchOriginTriPtr )
   {
       mSearchOriginTriPtr = 0;
       for( i=0; i<3; ++i )
           if( triPtr->tPtr(i) != 0 )
               mSearchOriginTriPtr = triPtr->tPtr(i);
   }

   return 1;
}


/**************************************************************************\
**
**  tMesh::RepairMesh
**
**  This function repairs the "hole" in the mesh that is created when
**  a node is deleted. Essentially, this function stiches the hole back
**  together by adding edges and triangles as needed, preserving
**  Delaunay-ness. The nodes around the hole are stored in the input
**  parameter nbrList. As each new triangle is added, its "interior"
**  point is removed from the neighbor list. The process of stitching
**  proceeds iteratively until the hole itself is a Delaunay triangle.
**
**  For each set of 3 successive counter-clockwise points along the rim
**  of the whole, the function calls Next3Delaunay to compare the
**  potential triangle p0, p1, p2 with other potential triangles
**  p0, p1, ptest (where ptest is each of the other nodes along the rim).
**  When a Delaunay triangle is found, AddEdgeAndMakeTriangle is called
**  to create the necessary edge and triangle objects, and the interior
**  node is removed from the neighbor list. (or at least that's what
**  gt is able to deduce...)
**
**  Inputs: nbrList -- list of nodes surrounding the "hole"
**  Returns: 1 if successful, 0 if not
**  Calls: Next3Delaunay, AddEdgeAndMakeTriangle, MakeTriangle
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
RepairMesh( tPtrList< tSubNode > &nbrList )
{
   assert( &nbrList != 0 );
   //cout << "RepairMesh: " << endl;
   if( nbrList.getSize() < 3 ) return 0;
   tSubNode * meshnodePtr = 0;
   nbrList.makeCircular();
   tPtrListIter< tSubNode > nbrIter( nbrList );
   
   // Keep stitching until only 3 nodes are left
   while( nbrList.getSize() > 3 )
   {
      //cout << "in loop, nbr size = " << nbrList.getSize() << endl;
      //Xflowflag = 1;  // PURPOSE??
      if( Next3Delaunay( nbrList, nbrIter ) ) //checks for ccw and Del.
      {
           //cout << "found 3 Delaun!\n";
         AddEdgeAndMakeTriangle( nbrList, nbrIter );
           //remove "closed off" pt
         nbrList.removeNext( meshnodePtr, nbrIter.NodePtr() );
      }
        //else cout << "Not delaun\n";
      
      nbrIter.Next();                    //step forward once in nbrList
   }
   assert( nbrList.getSize() == 3 );
   assert( ntri == triList.getSize() );
   assert( nedges == edgeList.getSize() );
   assert( nnodes == nodeList.getSize() );       //make sure numbers are right
   MakeTriangle( nbrList, nbrIter );             //make final triangle
     //do some checking?
     //cout << "done" << endl;
   return 1;
}


/**************************************************************************\
**
**  tMesh::AddEdge
**
**  Function to add edge pair between two nodes. Resets edge IDs.
**
**  Inputs: three nodes; edge is added between first two, and third
**   should be CCW 3rd member of triangle.
**  Returns: 1 if successful, 0 if not
**
**  Created: SL fall, '97
**  Modified: SL 10/98--routine sometimes failed when node1 (or node2) had
**   had edges to neither node2 (or node1) nor node3; to fix, replaced
**   the "assert( !( spokIter.AtEnd() ) )"'s with new algorithm to find
**   where new spoke should be inserted: finds where the sequence of 3 spoke
**   unit vectors, including the new one in the middle, are CCW; calls new
**   global function, tArray< double > UnitVector( tEdge* ).
**   - GT 1/99 -- to avoid compiler warning, now stores output of 
**     UnitVector calls in arrays p1, p2, p3, which are then sent as
**     arguments to PointsCCW.
**   - GT 2/99 -- added calls to WelcomeCCWNeighbor and AttachNewSpoke
**     to update CCW edge connectivity
**
\**************************************************************************/
//vertices of tri in ccw order; edges are added between node1 and node2
//TODO: comment/document this fn
template< class tSubNode >
int tMesh< tSubNode >::
AddEdge( tSubNode *node1, tSubNode *node2, tSubNode *node3 ) 
{
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   /*cout << "AddEdge" << 
     "between nodes " << node1->getID()
     << " and " << node2->getID() << " w/ ref to node " << node3->getID() << endl;*/
   int flowflag = 1;  // Boundary code for new edges
   tEdge tempEdge1, tempEdge2;  // The new edges
   tEdge *ce, *le;
   tMeshListIter< tEdge > edgIter( edgeList );
   tMeshListIter< tSubNode > nodIter( nodeList );
   tPtrListIter< tEdge > spokIter;
   //XtPtrListNode< tEdge >* spokeNodePtr;
   tArray<double> p1, p2, p3;           // Used to store output of UnitVector

   // Set origin and destination nodes and find boundary status
   tempEdge1.setOriginPtr( node1 );   //set edge1 ORG
   tempEdge2.setDestinationPtr( node1 );//set edge2 DEST
   if( node1->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   tempEdge2.setOriginPtr( node2 );   //set edge2 ORG
   tempEdge1.setDestinationPtr( node2 );//set edge1 DEST
   if( node2->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   if( node1->getBoundaryFlag()==kOpenBoundary     // Also no-flow if both
       && node2->getBoundaryFlag()==kOpenBoundary ) //  nodes are open bnds
       flowflag = 0;
   
   // Set boundary status and ID
   /*Xif( !( edgeList.isEmpty() ) )
       newid = edgIter.LastP()->getID() + 1;
       else newid = 0;*/
   tempEdge1.setID( miNextEdgID );                     //set edge1 ID
   miNextEdgID++;
   tempEdge2.setID( miNextEdgID );                     //set edge2 ID
   miNextEdgID++;
   tempEdge1.setFlowAllowed( flowflag );         //set edge1 FLOWALLOWED
   tempEdge2.setFlowAllowed( flowflag );         //set edge2 FLOWALLOWED

   // Place new edge pair on the list: active back if not a boundary
   // edge, back otherwise
   if( flowflag == 1 )
   {
      edgeList.insertAtActiveBack( tempEdge1 );  //put edge1 active in list
      edgeList.insertAtActiveBack( tempEdge2 );  //put edge2 active in list
      le = edgIter.LastActiveP();                      //set edgIter to lastactive
   }
   else
   {
      edgeList.insertAtBack( tempEdge1 );        //put edge1 in list
      edgeList.insertAtBack( tempEdge2 );        //put edge2 in list
      le = edgIter.LastP();                            //set edgIter to last
   }

   // Add pointers to the new edges to nodes' spokeLists
   // Three possible cases: (1) there aren't any spokes currently attached,
   // so just put the new one at the front of the list and make it circ'r;
   // (2) there is only one spoke, so it doesn't matter where we attach
   // (3) there is already >1 spoke (the general case)
   spokIter.Reset( node2->getSpokeListNC() );
     //cout << "node " << node2->getID() << ": ";
   if( node2->getSpokeListNC().isEmpty() )
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " in otherwise empty list" << endl;
      node2->insertFrontSpokeList( le);
      node2->getSpokeListNC().makeCircular();
      node2->AttachFirstSpoke( le ); // gt added to update ccwedg 2/99
   }
   else if( spokIter.ReportNextP() == spokIter.DatPtr() )
   {
      node2->insertFrontSpokeList( le);
      ce = node2->getEdg();  // these 2 lines added by gt 2/99
      assert( ce!=0 );
      ce->WelcomeCCWNeighbor( le );
        //node2->getSpokeListNC().makeCircular();
   }
   else // general case: figure out where to attach spoke
   {
      //NB: I (GT) wonder whether we could speed this up. If you knew what
      // triangle you were falling in, wouldn't you be easily able to find
      // the right spoke? TODO -- investigate
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " w/ reference to node " << node3->getID() << endl;
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      //assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
      if( spokIter.AtEnd() )
      {
         //cout << "AddEdge: using new algorithm\n";
         for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
         {
            /*Xf( PointsCCW( UnitVector( ce ),
                           UnitVector( le ),
                           UnitVector( spokIter.ReportNextP() ) ) )*/
            p1 = UnitVector( ce );
            p2 = UnitVector( le );
            p3 = UnitVector( spokIter.ReportNextP() );
            if( PointsCCW( p1, p2, p3 ) )
                break;
         }
      }
      node2->getSpokeListNC().insertAtNext( le,
                                            spokIter.NodePtr() ); //put edge2 in SPOKELIST
      assert( ce!=0 );
      ce->WelcomeCCWNeighbor( le );
   }
   spokIter.Reset( node1->getSpokeListNC() );
   le = edgIter.PrevP();                     //step backward once in edgeList
     //cout << "node " << node1->getID() << ": ";
   if( node1->getSpokeListNC().isEmpty() )
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " in otherwise empty list" << endl;
      node1->insertFrontSpokeList( le );
      node1->getSpokeListNC().makeCircular();
      node1->AttachFirstSpoke( le ); // Tell node it's getting a spoke
   }
   else if( spokIter.ReportNextP() == spokIter.DatPtr() )
   {
      node1->insertFrontSpokeList( le );
        //node2->getSpokeListNC().makeCircular();
      ce = node1->getEdg();  // these 2 lines added by gt 2/99
      assert( ce!=0 );
      ce->WelcomeCCWNeighbor( le ); // Tell node it has a new neighbor!
   }
   else
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " w/ reference to node " << node3->getID() << endl;
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      //assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
      if( spokIter.AtEnd() )
      {
         //cout << "AddEdge: using new algorithm\n";
         for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
         {
            /*Xf( PointsCCW( UnitVector( ce ),
                           UnitVector( le ),
                           UnitVector( spokIter.ReportNextP() ) ) )*/
            p1 = UnitVector( ce );
            p2 = UnitVector( le );
            p3 = UnitVector( spokIter.ReportNextP() );
            if( PointsCCW( p1, p2, p3 ) )
            {
               spokIter.Next();
               break;
            }
         }
      }
      node1->getSpokeListNC().insertAtPrev( le,
                                            spokIter.NodePtr() );//put edge1 in SPOKELIST
      assert( ce!=0 );
      ce->WelcomeCCWNeighbor( le );  // Tell node it has a new neighbor!
   }
   
   nedges+=2;

   // Reset edge id's
   for( ce = edgIter.FirstP(), miNextEdgID = 0; !( edgIter.AtEnd() ); 
        ce = edgIter.NextP(), miNextEdgID++ )
   {
      ce->setID( miNextEdgID );
      /*cout << "    Edg " << i << " (" << ce->getOriginPtr()->getID() << "->"
           << ce->getDestinationPtr()->getID() << ")\n";*/
   }
   //cout << "AddEdge() done\n" << flush;
   return 1;
}


/**************************************************************************\
**
**  tMesh::AddEdgeAndMakeTriangle 
**
**  Function to add the "missing" edge and make
**  the triangle. Formerly more complicated than AddEdge() and
**  MakeTriangle(); now simply calls these functions.
**
**  Inputs: a tPtrList<tSubNode> of nodes in triangle; a tPtrListIter
**   object, the iterator of the latter list; edge is added between node
**   currently pointed to by iterator and the node-after-next. List should
**   be circular.
**  Calls: AddEdge(), MakeTriangle()
**  Created: SL fall, '97
**  Modified: SL 10/98 to call AddEdge() and MakeTriangle()
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
AddEdgeAndMakeTriangle( tPtrList< tSubNode > &nbrList,
                        tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   //cout << "AddEdgeAndMakeTriangle" << endl << flush;
   //Xint i, j, newid;
   tSubNode *cn, *cnn, *cnnn;
   tPtrList< tSubNode > tmpList;
   tPtrListIter< tSubNode > tI( tmpList );

   cn = nbrIter.DatPtr();
   cnn = nbrIter.NextP();
   cnnn = nbrIter.ReportNextP();
   nbrIter.Prev();
   if( !AddEdge( cnnn, cn, cnn ) ) return 0;
   tmpList.insertAtBack( cn );
   tmpList.insertAtBack( cnn );
   tmpList.insertAtBack( cnnn );
   tmpList.makeCircular();
   if( !MakeTriangle( tmpList, tI ) ) return 0;
   return 1;
}


/**************************************************************************\
**  
**  tMesh::MakeTriangle
**
**   Function to make triangle and add it to mesh; called
**   when all necessary nodes and edges already exist, i.e., the triangle
**   exists geometrically but not as a "triangle" member of the data
**   structure. Checks to make sure points are CCW. Resets triangle IDs at
**   end (necessary?). This function is relatively messy and complicated
**   but is extensively commented below.
**
**  Inputs: a tPtrList<tSubNode> of nodes in triangle; a tPtrListIter
**   object, the iterator of the latter list; edge is added between node
**   currently pointed to by iterator and the node-after-next. List must
**   contain three, and only three, members and be circular.
**  Created: SL fall, '97
**
**  Modifications:
**   - mSearchOriginTriPtr is reset to point to the newly added triangle
**     in an attempt to speed up triangle searches, especially during
**     mesh creation. GT 1/2000
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
MakeTriangle( tPtrList< tSubNode > &nbrList,
              tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   assert( nbrList.getSize() == 3 );
   //cout << "MakeTriangle" << endl << flush;
   int i, j;
   //Xint newid;                          // ID of new triangle
   //XtTriangle tempTri;
   tTriangle *nbrtriPtr;
   tSubNode *cn, *cnn, *cnnn;
   tEdge *ce;
   tTriangle *ct;
   tListIter< tTriangle > triIter( triList );
   tMeshListIter< tEdge > edgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   assert( nbrList.getSize() == 3 );

   //Xcn = nbrIter.DatPtr(); Below is bug fix:
   cn = nbrIter.FirstP();      // cn, cnn, and cnnn are the 3 nodes in the tri
   cnn = nbrIter.NextP();
   cnnn = nbrIter.NextP();
   nbrIter.Next();
   tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
       p2( cnnn->get2DCoords() );

   // error check
   if( !PointsCCW( p0, p1, p2 ) )
   {
      cerr << "in MT nodes not CCW: " << cn->getID() << ", "
           << cnn->getID() << ", " << cnnn->getID();
      /*if( cn->Meanders() ) p0 = cn->getNew2DCoords();
      else p0 = cn->get2DCoords();
      if( cnn->Meanders() ) p1 = cnn->getNew2DCoords();
      else p1 = cnn->get2DCoords();
      if( cnnn->Meanders() ) p2 = cnnn->getNew2DCoords();
      else p2 = cnnn->get2DCoords();
      if( !PointsCCW( p0, p1, p2 ) )
      cerr << "; nor are new coords CCW ";*/
      cerr << endl;
   }

   /*cout << "In MT, the 3 nodes are: " << cn->getID() << " " << cnn->getID()
        << " " << cnnn->getID() << endl;*/
   
   //X OBSOLETE: Now uses miNextTriID, GT 1/2000
   // set the ID for the new triangle based on the ID of the last triangle
   // on the list plus one, or if there are no triangles on the list yet
   // (which happens when we're creating an initial "supertriangle" as in
   // MakeMeshFromPoints), set the ID to zero.
   /*ct = triIter.LastP();
   if( ct ) newid = ct->getID() + 1;
   else newid = 0;
   tempTri.setID( newid );*/

   /* The following block is made obsolete through the use of a new
      triangle constructor that automatically initializes vertex and
      edge pointers, 1/2000 GT */

/*X   // set edge and vertex ptrs & add to triList: We go through each point,
   // p0, p1, and p2. At each step, we assign p(i) to the triangle's
   // corresponding pPtr(i), and get the spoke list for node p(i). We then
   // advance such that cn points to p(i+1) and cnn points to p(i+2), and
   // find the edge that connects p(i) and p(i+2). This edge is the clw
   // edge #i for the triangle, so we assign it as such. After this loop is
   // done, tempTri points to the three vertices and to the three clockwise-
   // oriented edges p0->p2 (e0), p1->p0 (e1), and p2->p1 (e2). The nbr
   // triangle pointers are initialized to zero.
   for( i=0; i<3; i++ )
   {
      tempTri.setPPtr( i, cn );       //set tri VERTEX ptr i
      //cout << "cn: " << cn->getID() << endl;
      spokIter.Reset( cn->getSpokeListNC() );
      cn = nbrIter.NextP();           //step forward once in nbrList to p(i+1)
      cnn = nbrIter.ReportNextP();    //get p(i+2) (but don't step forward)

      // Find edge p(i)->p((i+2)%3)
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );

      assert( !( spokIter.AtEnd() ) );

      // Assign edge p(i)->p(i+2) as the triangle's clockwise edge #i
      // (eg, ePtr(0) is the edge that connects p0->p2, etc)
      tempTri.setEPtr( i, ce );      //set tri EDGE ptr i
      tempTri.setTPtr( i, 0 );       //initialize tri TRI ptrs to nul
      }*/

   // Create the new triangle and insert a pointer to it on the list.
   // Here, the triangle constructor takes care of setting pointers to
   // the 3 vertices and 3 edges. The neighboring triangle pointers are
   // initialized to zero.
   triList.insertAtBack( tTriangle( miNextTriID++, cn, cnn, cnnn ) );//put 

   ct = triIter.LastP();            //ct now points to our new triangle
   assert( cn == ct->pPtr(0) );     //make sure we're where we think we are

   // To speed up future searches in calls to LocateTriangle, assign the
   // starting triangle, mSearchOriginTriPtr, to the our new triangle.
   // The idea here is that there's a good chance that the next point
   // to be added will be close to the current location. (added 1/2000)
   mSearchOriginTriPtr = ct;
   
   /*cout << "IN MT, created triangle:\n";
   ct->TellAll();*/
   
   // Now we assign the neighbor triangle pointers. The loop successively
   // gets the spokelist for (p0,p1,p2) and sets cn to the next ccw point
   // (p1,p2,p0). It then finds the edge (spoke) that joins the two points
   // (p0->p1, p1->p2, p2->p0). These are the edges that are shared with
   // neighboring triangles (t2,t0,t1) and are pointed to by the neighboring
   // triangles. This means that in order to find neighboring triangle t2,
   // we need to find the triangle that points to edge (p0->p1), and so on.
   // In general, t((j+2)%3) is the triangle that points to edge
   // p(j)->p((j+1)%3).
   //Xdce = 0;
   nbrtriPtr = 0;
   cn = nbrIter.FirstP();
   //cout << "starting w/ node " << cn->getID();
   for( j=0; j<3; j++ )
   {
      // get spokelist for p(j) and advance to p(j+1)
      spokIter.Reset( cn->getSpokeListNC() );
      cn = nbrIter.NextP();               //step forward once in nbrList
      //Xif( j>0 ) dce = ce;

      // Find edge ce that connects p(j)->p(j+1)
      for( ce = spokIter.FirstP();
           ce->getDestinationPtrNC() != cn && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );
        //********BUG: following assertion failed; called from FlipEdge,
        //from CheckForFlip, from CheckLocallyDelaunay, from MoveNodes***************
      if( !( TriWithEdgePtr( ce ) != nbrtriPtr || nbrtriPtr == 0 ) )
      {
         p0 = cn->get2DCoords();
         p1 = cnn->get2DCoords();
         p2 = cnnn->get2DCoords();
         
         if( PointsCCW( p0, p1, p2 ) )
             cerr << "something FUNNY going on";
         else cerr << "tri not CCW: " << nbrtriPtr->getID() << endl;
      }

      // Find the triangle, if any, that shares (points to) this edge
      // and assign it as the neighbor triangle t((j+2)%3).
      nbrtriPtr = TriWithEdgePtr( ce );
      /*cout << "The following tri shares edg the following edge:\n";
      ce->TellCoords();
      if( nbrtriPtr )
          nbrtriPtr->TellAll();
      else cout << "(none)\n";*/
         
      ct->setTPtr( (j+2)%3, nbrtriPtr );      //set tri TRI ptr (j+2)%3

      //cout << "This is our nbr #" << (j+2)%3 << endl << endl;
      
      // If a neighboring triangle was found, tell it that the current
      // new triangle is its neighbor too. We need to tell it which
      // neighbor we are (0, 1, or 2), and the mapping is like this:
      // if the nbr tri calls the shared edge (0,1,2) then we are its
      // nbr (1,2,0). (ie, tri_number = (edg_number+1)%3 )
      if( nbrtriPtr != 0 )
      {
         for( i=0; i<3; i++ )
         {
            assert( nbrtriPtr->ePtr(i) != 0 );
            assert( ce != 0 );
            if( nbrtriPtr->ePtr(i) == ce ) break;
         }
         assert( i < 3 );
         nbrtriPtr->setTPtr( (i+1)%3, ct );  //set NBR TRI ptr to tri
      }
   }   
   ntri++;
   
   //reset triangle id's (why needed??) because when we make a new item of any kind we
     //give it an id; how do we know what id to use (i.e., what's large enough but not
     //too large)? we find the id of the last item in the list and add one; if the items
     //in the list have been "mixed up", then we could assign an id already in use;
     //also, if for some reason numbers are systematically skipped, the id could blow up;
     //this step
     //may not be strictly necessary for triangles (it is for nodes and edges where
     //we have active and inactive members), but I'm sure it doesn't hurt; better safe
     //than sorry...
   for( ct = triIter.FirstP(), miNextTriID=0; !( triIter.AtEnd() );
        ct = triIter.NextP(), miNextTriID++ )
   {
      ct->setID( miNextTriID );
   }
   
   return 1;
}


/**************************************************************************\
**
**   tMesh::AddNode ( tSubNode nodeRef& )
**
**   Adds a new node with the properties of nodRef to the mesh.
**
**   Calls: tMesh::LocateTriangle, tMesh::DeleteTriangle, tMesh::AddEdge,
**            tMesh::AddEdgeAndMakeTriangle, tMesh::MakeTriangle,
**            tMesh::CheckForFlip; various member functions of tNode,
**            tMeshList, tMeshListIter, tPtrList, etc. Also tLNode
**            functions (TODO: this needs to be removed somehow),
**            and temporarily, tMesh::UpdateMesh
**   Parameters: nodeRef -- reference to node to be added (really,
**                          duplicated)
**   Returns:  (always TRUE: TODO make void return type)
**   Assumes:
**   Created: SL fall, '97
**   Modifications:
**        - 4/98: node is no longer assumed to be a non-boundary (GT)
**        - 7/98: changed return type from int (0 or 1) to ptr to
**                the new node (GT)
**        -10/98: if node is open boundary,
**                added with tMeshList::insertAtBoundFront() (SL)
**        -5/99: removed unreferenced vars tedg1, tedg3 (GT)
**
\**************************************************************************/
#define kLargeNumber 1000000000
template< class tSubNode >
tSubNode * tMesh< tSubNode >::
AddNode( tSubNode &nodeRef, int updatemesh, double time )
{
   int i, ctr;
   tTriangle *tri;
   tSubNode *cn;
   tArray< double > xyz( nodeRef.get3DCoords() );
   tMeshListIter< tSubNode > nodIter( nodeList );
   assert( &nodeRef != 0 );

   //cout << "AddNode at " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << " time "<<time<<endl;

   //cout << "locate tri & layer interp" << endl << flush;
   tri = LocateTriangle( xyz[0], xyz[1] );
   //assert( tri != 0 );
   if( tri == 0 )
   {
      cerr << "tMesh::AddNode(...): node coords out of bounds: " << xyz[0] << " " << xyz[1] << endl;
      return 0;
   }
   if( layerflag && time > 0 ) 
       nodeRef.LayerInterpolation( tri, xyz[0], xyz[1], time );
   
   // Assign ID to the new node and insert it at the back of either the active
   // portion of the node list (if it's not a boundary) or the boundary
   // portion (if it is)
   //cout << "inserting node on list\n";
   //Xint newid = nodIter.LastP()->getID() + 1;
   nodeRef.setID( miNextNodeID );
   miNextNodeID++;
   if( nodeRef.getBoundaryFlag()==kNonBoundary ){
       nodeList.insertAtActiveBack( nodeRef );
       //cout<<"knonboundary"<<endl;
   }
   else if( nodeRef.getBoundaryFlag() == kOpenBoundary ){
       nodeList.insertAtBoundFront( nodeRef );
       //cout<<"kOpenBoundary"<<endl;
   }
   else{
       nodeList.insertAtBack( nodeRef );
       //cout<<"other"<<endl;
   }

   //assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
   
   // Retrieve a pointer to the new node and flush its spoke list
   if( nodeRef.getBoundaryFlag() == kNonBoundary )
       cn = nodIter.LastActiveP();
   else if( nodeRef.getBoundaryFlag() == kOpenBoundary )
       cn = nodIter.FirstBoundaryP();
   else
       cn = nodIter.LastP();
   assert( cn!=0 );
   cn->getSpokeListNC().Flush();
   
   //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = static_cast<tSubNode *>(tri->pPtr(i));
      bndyList.insertAtBack( tmpPtr );
   }
   bndyList.makeCircular();


   // Delete the triangle in which the node falls
   i = DeleteTriangle( tri );
   assert( i != 0 );  //if ( !DeleteTriangle( tri ) ) return 0;

   //make 3 new triangles
   //cout << "creating new triangles\n" << flush;
   tPtrListIter< tSubNode > bndyIter( bndyList );
   tSubNode *node3 = bndyIter.FirstP();     // p0 in original triangle
   //XtSubNode *node2 = nodIter.LastActiveP(); // new node
   tSubNode *node2 = cn;                    // new node
   tSubNode *node1 = bndyIter.NextP();      // p1 in orig triangle
   tSubNode *node4 = bndyIter.NextP();      // p2 in orig triangle
   tArray< double > p1( node1->get2DCoords() ),
       p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
       p4( node4->get2DCoords() );

   if( xyz.getSize() == 3) //why would this ever not be the case? If we need to access new coords:
       //size of xyz is basically the flag; the 4th element is never used o.w.
   {
      //cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
      //     << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   else
   {
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) p3 = node3->getNew2DCoords();
      if( node4->Meanders() ) p4 = node4->getNew2DCoords();  
      /*cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
           << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;*/
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }

   // Here's how the following works. Let the old triangle vertices be A,B,C
   // and the new node N. The task is to create 3 new triangles ABN, NBC, and
   // NCA, and 3 new edge-pairs AN, BN, and CN.
   // First, edge pair BN is added. Then AEMT is called to create triangle
   // ABN and edge pair AN. AEMT is called again to create tri NBC and edge
   // pair CN. With all the edge pairs created, it remains only to call
   // MakeTriangle to create tri NCA.
   //cout << "calling AE, AEMT, AEMT, and MT\n" << flush;
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   AddEdge( node1, node2, node3 );  //add edge between node1 and node2
   tPtrList< tSubNode > tmpList;
   tmpList.insertAtBack( node3 );  // ABN
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node2 );
   tPtrListIter< tSubNode > tmpIter( tmpList );
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );  // NBC
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node4 );
   tmpIter.First();
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );  // NCA
   tmpList.insertAtBack( node4 );
   tmpList.insertAtBack( node3 );
   tmpList.makeCircular();
   tmpIter.First();
   MakeTriangle( tmpList, tmpIter );

   //hasn't changed yet
   //put 3 resulting triangles in ptr list
   //cout << "Putting tri's on list\n" << flush;
   if( xyz.getSize() == 3 )
   {
      //Xcout << "flip checking in addnode" << endl;
      tPtrList< tTriangle > triptrList;
      tListIter< tTriangle > triIter( triList );
      tPtrListIter< tTriangle > triptrIter( triptrList );
      tTriangle *ct;
      triptrList.insertAtBack( triIter.LastP() );
      triptrList.insertAtBack( triIter.PrevP() );
      triptrList.insertAtBack( triIter.PrevP() );

      //check list for flips; if flip, put new triangles at end of list
      int flip = 1;
      ctr = 0;
      while( !( triptrList.isEmpty() ) )
      {
         ctr++;
         if( ctr > kLargeNumber ) // Make sure to prevent endless loops
         {                        // TODO: remove for release ver
            cerr << "Mesh error: adding node " << node2->getID()
                 << " flip checking forever"
                 << endl;
            ReportFatalError( "Bailing out of AddNode()" );
         }
         ct = triptrIter.FirstP();
         
         for( i=0; i<3; i++ )
         {
            if( ct->tPtr(i) != 0 )
            {
               if( CheckForFlip( ct, i, flip ) )
               {
                  triptrList.insertAtBack( triIter.LastP() );
                  triptrList.insertAtBack( triIter.PrevP() );
                  break;
               }
            }
         }

         triptrList.removeFromFront( ct );
      }
   }

   //reset node id's
   //cout << "resetting ids\n" << flush;
   for( cn = nodIter.FirstP(), miNextNodeID=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), miNextNodeID++ )
   {
      cn->setID( miNextNodeID );
   }
   node2->makeCCWEdges();
   node2->InitializeNode();

   if( updatemesh ) UpdateMesh();
   //cout << "AddNode finished" << endl;

   tEdge *ce, *fe;
   //fe = node2->getFlowEdg(); GT changed to Edg 'cus was crashing
   fe = node2->getEdg();
   ce = fe;

   //Xint hlp=0;
   //do{
   // ce=ce->getCCWEdg();
   // hlp++;
   //}while(ce != fe );
   //  if(hlp !=  node2->getSpokeListNC().getSize() ){
//        cout<<"AddNode  number of spokes "<<node2->getSpokeListNC().getSize()<<" number of ccwedges "<<hlp<<endl<<flush;
//     }

   return node2;  // Return ptr to new node
}


/**************************************************************************\
**
**  tMesh::AddNodeAt
**
**   add a node with referenced coordinates to mesh;
**   this fn duplicates functionality of AddNode
**
**  Created: SL fall, '97
**  Modified: NG summer, '98 to deal with layer interpolation
**
\**************************************************************************/
//TODO: ; just assign coords
// to a dummy new node and call AddNode
template< class tSubNode >
tSubNode *tMesh< tSubNode >::
AddNodeAt( tArray< double > &xyz, double time )
{
   assert( &xyz != 0 );
   // cout << "AddNodeAt " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] <<" time "<<time<< endl;
   tTriangle *tri;
   //cout << "locate tri" << endl << flush;
   if( xyz.getSize() == 3 ) tri = LocateTriangle( xyz[0], xyz[1] );
   else tri = LocateNewTriangle( xyz[0], xyz[1] );
   if( tri == 0 )      return 0;
   
   int i, ctr;
   tMeshListIter< tSubNode > nodIter( nodeList );
   tSubNode tempNode, *cn;
   tempNode.set3DCoords( xyz[0], xyz[1], xyz[2]  );
   if( layerflag && time > 0.0) tempNode.LayerInterpolation( tri, xyz[0], xyz[1], time );
   if( xyz.getSize() != 3 ) tempNode.setNew2DCoords( xyz[0], xyz[1] );
   tempNode.setBoundaryFlag( 0 );

   // Assign ID to the new node and insert it at the back of the active
   // portion of the node list (NOTE: node is assumed NOT to be a boundary)
   //Xint newid = nodIter.LastP()->getID() + 1;
   tempNode.setID( miNextNodeID );
   miNextNodeID++;
   cout << miNextNodeID << endl;
   nodeList.insertAtActiveBack( tempNode );
   assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
   
     //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = static_cast<tSubNode *>(tri->pPtr(i));
      bndyList.insertAtBack( tmpPtr );
   }
   bndyList.makeCircular();
   //delete triangle
   //Xcout << "calling deletetri from addnodeat\n";
   if ( !DeleteTriangle( tri ) ) return 0;
   //make 3 new triangles
   tPtrListIter< tSubNode > bndyIter( bndyList );
   tSubNode *node3 = bndyIter.FirstP();
   tSubNode *node2 = nodIter.LastActiveP();
   tSubNode *node1 = bndyIter.NextP();
   tSubNode *node4 = bndyIter.NextP();
   tArray< double > p1( node1->get2DCoords() ),
       p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
       p4( node4->get2DCoords() );
   if( xyz.getSize() == 3)
   {
      //cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
      // << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   else
   {
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) p3 = node3->getNew2DCoords();
      if( node4->Meanders() ) p4 = node4->getNew2DCoords();  
      //cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
      //   << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) || !PointsCCW( p2, p1, p4 ) || !PointsCCW( p2, p4, p3 ) )
          cout << "new tri not CCW" << endl;
   }
   
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   AddEdge( node1, node2, node3 );  //add edge between node1 and node2
   tPtrList< tSubNode > tmpList;
   tmpList.insertAtBack( node3 );
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node2 );
   tPtrListIter< tSubNode > tmpIter( tmpList );
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );
   tmpList.insertAtBack( node1 );
   tmpList.insertAtBack( node4 );
   tmpIter.First();
   AddEdgeAndMakeTriangle( tmpList, tmpIter );
   tmpList.Flush();
   tmpList.insertAtBack( node2 );
   tmpList.insertAtBack( node4 );
   tmpList.insertAtBack( node3 );
   tmpList.makeCircular();
   tmpIter.First();
   MakeTriangle( tmpList, tmpIter );
   //put 3 resulting triangles in ptr list
   if( xyz.getSize() == 3 )
   {
      //cout << "flip checking" << endl;
      tPtrList< tTriangle > triptrList;
      tListIter< tTriangle > triIter( triList );
      tPtrListIter< tTriangle > triptrIter( triptrList );
      tTriangle *ct;
      triptrList.insertAtBack( triIter.LastP() );
      triptrList.insertAtBack( triIter.PrevP() );
      triptrList.insertAtBack( triIter.PrevP() );
        //check list for flips; if flip, put new triangles at end of list
      int flip = 1;
      ctr = 0;
      while( !( triptrList.isEmpty() ) )
      {
         ctr++;
         if( ctr > kLargeNumber ) // Make sure to prevent endless loops
         {
            cerr << "Mesh error: adding node " << node2->getID()
                 << " flip checking forever"
                 << endl;
            ReportFatalError( "Bailing out of AddNodeAt()" );
         }
         ct = triptrIter.FirstP();
         for( i=0; i<3; i++ )
         {
            if( ct->tPtr(i) != 0 )
            {
               if( CheckForFlip( ct, i, flip ) )
               {
                  triptrList.insertAtBack( triIter.LastP() );
                  triptrList.insertAtBack( triIter.PrevP() );
                  break;
               }
            }
         }
         triptrList.removeFromFront( ct );
      }
   }
   //reset node id's
   cout << "reset ids\n";
   for( cn = nodIter.FirstP(), miNextNodeID=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), miNextNodeID++ )
   {
      cn->setID( miNextNodeID );
   }
   //nmg uncommented line below and added initialize line
   node2->makeCCWEdges();
   node2->InitializeNode();
   
   UpdateMesh();

   tEdge *ce, *fe;
   fe = node2->getFlowEdg();
   ce = fe;

   int hlp=0;
   do{
      ce=ce->getCCWEdg();
      hlp++;
   }while(ce != fe );
   //  if(hlp !=  node2->getSpokeListNC().getSize() ){
//        cout<<"AddNodeAt  number of spokes "<<node2->getSpokeListNC().getSize()<<" number of ccwedges "<<hlp<<endl<<flush;
//     }

   //cout << "AddNodeAt finished, " << nnodes << endl;
   return node2;
}
#undef kLargeNumber


/**************************************************************************\
**
**  tMesh "get" functions
**
\**************************************************************************/

template <class tSubNode>
tMeshList<tEdge> * tMesh<tSubNode>::
getEdgeList() {return &edgeList;}

template <class tSubNode>
tMeshList<tSubNode> * tMesh<tSubNode>::
getNodeList() {return &nodeList;}

template <class tSubNode>
tList< tTriangle > * tMesh<tSubNode>::
getTriList() {return &triList;}


/**************************************************************************\
**
**  tMesh::getEdgeComplement
**
**  Returns the complement of _edge_ (i.e., the edge that shares the same
**  endpoints but points in the opposite direction). To find the complement,
**  it exploits the fact that complementary pairs of edges are stored 
**  together on the edge list, with the first of each pair having an 
**  even-numbered ID and the second having an odd-numbered ID.
**
**  Modifications: gt replaced 2nd IF with ELSE to avoid compiler warning
**
\**************************************************************************/
template< class tSubNode >
tEdge *tMesh< tSubNode >::
getEdgeComplement( tEdge *edge )
{
   tMeshListIter< tEdge > edgIter( edgeList );
   int edgid = edge->getID();

   assert( edgIter.Get( edgid ) );
   edgIter.Get( edgid ); // TODO: why necessary?
   if( edgid%2 == 0 ) return edgIter.GetP( edgid + 1 );
   else return edgIter.GetP( edgid - 1 );
}


/**************************************************************************\
**
**  tMesh::UpdateMesh
**
**  Updates mesh geometry:
**   - computes edge lengths
**   - finds Voronoi vertices
**   - computes Voronoi edge lengths
**   - computes Voronoi areas for interior (active) nodes
**   - updates CCW-edge connectivity
**
**  Note that the call to CheckMeshConsistency is for debugging
**  purposes and should be removed prior to release.
**
**  Calls: MakeCCWEdges(), setVoronoiVertices(), CalcVoronoiEdgeLengths(),
**   CalcVAreas(), CheckMeshConsistency()
**  Assumes: nodes have been properly triangulated
**  Created: SL fall, '97
**
\**************************************************************************/
template <class tSubNode>
void tMesh<tSubNode>::
UpdateMesh()
{
   //cout << "UpdateMesh()" << endl << flush;
   
   //tListIter<tTriangle> tlist( triList );
   tMeshListIter<tEdge> elist( edgeList );
   //tMeshListIter<tSubNode> nlist( nodeList );
   tEdge * curedg = 0;
   double len;
   
   // Edge lengths
   curedg = elist.FirstP();
   do
   {
      len = curedg->CalcLength();
      if( len<=0.0 ) {
	cout << "Edge " << curedg->getID() << " length: " << curedg->getLength() << endl;
	curedg->TellCoords();
      }
      assert( len>0.0 );
      curedg = elist.NextP();
      assert( curedg > 0 ); // failure = complementary edges not consecutive
      curedg->setLength( len );
   } while( curedg=elist.NextP() );

   MakeCCWEdges();

   setVoronoiVertices();
   CalcVoronoiEdgeLengths();
   CalcVAreas();
   CheckMeshConsistency( 0 );  // debug only -- remove for release

// Triangle areas
/*   for( tlist.First(); !tlist.AtEnd(); tlist.Next() )
   {
      curtri = tlist.DatPtr();
      curtri->length_sides();
      curtri->CalcArea();
      curtri = curtri->next;
   }
   */
}


/*****************************************************************************\
**
**  tMesh::CheckForFlip
**
**  Checks whether edge between two triangles should be
**  flipped; may either check, flip, and report, or just check and report.
**  Checks whether the present angle or the possible angle
**  is greater. Greater angle wins. Also uses flip variable
**  to determine whether to use newx, newy, or x, y.
**
**      Inputs: tri -- ptr to the triangle to be tested
**              nv -- the number of the vertex opposite the edge that
**                    might be flipped (0, 1, or 2)
**              flip -- flag indicating whether we want to actually flip
**                      the edge if needed (TRUE) or simply test the flip
**                      condition for a point that is about to be moved to
**                      a new position (FALSE)
**      Returns: 1 if flip is needed, 0 otherwise
**      Modifies: edge may be flipped
**      Called by: AddNode, AddNodeAt, CheckLocallyDelaunay,
**                 tStreamMeander::CheckBrokenFlowedge
**      Calls: PointsCCW, FlipEdge, TriPasses
**        
**      Created: 8/28/97 SL
**      Modified: 12/16/97 SL                                               
**                                                  
\*****************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
CheckForFlip( tTriangle * tri, int nv, int flip )
{
   if( tri == 0 )  // TODO: is this just a bug check?
   {
      cout << "CheckForFlip: tri == 0" << endl;
      return 0;
   }
   assert( nv < 3 );
   //cout << "THIS IS CheckForFlip(...) " << tri->getID() << endl;
   tSubNode *node0, *node1, *node2, *node3;
   node0 = static_cast< tSubNode * >(tri->pPtr(nv));
   //cout<<"node0 id "<<node0->getID()<<endl;
   node1 = static_cast< tSubNode * >(tri->pPtr((nv+1)%3));
   //cout<<"node1 id "<<node1->getID()<<endl;
   node2 = static_cast< tSubNode * >(tri->pPtr((nv+2)%3));
   //cout<<"node2 id "<<node2->getID()<<endl;
   tTriangle *triop = tri->tPtr(nv);
   int nvop = triop->nVOp( tri );
   node3 = static_cast< tSubNode * >(triop->pPtr( nvop ));
   tArray< double > ptest( node3->get2DCoords() ), p0( node0->get2DCoords() ),
       p1( node1->get2DCoords() ), p2( node2->get2DCoords() );

   // If "flip" flag isn't set and the node is a moving node, use "new"
   // coordinates rather than current coordinates
   // TODO: decouple this from meandering -- use a "moving" flag instead
   // for generality?
   if( !flip )
   {
      if( node0->Meanders() ) p0 = node0->getNew2DCoords();
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) ptest = node3->getNew2DCoords();
   }

   // If p0-p1-p2 passes the test, no flip is necessary
   if( TriPasses( ptest, p0, p1, p2 ) ) return 0;

   // Otherwise, a flip is needed, provided that the new triangles are
   // counter-clockwise (would this really ever happen??) and that the
   // node isn't a moving node (ie "flip" is true)
   if( flip )                     //and make sure there isn't already an edge?
   {
      if( !PointsCCW( p0, p1, ptest ) || !PointsCCW( p0, ptest, p2 ) ) 
          return 0;
      //cout << "calling Flip edge from cff" << endl;
      FlipEdge( tri, triop, nv, nvop );
        /*
      assert( DeleteEdge( tri->ePtr( (nv+1)%3 ) ) );
      tPtrList< tSubNode > nbrList;
      nbrList.insertAtBack( node0 );
      nbrList.insertAtBack( node1 );
      nbrList.insertAtBack( node3 );
      nbrList.insertAtBack( node2 );
      assert( RepairMesh( nbrList ) );*/
   }
     //cout << "finished" << endl;
   return 1;
}


/******************************************************************\
**
**  tMesh::FlipEdge
**
**  Flips the edge pair between two adjacent triangle to
**  re-establish Delaunay-ness.
**
**  Note on notation in flip edge:
**
**                d
**               /|\
**       tri->  / | \ <-triop
**             /  |  \
**            a   |   c
**             \  |  /
**              \ | /
**               \|/
**                b
**        Edge bd will be removed
**        and an edge ac will be made.
**        nbrList contains the points a, b, c, d
**
**    Inputs:  tri, triop -- the triangles sharing the edge to be
**                           flipped
**             nv -- the number of tri's vertex (0, 1 or 2) opposite
**                   the edge (ie, point a)
**             nvop -- the number of triop's vertex (0, 1 or 2)
**                     opposite the edge (ie, point c)
**    Calls: DeleteEdge, AddEdgeAndMakeTriangle, MakeTriangle
**    Called by: CheckForFlip, CheckTriEdgeIntersect
**
\*******************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
FlipEdge( tTriangle * tri, tTriangle * triop ,int nv, int nvop )
{
   //cout << "FlipEdge(...)..." << endl;
   tSubNode *cn = 0;
   tPtrList< tSubNode > nbrList;
   //DumpTriangles();

   // Place the four vertices of the two triangles on a list
   nbrList.insertAtBack( static_cast<tSubNode *>(tri->pPtr(nv)) );
   nbrList.insertAtBack( static_cast<tSubNode *>(tri->pPtr((nv+1)%3)) );
   nbrList.insertAtBack( static_cast<tSubNode *>(triop->pPtr( nvop )) );
   nbrList.insertAtBack( static_cast<tSubNode *>(tri->pPtr((nv+2)%3)) );
   nbrList.makeCircular();

   // Delete the edge pair between the triangles, along with the tri's
   //cout << "calling deleteedge from flipedge\n";
   //XDeleteEdge( tri->ePtr( (nv+1)%3 ) );
   DeleteEdge( tri->ePtr( (nv+2)%3 ) );  // Changed for right-hand data struc

   // Recreate the triangles and the edges in their new orientation
   tPtrListIter< tSubNode > nbrIter( nbrList );
   AddEdgeAndMakeTriangle( nbrList, nbrIter );
   nbrIter.First();
   nbrList.removeNext( cn, nbrIter.NodePtr() );
   MakeTriangle( nbrList, nbrIter );
   //cout << "finished" << endl;

   // TODO: why not just change the endpoints of the edges rather than
   // deleting and recreating them?
}


/*****************************************************************************\
**
**  tMesh::CheckLocallyDelaunay
**
**  Updates the triangulation after moving some points.
**  Only uses x and y values, which have already been updated in
**  MoveNodes (frmr PreApply).
**  MoveNodes SHOULD BE CALLED BEFORE THIS FUNCTION IS CALLED
**
**  The logic here is somewhat complicated. Here is GT's understanding
**  of it (Stephen, can you confirm?):
**
**  1. We create a list of triangles that have at least one vertex that has
**     moved (triPtrList) and which therefore might no longer be
**     Delaunay.
**  2. For each of these, we do a flip check across each face. Before
**     doing so, however, we find the triangle on the triPtrList, if any,
**     that comes just before this neighboring triangle. If the edge between
**     the triangles gets flipped, both the triangles will be deleted and
**     recreated on the master triangle list; thus, we will need to delete
**     both affected triangles from triPtrList and re-add the new ones.
**  3. If a flip occurs, remove the opposite triangle pointer from the
**     list if needed in order to prevent a dangling pointer. The two
**     affected triangles will have been replaced by new triangles which
**     are now at the back of the master triangle list; add these two to
**     the triPtrList to be rechecked, and break out of the vertex loop.
**  4. Remove the triangle in question from the head of the triPtrList
**     (regardless of whether it was flipped or not; if it was, its a 
**     dangling pointer; if not, it is Delaunay and we no longer need
**     worry about it)
**  5. Continue until there are no more triangles to be checked.
**
**      Data members updated: Mesh
**      Called by: MoveNodes
**      Calls: CheckForFlip
**      Created: SL fall, '97
**        
\*****************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
CheckLocallyDelaunay()
{
   //cout << "CheckLocallyDelaunay()" << endl;
   tTriangle *at;
   tPtrList< tTriangle > triPtrList;
   tPtrListIter< tTriangle > triPtrIter( triPtrList );
   tListIter< tTriangle > triIter( triList );
   int i, change;
   //Xint id0, id1, id2;
   tArray< int > npop(3);
   tSubNode *nodPtr;
   int flip = 1;

   // Search through tri list to find triangles with at least one
   // moving vertex, and put these on triPtrList
   //put each triangle into the stack
   for( at = triIter.FirstP(); !( triIter.AtEnd() ); at = triIter.NextP() )
   {
      change = FALSE;
      for( i = 0; i < 3; i++ )
      {
         nodPtr = static_cast< tSubNode * >(at->pPtr(i));
         if( nodPtr->Meanders() ) change = TRUE;
      }
      if( change ) triPtrList.insertAtBack( at );
   }

   // Check list for flips; if flip, put new triangles at end of list
   tPtrListIter< tTriangle > duptriPtrIter( triPtrList );
   tTriangle *tn, *tp;
   while( !( triPtrList.isEmpty() ) ) // keep going 'til we've done em all
   {
      at = triPtrIter.FirstP();
      for( i=0; i<3; i++ )
      {
         // If a neighboring triangle exists across this face, check for flip
         if( at->tPtr(i) != 0 )
         {
            tp = at->tPtr(i);

            // If the neighboring triangle is also on the triPtrList, find
            // the item on the list that comes just before it. This is done
            // so that we can avoid a dangling pointer to a non-existent
            // triangle if the edge is flipped (in which case, both tri's
            // are deleted and new ones created). At the end of the for loop
            // duptriPtrIter will either point to the item just before the
            // neighbor triangle, or it will point to the last item;
            // if tn is nonzero, it means the former is true.
            for( tn = duptriPtrIter.FirstP();
                 duptriPtrIter.ReportNextP() != tp &&
                     !( duptriPtrIter.AtEnd() );
                 tn = duptriPtrIter.NextP() );
            tn = 0;
            if( !( duptriPtrIter.AtEnd() ) )
            {
               // doesn't this just mean tn == tp?? seems that tn just
               // acts as a flag here
               tn = duptriPtrIter.ReportNextP();
            }
            
            /*if( at->tPtr(0) != 0 ) id0 = at->tPtr(0)->getID();
            else id0 = -1;
            if( at->tPtr(1) != 0 ) id1 = at->tPtr(1)->getID();
            else id1 = -1;
            if( at->tPtr(2) != 0 ) id2 = at->tPtr(2)->getID();
            else id2 = -1;
            cout << "check tri " << at->getID() << " with nbrs "
                 << id0 << ", " << id1
                 << ", and " << id2;  
                 if( tp->tPtr(0) != 0 ) id0 = tp->tPtr(0)->getID();
            else id0 = -1;
            if( tp->tPtr(1) != 0 ) id1 = tp->tPtr(1)->getID();
            else id1 = -1;
            if( tp->tPtr(2) != 0 ) id2 = tp->tPtr(2)->getID();
            else id2 = -1;
            cout << " against tri " << tp->getID() << " with nbrs "
                 << id0 << ", " << id1
                 << ", and " << id2 << endl;*/
            //cout << "call cff from cld\n";

            // Check triangle _at_ for a flip across face opposite vertex i,
            // and do the flip if needed
            if( CheckForFlip( at, i, flip ) )
            {
               // If the neighboring triangle (tp) is also on the list,
               // we'll need to delete the pointer to it, which is now
               // dangling (both triangles having been deleted and
               // recreated during CheckForFlip)
               //cout << "flipped tri's, got tri ";
               if( tn != 0 )
                   triPtrList.removeNext( tn, duptriPtrIter.NodePtr() );

               // Now put the two recreated triangles on the triPtrList
               // We assume that CheckForFlip has put them at the back of
               // the triangle list, so we simply add the last two triangles
               // to triPtrList, then break out of the for loop
               tn = triIter.LastP();

               /*if( tn->tPtr(0) != 0 ) id0 = tn->tPtr(0)->getID();
               else id0 = -1;
               if( tn->tPtr(1) != 0 ) id1 = tn->tPtr(1)->getID();
               else id1 = -1;
               if( tn->tPtr(2) != 0 ) id2 = tn->tPtr(2)->getID();
               else id2 = -1;
               cout << tn->getID() << " with nbrs "
                    << id0 << ", " << id1
                    << ", and " << id2;*/

               triPtrList.insertAtBack( tn );
               tn = triIter.PrevP();
               /*if( tn->tPtr(0) != 0 ) id0 = tn->tPtr(0)->getID();
               else id0 = -1;
               if( tn->tPtr(1) != 0 ) id1 = tn->tPtr(1)->getID();
               else id1 = -1;
               if( tn->tPtr(2) != 0 ) id2 = tn->tPtr(2)->getID();
               else id2 = -1;
               cout << " and tri " << tn->getID() << " with nbrs "
                    << id0 << ", " << id1
                    << ", and " << id2 << endl;*/

               triPtrList.insertAtBack( tn );
               break;
            }
         }
      }
      // Whether or not we did the flip, we remove the triangle from
      // triPtrList (if it wasn't flipped, we no longer need to consider it;
      // if it was, it will have been deleted and its replacement added to
      // the back of the triPtrList)
      triPtrList.removeFromFront( at );
   }
      //for each triangle in the stack
/*      for( at = triPtrIter.FirstP(); !( triPtrIter.AtEnd() );
           at = triPtrIter.NextP() )
      {
         for( i=0; i<3; i++ )
         {
            trop[i] = at->tPtr(i);
            if( trop[i] ) npop[i] = trop[i]->nVOp( at );
            else npop[i] = NULL;
         }
         for( i=0; i<3; i++ )
         {
            if( CheckForFlip( trop[i], npop[i], flip ) )
            {
               flipped = TRUE;
            }
         }
      }
   } while( flipped );*/
   //cout << "finished" << endl;
}


/*****************************************************************************\
**
**  tMesh::IntersectsAnyEdge
**
**  Returns the first edge in the list which intersects
**  "edge" or NULL if "edge" intersects no other edges
**
**      Data members updated: Mesh
**      Called by: APPARENTLY NEVER CALLED! Replaced by global
**                 function IntersectsAnyEdgeInList. Retain code in case
**                 it needs to be revived.
**      Calls: Intersect
**      Created: SL fall, '97
**
\*****************************************************************************/
/*template< class tSubNode >
tEdge *tMesh< tSubNode >::
IntersectsAnyEdge( tEdge * edge )
{
   //cout << "IntersectsAnyEdge( tEdge * edge )..." << endl;
   int i;
   tEdge * ce;
   tMeshListIter< tEdge > edgIter( edgeList );

   if( !edge )
   {
      cout<<"IntersectsAnyEdge: Warning: invalid edge"<<endl<<flush;
      return( NULL );
   }
     //cout<<"IAE: nedges "<<nedges<<endl<<flush;
     //cout << "call Intersect for edges " << edge->getID()
     //   << " from nodes " << edge->getOriginPtr()->getID()
     //   << " to " << edge->getDestinationPtr()->getID() << "; " << endl;

   // For every other edge on the list, call Intersect to test
   for( ce = edgIter.FirstP(); !(edgIter.AtEnd());
        edgIter.Next(), ce = edgIter.NextP() )
   {
      assert( edgIter.NodePtr()->getNext() != 0 );
      if( edge->getID() != ce->getID() &&
          edge->getID() != getEdgeComplement( edge )->getID() )
      {
           //cout  << " and " << ce->getID() << " from nodes "
           //    << ce->getOriginPtr()->getID()
           //    << " to " << ce->getDestinationPtr()->getID() << endl;
         if( Intersect( edge, ce ) ) return( ce );
      }
      
   }
   assert( edgIter.AtEnd() );
   //if( i < nedges - 1 )
   //  cout<<"IntersectsAnyEdge: Warning: whole list not checked"<<endl<<flush;
   return( NULL );
}*/


/*****************************************************************************\
**
**  tMesh::CheckTriEdgeIntersect
**
**        This function implements node movement.
**        We want to know if the moving point has passed beyond the polygon
**        defined by its spoke edges; if it has, then we will have edges
**        intersecting one another. In the case where the point has simply
**        passed into one of the 'opposite' triangles, then we can just do a
**        flip operation. In the other case, the remedial action is much more
**        complicated, so we just delete the point and add it again.
**
**  Created: SL fall, '97
**  Modifications:
**   - minor change from i = AddNode to cn = AddNode to handle changed
**     return type (GT 7/98)
**
\*****************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
CheckTriEdgeIntersect()
{
   //cout << "CheckTriEdgeIntersect()..." << flush << endl;
     //DumpNodes();
   int i, j, nv, nvopp;
   int flipped = TRUE;
   int crossed;
   tSubNode *subnodePtr, tempNode, newNode;  
   tEdge * cedg, *ce;
   tTriangle * ct, * ctop, *rmtri/* *tri*/;
   tListIter< tTriangle > triIter( triList );
   tMeshListIter< tEdge > edgIter( edgeList );
   tMeshListIter< tSubNode > nodIter( nodeList );
   tMeshListIter< tEdge > xedgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   tMeshList< tSubNode > tmpNodeList;
   tMeshListIter< tSubNode > tmpIter( tmpNodeList );
   tArray< double > p0, p1, p2, xy, xyz, xy1, xy2;
   tSubNode *cn;
   tPtrList< tTriangle > triptrList;
   tPtrListNode< tTriangle > *tpListNode;
   tPtrListIter< tTriangle > tpIter( triptrList );

   //check for triangles with edges which intersect (an)other edge(s)
   //newedg = new tEdge;
   while( flipped )
   {
      flipped = FALSE;

      // Make a list of triangles containing at least one moving vertex
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      {
         for( i=0; i<3; i++ )
         {
            cn = static_cast<tSubNode *>(ct->pPtr(i));
            if( cn->Meanders() ) break;
         }
         if( i!=3 ) triptrList.insertAtBack( ct );
      }
        //for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      for( ct = tpIter.FirstP(); !(triptrList.isEmpty());
           triptrList.removeFromFront( ct ), ct = tpIter.FirstP() )
      {
           //cout<<"PA: check triangle "<<ct->id<<", w edges "
           //<<ct->e[0]->id<<", "<<ct->e[1]->id<<", "<<ct->e[2]->id<<endl<<flush;
         if( !NewTriCCW( ct ) )
         {
            flipped = TRUE;
            for( i=0, j=0; i<3; i++ )
            {
               if( ct->pPtr(i)->getBoundaryFlag() != kNonBoundary ) j++;
            }
            if( j > 1 )
            {
               for( i=0, j=0; i<3; i++ )
               {
                  subnodePtr = static_cast<tSubNode *>(ct->pPtr(i));
                  subnodePtr->RevertToOldCoords();
               }
            }
            else
            {   
               crossed = FALSE;
               for( i=0; i<3; i++ )
               {
                  cn = static_cast<tSubNode *>(ct->pPtr(i));
                  if( cn->Meanders() )
                  {
                     cedg = ct->ePtr( (i+2)%3 );
                     spokIter.Reset( cn->getSpokeListNC() );
                     for( ce = spokIter.FirstP(); !( spokIter.AtEnd() );
                          ce = spokIter.NextP() )
                     {
                        if( Intersect( ce, cedg ) )
                        {
                           if( ct->tPtr(i) == 0 ) //boundary has been crossed
                           {
                              subnodePtr = static_cast<tSubNode *>(ct->pPtr(i));
                              subnodePtr->RevertToOldCoords();
                           }
                           else
                           {
                              crossed = TRUE;
                              ctop = ct->tPtr(i);
                              xy = cn->getNew2DCoords();
                                //check to make sure the opposite tri is still CCW;
                                //if so, check whether the point has moved into it;
                                //otherwise delete node and re-add it
                              if( NewTriCCW( ctop ) && InNewTri( xy, ctop ) )
                              {
                                   //if node has simply moved into 'opposite' triangle;
                                   //remove opposite tri from ptr list, flipedge,
                                   //add two new tri's to ptr list.
                                 for( rmtri = tpIter.FirstP();
                                      tpIter.ReportNextP() != ctop && !(tpIter.AtEnd());
                                      rmtri = tpIter.NextP() );
                                 if( !(tpIter.AtEnd()) ) //ctop is in tri ptrlist
                                 {
                                    tpListNode = tpIter.NodePtr();
                                    triptrList.removeNext( rmtri, tpListNode );
                                 }                           
                                 nv = ct->nVOp( ctop );
                                 nvopp = ctop->nVOp( ct );
                                 //cout << "call FlipEdge from CTEI for edge between nodes "
                                 //     << ct->pPtr( (nv+1)%3 )->getID() << " and "
                                 //     << ct->pPtr( (nv+2)%3 )->getID() << endl;
                                 FlipEdge( ct, ctop, nv, nvopp );
                                 rmtri = triIter.LastP();
                                 triptrList.insertAtBack( rmtri );
                                 rmtri = triIter.PrevP();
                                 triptrList.insertAtBack( rmtri );
                              }
                              else
                                    //things have gotten complicated and it's probably
                                    //easier to just delete the node and add it again
                                    //at the new location
                              {
                                 if( LocateTriangle( xy[0], xy[1] ) != 0 )
                                 {
                              
                                      //tempNode = *cn;
                                      //find spoke tri's in tri ptr list and remove them
                                    for( ce = spokIter.FirstP(); !(spokIter.AtEnd());
                                         ce = spokIter.NextP() )
                                    {
                                       rmtri = TriWithEdgePtr( ce );
                                       for( /*tri =*/ tpIter.FirstP();
                                            tpIter.ReportNextP() != rmtri &&
                                                !(tpIter.AtEnd());
                                                      /*tri =*/ tpIter.NextP() );
                                       if( !(tpIter.AtEnd()) ) //rmtri is in tri ptrlist
                                       {
                                          tpListNode = tpIter.NodePtr();
                                          triptrList.removeNext( rmtri, tpListNode );
                                       }
                                    }
                                      //delete the node;
                                    xyz = cn->getNew3DCoords();
                                    //cout << "delete node at " << xyz[0] << ", " << xyz[1]
                                    //     << ", " << xyz[2] << endl << flush;
                                    tmpNodeList.insertAtBack( *cn );
                                    DeleteNode( cn, kRepairMesh );
                                 }
                                 else
                                 {
                                    subnodePtr = static_cast<tSubNode *>(ct->pPtr(i));
                                    subnodePtr->RevertToOldCoords();
                                 }
                              }
                           }
                           break;
                        }
                     }
                  }
                  if( crossed ) break;
               }
            }      
         }
      }
   }
   
   // Update coordinates of moving nodes. TODO: make general
   for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
       if ( cn->Meanders() ) cn->UpdateCoords();//Nic, here is where x&y change
   for( cn = tmpIter.FirstP(); !(tmpIter.AtEnd()); cn = tmpIter.NextP() )
   {
      if ( cn->Meanders() ) cn->UpdateCoords();//Nic, here is where x&y change
      //cout << "add node at " << cn->getX() << ", " << cn->getY() << ", "
      //     << cn->getZ() << endl << flush;
      cn->getSpokeListNC().Flush();
      cn = AddNode( *cn );
      assert( cn!=0 );
   }
   
/*   for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
   {
      if( ct->tPtr(0) != 0 ) id0 = ct->tPtr(0)->getID();
      else id0 = -1;
      if( ct->tPtr(1) != 0 ) id1 = ct->tPtr(1)->getID();
      else id1 = -1;
      if( ct->tPtr(2) != 0 ) id2 = ct->tPtr(2)->getID();
      else id2 = -1;
      cout << "end of CTEI tri " << ct->getID() << " with nbrs "
           << id0 << ", " << id1 << ", and " << id2 << endl;
   }*/
   //cout << "finished, " << nnodes << endl;

}//end CheckTriEdgeIntersect()

               
/*****************************************************************************\
**
**  tMesh::MoveNodes (formerly PreApply)
**
**  Once the new coordinates for moving nodes have been established, this
**  function is called to update the node coordinates, modify the 
**  triangulation as needed, and update the mesh geometry (Voronoi areas,
**  edge lengths, etc) through a series of calls to helper functions.
**  
**  Interpolation is performed on nodes with layering (3D vertical 
**  component) here. TODO: make interpolation general, perhaps by
**  defining a virtual tNode function called "AlertNodeMoving" or
**  some such.
**
**      Inputs: time -- simulation time (for layer updating)
**      Data members updated: Mesh elements & their geometry
**      Called by:  called outside of tMesh by routines that compute
**                  node movement (e.g., stream meandering, as implemented
**                  by tStreamMeander)
**      Calls: CheckTriEdgeIntersect, CheckLocallyDelaunay, UpdateMesh,
**             LocateTriangle, tLNode::LayerInterpolation
**      Created: SL
**      Modifications:
**       - added interpFlag parameter to make layer interpolation optional,
**         so it needn't be called for tectonic motions (GT 4/00)     
**
\*****************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MoveNodes( double time, int interpFlag )
{
   cout << "MoveNodes()... time " << time <<flush << endl;
   tSubNode * cn;  
   tMeshListIter< tSubNode > nodIter( nodeList );

   //Before any edges and triangles are changed, layer interpolation
   //must be performed.
   if( layerflag && time > 0.0 ) {   
      tTriangle *tri;
      tArray<double> newxy(2);
      for(cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP()){
         newxy=cn->getNew2DCoords();
         if( (cn->getX()!=newxy[0]) || (cn->getY()!=newxy[1]) ){
            //Nic - there may be some issues if boundary nodes make up
            //the triangle.
            //cout<<"a point will be moved in movenodes"<<endl;
            tri = LocateTriangle( newxy[0], newxy[1] );
            if( interpFlag )
                cn->LayerInterpolation( tri, newxy[0], newxy[1], time );
            // TODO: is there a way to make this general, e.g. virtual fn?
         }
      }
   }

   //check for triangles with edges which intersect (an)other edge(s)
   CheckTriEdgeIntersect(); //calls tLNode::UpdateCoords() for each node
   //resolve any remaining problems after points moved
   CheckLocallyDelaunay();
   UpdateMesh();
   CheckMeshConsistency();  // TODO: remove this debugging call for release
   //cout << "MoveNodes() finished" << endl;
}


/*****************************************************************************\
**
**  tMesh::AddNodesAround
**
**  Densifies the mesh in the vicinity of a given node (centerNode) by
**  adding new nodes at the coordinates of the centerNode's Voronoi
**  vertices.
**
**  Properties of each node are initially those of the centerNode, except
**  z which is computed using interpolation by getVoronoiVertexXYZList.
**
**      Inputs: centerNode -- the node around which to add new nodes
**              time -- simulation time (for layer updating)
**      Data members updated: Mesh elements & their geometry
**      Called by:  called outside of tMesh by routines that handle
**                  adaptive meshing
**      Calls: AddNode, UpdateMesh, tNode::getVoronoiVertexXYZList
**      Created: GT, for dynamic mesh updating, Feb 2000  
**
\*****************************************************************************/
template<class tSubNode>
void tMesh< tSubNode >::
AddNodesAround( tSubNode * centerNode, double time )
{
   tList< Point3D > vvtxlist;  // List of V. vertex (x,y,z) coords at a node
   tListIter< Point3D > vtxiter( vvtxlist );

   assert( centerNode!=0 );

   // Get a list of Voronoi vertex coords and add a new node at each
   // (note: we get the list first because the vertices will change
   // as soon as we add the first node)
   centerNode->getVoronoiVertexXYZList( &vvtxlist );
   tLNode tmpnode = *centerNode;  // New node to be added -- passed to AddNode
   Point3D *xyz;  // Coordinates of current vertex

   // Here we add a new node at each vertex. Note that the call to
   // getVoronoiVertexListXYZList will compute a z value at each vertex
   // using plane (linear) interpolation.
   for( xyz=vtxiter.FirstP(); !(vtxiter.AtEnd()); xyz=vtxiter.NextP() )
   {
      //cout << "COORDS: x " << xyz->x << " y " << xyz->y << " z " << xyz->z << endl;
      tmpnode.set3DCoords( xyz->x, xyz->y, xyz->z );  // Assign to tmpnode
      //cn->TellAll();
      //cout << "Before addition\n";
      AddNode( tmpnode, FALSE, time );  // Add the node
   }
   UpdateMesh();

}



#ifndef NDEBUG
/*****************************************************************************\
**
**      DumpEdges(), DumpSpokes(), DumpTriangles(), DumpNodes(): debugging
**         routines which simply write out information pertaining to the mesh;
**      DumpNodes() calls DumpSpokes for each node;
**      DumpSpokes() takes a pointer to a node as an argument.
**
**      Created: SL 1/98
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
DumpEdges()
{
   tMeshListIter< tEdge > edgIter( edgeList );
   tEdge *ce;
   tTriangle *ct;
   int tid;
   cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      ct = TriWithEdgePtr( ce );
      tid = ( ct != 0 ) ? ct->getID() : -1;
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << "; in tri "
           << tid << " (flw " << ce->getBoundaryFlag() << ")" << endl;
   }
}


template<class tSubNode>
void tMesh<tSubNode>::
DumpSpokes( tSubNode *cn )
{
   tEdge *ce;
   tPtrListIter< tEdge > spokIter( cn->getSpokeListNC() );
   cout << "node " << cn->getID() << " with spoke edges " << endl;
   for( ce = spokIter.FirstP(); !( spokIter.AtEnd() ); ce = spokIter.NextP() )
   {
      cout << "   " << ce->getID()
          << " from node " << ce->getOriginPtrNC()->getID()
              << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }
}


template<class tSubNode>
void tMesh<tSubNode>::
DumpTriangles()
{
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct, *nt;
   int tid0, tid1, tid2;
   cout << "triangles:" << endl;
   for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
   {
      nt = ct->tPtr(0);
      tid0 = ( nt != 0 ) ? nt->getID() : -1;
      nt = ct->tPtr(1);
      tid1 = ( nt != 0 ) ? nt->getID() : -1;
      nt = ct->tPtr(2);
      tid2 = ( nt != 0 ) ? nt->getID() : -1;
      cout << ct->getID() << " with vertex nodes "
           << ct->pPtr(0)->getID() << ", "
           << ct->pPtr(1)->getID() << ", and "
           << ct->pPtr(2)->getID() << "; edges "
           << ct->ePtr(0)->getID() << ", "
           << ct->ePtr(1)->getID() << ", and "
           << ct->ePtr(2)->getID() << "; nbr triangles "
           << tid0 << ", "
           << tid1 << ", and "
           << tid2 << endl;
   }
}

template<class tSubNode>
void tMesh<tSubNode>::
DumpNodes()
{
   tMeshListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn;
   cout << "nodes: " << endl;
   for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
   {
      cout << " at " << cn->getX() << ", " << cn->getY() << ", " << cn->getZ()
           << "; bndy: " << cn->getBoundaryFlag() << "; ";
      DumpSpokes( cn );
   }
}
#endif
