/***************************************************************************/
/**
**  @file tMesh.cpp
**  @brief Functions for class tMesh (see tMesh.h) plus global
**         functions used by tMesh methods (formerly tGrid)
**
**  Summary of recent changes:
**    - modified LocateTriangle to implement triangle search
**      starting from a given location; modified constructors to set
**      initial value of mSearchOriginTriPtr, and modified ExtricateTri...
**      to avoid dangling ptr. GT, 1/2000
**    - added initial densification functionality, GT Sept 2000
**
**  $Id: tMesh.cpp,v 1.165 2003-06-19 15:03:02 childcvs Exp $
*/
/***************************************************************************/

#include "tMesh.h"

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
int Next3Delaunay( tPtrList< tSubNode > & /*nbrList*/,
                   tPtrListIter< tSubNode > &nbrIter )
{
   if (0) { //DEBUG
     static int ncalls = 0;
     ncalls++;
     cout << "Next3Delaunay? no calls =  " << ncalls << endl;
   }

   tSubNode *nbrnd = nbrIter.DatPtr();
   tPtrListIter< tSubNode > nbrIterCopy( nbrIter );
   const tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   const tArray< double > p1( nbrIterCopy.NextP()->get2DCoords() );
   const tArray< double > p2( nbrIterCopy.NextP()->get2DCoords() );

   // If points aren't counter-clockwise, we know it's not Delaunay
   if( !PointsCCW( p0, p1, p2 ) ) return 0;

   // Otherwise, compare it to each of the other potential triangles
   // p0-p1-ptest (?) where ptest is one of the other points in the
   // ring
   tSubNode *cn;
   // Keep testing 'til we're back to p0
   while( ( cn  = nbrIterCopy.NextP() ) != nbrnd )
   {
     const tArray< double > ptest = cn->get2DCoords();
     if( !TriPasses( ptest, p0, p1, p2 ) )
       return 0;
   }
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
  miNextNodeID(0),
  miNextEdgID(0),
  miNextTriID(0),
  seed(0),
  layerflag(false),
  mSearchOriginTriPtr(0)
{
  if (0) //DEBUG
   cout<<"tMesh()"<<endl;
}

//copy constructor (created 11/99, GT)
//WARNING: this constructor relies on the behavior of assignment operations
//in tMeshList, tList, and tPtrList, which may or not give the desired
//results! Caveat emptor! -GT
template< class tSubNode >
tMesh<tSubNode>::tMesh( tMesh const *originalMesh )
  :
  nnodes(originalMesh->nnodes),
  nedges(originalMesh->nedges),
  ntri(originalMesh->ntri),
  nodeList(originalMesh->nodeList),
  edgeList(originalMesh->edgeList),
  triList(originalMesh->triList),
  miNextNodeID(originalMesh->miNextNodeID),
  miNextEdgID(originalMesh->miNextEdgID),
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
  layerflag(false),
  mSearchOriginTriPtr(0)
{
   // mSearchOriginTriPtr:
   // initially set search origin (tTriangle*) to zero:
   // in initialisation list

   // As "layerflag" is used in this constructor, we compute it now.
   {
     int help;
     help = infile.ReadItem( help, "OPTINTERPLAYER" );
     layerflag = BOOL(help>0);
   }
   // option for reading/generating initial mesh
   int read;
   read = infile.ReadItem( read, "OPTREADINPUT" );
   switch (read){
   case 0:
     MakeMeshFromScratch( infile ); //create new mesh with parameters
     break;
   case 10:
     MakeMeshFromScratchTipper( infile ); //create new mesh with parameters
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
   case 12:
     MakeMeshFromPointsTipper( infile ); //create new mesh from list of points
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
	  << " 10 -- idem, using Tipper triangulator\n"
	  << "  1 -- read mesh from input data files\n"
	  << "  2 -- create mesh from a list of (x,y,z,b) points\n"
	  << " 12 -- idem, using Tipper triangulator\n"
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
~tMesh() {
  mSearchOriginTriPtr = 0;
  if (1)//DEBUG
    cout << "    ~tMesh()" << endl;
}


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

   if (0) //DEBUG
     cout<<"in MakeLayersFromInputData..."<<endl;

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

   assert(sumbr>0.999 && sumbr<1.001);

   layhelp.setCtime(0.);
   layhelp.setRtime(0.);
   //layhelp.setFlag(0);
   layhelp.setErody(0.);
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
**   tMesh::MeshDensification
**
**   The user may wish to densify the starting mesh uniformly by adding
**   a new node at the circumcenter of every triangle. That option is
**   implemented here. We simply iterate through the list of triangles
**   and add a node at the circumcenter of each. In doing so we take
**   advantage of the fact that the circumcenter is also a Voronoi vertex,
**   and is pointed to by each of the clockwise-directed edges in the
**   triangle. The z value of each new node is obtained by linear (plane)
**   interpolation from the 3 triangle vertices.
**     "initMeshDensLevel" serves as both a flag indicating whether the
**   user wants densification, and as an indicator of the number of passes
**   (the "level") to make -- ie, the number of times we sweep through
**   adding a node in each triangle.
**     Added Sept. 2000, GT
**
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MeshDensification( tInputFile &infile )
{
   int initMeshDensLevel;
   initMeshDensLevel = infile.ReadItem( initMeshDensLevel, "OPTINITMESHDENS" );
   if( initMeshDensLevel)
   {
     tSubNode tempnode( infile );
     int j;  // Level counter
     int nnewpoints;  // No. of new points added in a given pass
     tArray<double> newx, newy, newz;   // Lists of new coords
     tArray<double> zvals(3);   // z values of a triangle's 3 nodes
     tempnode.setBoundaryFlag( 0 );  // assumed all interior points
     tListIter< tTriangle > triIter( triList );
     tTriangle * ct;
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
       {
	 for( int i=0; i<ntri; i++ )    // loop through the triangles
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
       }

       // Now loop through and add the nodes
       {
	 for( int i=0; i<nnewpoints; i++ )
	   {
	     tempnode.set3DCoords( newx[i], newy[i], newz[i] );  // assign them
	     miNextNodeID = nnodes+i;
	     tempnode.setID( miNextNodeID );
	     AddNode( tempnode );        // Add the new node
	   }
       }
     }  // end of current densification level
   } // end of optional mesh densification
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
   if (0) //DEBUG
     cout << "nnodes, nedges, ntri: "
	  << nnodes << " " << nedges << " " << ntri << endl;
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
      miNextNodeID = i;
      tempnode.setID( miNextNodeID );
      bound = input.boundflag[i];
      assert( bound >= 0 && bound <= 2 );
      tempnode.setBoundaryFlag( bound );
      if( bound == kNonBoundary )
          nodeList.insertAtActiveBack( tempnode );
      else if( bound == kOpenBoundary )
          nodeList.insertAtBoundFront( tempnode );
      else
          nodeList.insertAtBack( tempnode );       //kClosedBoundary
      if (0) { // DEBUG
	cout << input.x[i] << input.y[i] << input.z[i]
	     << input.boundflag[i] << endl;
	cout << tempnode.getBoundaryFlag() << " ";
	cout << nodeList.getLast()->getDataPtr()->getBoundaryFlag() << endl;
      }
   }
   cout << "done.\n";

   {
     const tIdArray< tSubNode > NodeTable(nodeList); // for fast lookup per ID
     // Create and initialize the edge list by creating two temporary edges
     // (which are complementary, ie share the same endpoints) and then
     // iteratively assigning values to the pair and inserting them onto the
     // back of the edgeList
     cout << "Creating edge list..." << flush;
     {
       for( miNextEdgID = 0; miNextEdgID < nedges-1; miNextEdgID+=2 )
	 {
	   if (0) //DEBUG
	     cout << input.orgid[miNextEdgID] << " "
		  << input.destid[miNextEdgID] << endl;
	   tSubNode *nodPtr1 = NodeTable[ input.orgid[miNextEdgID] ];
	   tSubNode *nodPtr2 = NodeTable[ input.destid[miNextEdgID] ];
	   if (0) //DEBUG
	     cout << nodPtr1->getID() << "->" << nodPtr2->getID() << endl;

	   // Assign values: ID, origin and destination pointers,
	   // "flowallowed" status
	   tEdge
	     tempedge1( miNextEdgID  , nodPtr1, nodPtr2),
	     tempedge2( miNextEdgID+1, nodPtr2, nodPtr1);

	   // insert edge pair onto the list --- active
	   // part of list if flow is allowed, inactive if not
	   if (0) //DEBUG
	     cout << "Setting edges " << tempedge1.getID() << " and "
		  << tempedge2.getID() <<
	       (tempedge1.FlowAllowed() ? " as OPEN" : " as no-flux")
		  << endl;
	   if( tempedge1.FlowAllowed() )
	     {
	       edgeList.insertAtActiveBack( tempedge1 );
	       tEdge *e1 = edgeList.getLastActive()->getDataPtrNC();
	       edgeList.insertAtActiveBack( tempedge2 );
	       tEdge *e2 = edgeList.getLastActive()->getDataPtrNC();
	       e1->setComplementEdge(e2);
	       e2->setComplementEdge(e1);
	     }
	   else
	     {
	       edgeList.insertAtBack( tempedge1 );
	       tEdge *e1 = edgeList.getLast()->getDataPtrNC();
	       edgeList.insertAtBack( tempedge2 );
	       tEdge *e2 = edgeList.getLast()->getDataPtrNC();
	       e1->setComplementEdge(e2);
	       e2->setComplementEdge(e1);
	     }
	 }
     }
     const tIdArray< tEdge > EdgeTable(edgeList); // for fast lookup per ID
     cout << "done.\n";

     //DEBUG
     if (0) {
       cout << "JUST ADDED EDGES:\n";
       tMeshListIter< tEdge > ei( edgeList );
       for( tEdge *ce=ei.FirstP(); !(ei.AtEnd()); ce=ei.NextP() )
	 {
	   ce->TellCoords();
	   cout << ce->FlowAllowed() << endl;
	 }
     }

     // set up the lists of edges (spokes) connected to each node
     // (GT added code to also assign the 1st edge to "edg" as an alternative
     // to spokelist implementation)
     cout << "Setting edg pointers and spoke configs..." << flush;
     // (no real lists: "insert" functions simply set appropriate ccwedg's and
     // cwedg's)
     {
       tMeshListIter< tSubNode > nodIter( nodeList );
       int ret = nodIter.First();
       assert( ret != 0 );
       tSubNode * curnode;
       do
	 {
	   curnode = nodIter.DatPtr();
	   tSpkIter spkI( curnode );
	   const int edgid1 = input.edgid[curnode->getID()];  //fix of above error
	   spkI.insertAtBack(  EdgeTable[edgid1] );
	   int ne;
	   for( ne = input.nextid[edgid1]; ne != edgid1; ne = input.nextid[ne] )
	     {
	       if( ne>=nedges )
		 {
		   cerr << "Warning: edge " << edgid1
			<< " has non-existant ccw edge " << ne << '\n';
		   cerr << "This is likely to be a problem in the edge "
		     "input file" << endl;
		 }
	       spkI.insertAtBack( EdgeTable[ne] );
	     }
	 }
       while( (curnode = nodIter.NextP()) != 0 );
     }
     cout << "done.\n";

     cout << "Setting up triangle connectivity..." << flush;
     {
       for ( int i=0; i<ntri; i++ )
	 {
	   if (0) // DEBUG
	     cout << "TRI " << i << endl;
	   tTriangle newtri( i,
			     NodeTable[ input.p0[i] ],
			     NodeTable[ input.p1[i] ],
			     NodeTable[ input.p2[i] ],
			     EdgeTable[ input.e0[i] ],
			     EdgeTable[ input.e1[i] ],
			     EdgeTable[ input.e2[i] ]
			     );
	   triList.insertAtBack( newtri );
	 }
     }
   }
   {
     const tIdArray< tTriangle > TriTable(triList); // for fast lookup per ID
     for( int i=0; i<ntri; ++i )
       {
	 tTriangle *ct = TriTable[ i ];
	 ct->setTPtr( 0,
		      ( input.t0[i]>=0 ) ? TriTable[ input.t0[i] ] : 0
		      );
	 ct->setTPtr( 1,
		      ( input.t1[i]>=0 ) ? TriTable[ input.t1[i] ] : 0
		      );
	 ct->setTPtr( 2,
		      ( input.t2[i]>=0 ) ? TriTable[ input.t2[i] ] : 0
		      );
       }
   }
   cout << "done.\n";

   MeshDensification( infile ); // optional mesh densification

   if (0) { // DEBUG
     tMeshListIter< tEdge > ei( edgeList );
     tEdge * ce;
     cout << "JUST BEFORE UPDATEMESH\n";
     for( ce=ei.FirstP(); !(ei.AtEnd()); ce=ei.NextP() )
       {
	 ce->TellCoords();
	 cout << ce->getVEdgLen() << " " << ce->getBoundaryFlag() << endl;
       }
   }

   UpdateMesh();
   CheckMeshConsistency();

   cout << "end of tMesh( input )" << endl;
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
   while( (ce = tmpbndList.removeFromFront()) != 0 ) // current bndy edge removed from list
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
            be = tmpbndList.removeFromBack();
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
            be = tmpbndList.removeFromBack();
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
               cn = tmpnodList.removeFromBack();
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
      tempnode.set3DCoords( 0., 0., 0. );
      tempnode.setID( miNextNodeID );
      n = ROUND( Param.xGrid / Param.delGrid );
      tempnode.setBoundaryFlag( kOpenBoundary );
      nodeList.insertAtBack( tempnode );
      bndList.insertAtBack( nodIter.LastP() );
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=1, miNextNodeID++; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0., 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0., dist, 0. );
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
         tempnode.set3DCoords( dist, 0., 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>=0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0., dist, 0. );
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
         tempnode.set3DCoords( dist, 0., 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<=n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0. );
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
         tempnode.set3DCoords( 0., dist, 0. );
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
         tempnode.set3DCoords( dist, 0., 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=0; i<n; i++, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( Param.xGrid, dist, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.xGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, Param.yGrid, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = ROUND( Param.yGrid / Param.delGrid );
      for( i=n; i>0; i--, miNextNodeID++ )
      {
         dist = i * Param.delGrid + 0.0001 * Param.delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0., dist, 0. );
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
         tempnode.set3DCoords( dist, 0., 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );

         // If user wants outlet on this side between this and the next pt,
         // create the outlet now
         if( Param.yout == 0 && Param.xout > dist && Param.xout < dist + Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0. );
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
         tempnode.set3DCoords( Param.xGrid, dist, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( Param.xout == Param.xGrid && Param.yout > dist && Param.yout < dist + Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0. );
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
         tempnode.set3DCoords( dist, Param.yGrid, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( Param.yout == Param.yGrid && Param.xout < dist && Param.xout > dist - Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0. );
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
         tempnode.set3DCoords( 0., dist, 0. );
         tempnode.setID( miNextNodeID );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( Param.xout == 0 && Param.yout < dist && Param.yout > dist - Param.delGrid )
         {
            tempnode.set3DCoords( Param.xout, Param.yout, 0. );
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
MakePointInterior( const ParamMMFS_t &Param, tInputFile &infile,
		   bool makeMesh )
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
            xyz[2] = Param.mElev + Param.randElev * ( ran3( &seed ) - 0.5 );
            if( Param.boundType == kOppositeSidesOpen || Param.kSloped)
            {
               slope = Param.upperZ / Param.yGrid;
               xyz[2] += slope * xyz[1] - Param.mElev;
            }
            tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
            tempnode.setID( miNextNodeID );
	    if (makeMesh)
	      AddNode( tempnode );
	    else
	      nodeList.insertAtActiveBack(tempnode);
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
         xyz[2] = Param.mElev + Param.randElev * ( ran3( &seed ) - 0.5 );
         if( xyz[0] != 0 && xyz[0] != Param.xGrid && xyz[1] != 0 && xyz[1] != Param.yGrid )
         {
            tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
            tempnode.setID( miNextNodeID );
	    if (makeMesh)
	      AddNode( tempnode );
	    else
	      nodeList.insertAtActiveBack(tempnode);
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
      tempnode.set3DCoords( Param.xout, Param.yout, 0. );
      tempnode.setID( miNextNodeID );
      miNextNodeID++;
      if (makeMesh)
	AddNode( tempnode );
      else
	nodeList.insertAtActiveBack(tempnode);
   }
}

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromScratch( tInputFile &infile )
{
   if (0) //DEBUG
     cout << "In MGFS, calling node constr w/ infile\n";

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
   MakePointInterior(Param, infile, true);

   // Now finalize the initialization by updating mesh properties
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

   cout<<"MakeMeshFromPoints"<<endl;

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
   if( !pointfile.good() ){
     cerr << "\nPoint file name: '" << pointFilenm << "'\n";;
     ReportFatalError( "I can't find a file by this name." );
   }
   x.setSize( numpts );
   y.setSize( numpts );
   z.setSize( numpts );
   bnd.setSize( numpts );
   for( i=0; i<numpts; i++ )
   {
      if( pointfile.eof() )
          ReportFatalError( "Reached end-of-file while reading points." );
      pointfile >> x[i] >> y[i] >> z[i] >> bnd[i];
      if( pointfile.fail() ) {
	cerr << "\nPoint file name: '" << pointFilenm
	     << "' - point " << i << endl;
	ReportFatalError( "I can't read the point above." );
      }
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

   cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: " << nedges << " NT: " << ntri << endl;
   cout << "c4\n";

   // Now add the points one by one to construct the mesh.
   for( i=0; i<numpts; i++ )
   {
      //cout << "IN MGFP c0, ADDING NODE " << i << endl;
      //Xtempnode.setID( miNextNodeID );
      tempnode.set3DCoords( x[i],y[i],z[i] );
      tempnode.setBoundaryFlag( bnd[i] );
      //if(bnd[i]==kNonBoundary && z[i]<0)
      //  cout<<"problem at x "<<x[i]<<" y "<<y[i]<<endl;
      AddNode( tempnode );
   }

   cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: " << nedges << " NT: " << ntri << endl;

   // We no longer need the supertriangle, so remove it by deleting its
   // vertices.
   DeleteNode( stp1, kNoRepair );
   DeleteNode( stp2, kNoRepair );
   DeleteNode( stp3, kNoRepair );

   cout << "3 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: " << nedges << " NT: " << ntri << endl;

   MeshDensification( infile ); // optional mesh densification

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
         cn = AddNode( tempnode, /*updatemesh =*/ false );
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
      tSpkIter sI( cn );
      for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
          if( ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary )
              ++numActNbrs;
      if( numActNbrs ) cn->setZ( minz );
      else deletelist.insertAtBack( cn );
   }
   for( cn = dI.FirstP(); !( dI.AtEnd() ); cn = dI.FirstP() )
   {
      DeleteNode( cn, kNoRepair );
      cn  = deletelist.removeFromFront();
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
      cn = AddNode( tempnode, /*updatemesh =*/ false );
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
      tSpkIter sI( cn );
      for( ce = sI.FirstP(); !( sI.AtEnd() ); ce = sI.NextP() )
          if( ce->getDestinationPtr()->getBoundaryFlag() == kNonBoundary )
              ++numActNbrs;
      if( numActNbrs ) cn->setZ( minz );
      else deletelist.insertAtBack( cn );
   }
   for( cn = dI.FirstP(); !( dI.AtEnd() ); cn = dI.FirstP() )
   {
      DeleteNode( cn, kNoRepair );
      cn = deletelist.removeFromFront();
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
CheckMeshConsistency( bool boundaryCheckFlag /* default: true */)
{
   tMeshListIter<tSubNode> nodIter( nodeList );
   tMeshListIter<tEdge> edgIter( edgeList );
   tListIter<tTriangle> triIter( triList );
   tNode * cn, * org, * dest;
   tEdge * ce, * cne, * ccwedg;
   tTriangle * ct, * optr;
   bool boundary_check_ok;
   const bool verbose = false;
   int i, nvop;

   // Edges: make sure complementary pairs are together in the list
   // (each pair Ei and Ei+1, for i=0,2,4,...nedges-1, should have the same
   // endpoints but the opposite orientation)
   for( ce=edgIter.FirstP(); !(edgIter.AtEnd()); ce=edgIter.NextP() )
   {
      assert( ce != NULL);
      cne = edgIter.NextP();
      assert( cne != NULL);
      if( ce->getOriginPtrNC() != cne->getDestinationPtrNC()
          || ce->getDestinationPtrNC() != cne->getOriginPtrNC() )
      {
          cerr << "EDGE #" << ce->getID()
               << " must be followed by its complement in the list\n";
          goto error;
      }
      if( ce->getComplementEdge() != cne
          || ce != cne->getComplementEdge() )
	{
          cerr << "EDGE #" << ce->getID() << " EDGE #" << cne->getID()
               << " do not point to each other\n";
          goto error;
	}
      if( ce->getID()%2 != 0 || cne->getID()%2 != 1 )
	{
          cerr << "EDGE #" << ce->getID()
	       << " and EDGE #" << cne->getID()
	       << " should have respectively an even and odd ID.\n";
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
      tEdge *cwedg;
      if( !(cwedg=ce->getCWEdg() ) )
      {
         cerr << "EDGE #" << ce->getID()
              << " does not point to a valid clockwise edge\n";
         goto error;
      }
      if( cwedg->getOriginPtrNC()!=org )
      {
         cerr << "EDGE #" << ce->getID()
              << " points to a CW edge with a different origin\n";
         goto error;
      }
      if( cwedg->getDestinationPtrNC()==dest )
      {
         cerr << "EDGE #" << ce->getID()
              << " points to a CCW edge with the same destination\n";
         goto error;
      }
      if ( ce->getCCWEdg()->getCWEdg() != ce ||
	   ce->getCWEdg()->getCCWEdg() != ce)
      {
	cerr << "EDGE #" << ce->getID()
	     << ": inconsistency related to pointers to CW and CCW edges\n";
	goto error;
      }

      if( org==dest )
      {
         cerr << "EDGE #" << ce->getID()
              << " has the same origin and destination nodes\n";
         goto error;
      }

   }
   // Edges: check active/boundary list
   if (CheckMeshListConsistency(edgeList, "edge"))
     goto error;
   if (verbose)
     cout << "EDGES PASSED\n";

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
                            boundaryCheckFlag ) ? false : true;
      i = 0;
      // Loop around the spokes until we're back at the beginning
      do
      {

         if( ce->getDestinationPtrNC()->getBoundaryFlag()!=kClosedBoundary )
             boundary_check_ok = true;  // OK, there's at least one open nbr
         i++;
         if( i>kMaxSpokes ) // Uh-oh, an infinite loop
         {
            cerr << "NODE #" << cn->getID()
                 << " has more than " << kMaxSpokes << " spokes.\n";
	    cerr << "This error can result from a very high differential "
		 << "mesh resolution.\n"
		 << "Check input parameters controlling mesh densification.\n";
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
      tSpkIter sIter( cn );
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
   // Nodes: check active/boundary list
   if (CheckMeshListConsistency(nodeList, "node"))
     goto error;
   if (verbose)
     cout << "NODES PASSED\n";

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
         if( (optr = ct->tPtr(i)) != 0 )
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
	 // check flip test
	 if( ct->tPtr(i) != 0 )
	   {
	     if ( CheckForFlip( ct, i, false, false ) ) {
               cerr << "TRIANGLE #" << ct->getID()
                    << ": flip test failed for edge opposite to vertex "
                    << cn->getID() << ".\n";
               goto error;
	     }
	   }
      }
   }
   if (verbose)
     cout << "TRIANGLES PASSED\n";
   if (verbose)
     cout << "MESH PASSED\n";
   return;

  error:
   ReportFatalError( "Error in mesh consistency." );

}
#undef kMaxSpokes

template< class tSubNode >
void tMesh< tSubNode >::
Print() const
{
   triList.print();
   nodeList.print();
   edgeList.print();
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
   if (0) //DEBUG
     cout << "setVoronoiVertices()..." << endl;
   tListIter< tTriangle > triIter( triList );
   tTriangle * ct;

   // Find the Voronoi vertex associated with each Delaunay triangle
   for( ct = triIter.FirstP(); !(triIter.AtEnd()); ct = triIter.NextP() )
   {
      const tArray< double > xy = ct->FindCircumcenter();
      //cout << "setVoronoiVertices(): " << xy[0] << " " << xy[1];
      // Assign the Voronoi point as the left-hand point of the three edges
      // associated with the current triangle
      ct->ePtr(0)->setRVtx( xy );
      ct->ePtr(1)->setRVtx( xy );
      ct->ePtr(2)->setRVtx( xy );

      if (0) {//DEBUG
	cout << "FOR edges: ";
	for( int i=0; i<=2; i++ )
          cout << ct->ePtr(i)->getID() << " ("
               << ct->ePtr(i)->getOriginPtr()->getID() << ","
               << ct->ePtr(i)->getDestinationPtr()->getID() << ") ";
	cout << ", v verts are:\n";
	const tArray< double > xy_ = ct->ePtr(0)->getRVtx();
	cout << "  setVoronoiVertices(): " << xy_[0] << " " << xy_[1] << endl;
      }
   }
   if (0) //DEBUG
     cout << "setVoronoiVertices() finished" << endl;
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
   if (0) //DEBUG
     cout << "CalcVAreas()..." << endl;
   tSubNode* curnode;
   tMeshListIter< tSubNode > nodIter( nodeList );

   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
      curnode->ComputeVoronoiArea();

   }
   //cout << "CalcVAreas() finished" << endl;
}


/**************************************************************************\
**
**  tMesh::DeleteNode( tSubNode *, int =1 )
**    (see DeleteNode( tListNode<tSubNode> *, int =1 ) below)
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
DeleteNode( tSubNode const *node, bool repairFlag )
{
  tMeshListIter< tSubNode > nodIter( nodeList );
  if( nodIter.Get( node->getID() ) )
    return DeleteNode( nodIter.NodePtr(), repairFlag );
  return 0;
}

/**************************************************************************\
**
**  tMesh::DeleteNode( tListNode<tSubNode> *, int =1 )
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
**  5/2003 SL, AD
**
\**************************************************************************/
template< class tSubNode >
int tMesh< tSubNode >::
DeleteNode( tListNode< tSubNode >* nodPtr, bool repairFlag)
{
   tPtrList< tSubNode > nbrList;
   tSubNode *node = nodPtr->getDataPtrNC();

   if (0) //DEBUG
     cout << "DeleteNode: " << node->getID() << " at " << node->getX() << " "
	  << node->getY() << " " << node->getZ() << endl;
   //assert( repairFlag || node->getBoundaryFlag()==kClosedBoundary );

   // extricate node from mesh and get list of its neighbors:
   if( !( ExtricateNode( node, nbrList ) ) ) return 0;

   // remove node from nodeList:
   tSubNode nodeVal;
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

   if (0) { //DEBUG
     tPtrListIter< tSubNode > nbrIter( nbrList );
     cout << "leaving hole defined by \n"
	  << "   Node  x  y " << endl;
     for( nbrIter.First(); (!nbrIter.AtEnd()); nbrIter.Next() )
       {
	 cout << "   " << nbrIter.DatPtr()->getID() << "     "
	      << nbrIter.DatPtr()->getX() << "  "
	      << nbrIter.DatPtr()->getY() << endl;
       }
   }

   if( repairFlag )
   {
      nbrList.makeCircular();
      if( !RepairMesh( nbrList ) ) return 0;
   }

   //reset node id's
   ResetNodeID();

   if (0) { //DEBUG
     cout << "Mesh repaired" << endl;
     tMeshListIter< tSubNode > nodIter( nodeList );
     tSubNode *cn;
     for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
       {
	 cout << "node " << cn->getID() << endl;
       }
   }

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
   if (0) //DEBUG
     cout << "ExtricateNode: " << node->getID() << endl;
   tSpkIter spokIter( node );
   tEdge *ce;
   tSubNode *nbrPtr;

     //cout << "Removing spokes: " << endl;
     //assert( ExtricateEdge( edgptrIter.DatPtr() ) );
   for( ce = spokIter.FirstP(); !(spokIter.AtEnd()); ce = spokIter.FirstP() )
   {
      nbrPtr = static_cast< tSubNode * >(ce->getDestinationPtrNC());
      nbrList.insertAtBack( nbrPtr );
      // If node is a bdy make sure nbrs are also boundaries:
      if( node->getBoundaryFlag()
          && nbrPtr->getBoundaryFlag()==kNonBoundary )
      {
         nbrPtr->ConvertToClosedBoundary();
         nodeList.moveToBack( nbrPtr );
      }
      if( !DeleteEdge( ce ) )
	return 0;
   }

   //cout<<"nnodes decremented now "<<nnodes<<endl;
   if( tSpkIter(node).isEmpty() ) return 1;
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
   if (0) //DEBUG
     cout << "DeleteEdge(...) " << edgePtr->getID() << endl;
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

   //if( &edgeVal1 == 0 || &edgeVal2 == 0 ) return 0;
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
   if (0) //DEBUG
     cout << "ExtricateEdge: " << edgePtr->getID() << endl;
   //edgePtr->TellCoords();
   assert( edgePtr != 0 );
   //temporary objects:
   tEdge *ce, *cce;
   tMeshListIter< tEdge > edgIter( edgeList );
   tSpkIter spokIter;
   tListNode< tEdge > *listnodePtr1, *listnodePtr2;
   tTriangle* triPtr0, * triPtr1;

   //cout << "find edge in list; " << flush;
   ce = edgIter.GetP( edgePtr->getID() );  //NB: why necessary? isn't ce the
   // same as edgePtr??  Yes, why???
   // Puts edgIter at edgePtr's node!

   // WarnSpokeLeaving is virtual:
   ce->getOriginPtrNC()->WarnSpokeLeaving( ce );
   // Remove the edge from it's origin's spokelist:
   spokIter.Reset( ce->getOriginPtrNC() );
   if( spokIter.Get( ce ) )
       if( spokIter.Next() )
	 spokIter.removePrev();

   // Find the triangle that points to the edge
   //cout << "find triangle; " << flush;
   triPtr0 = TriWithEdgePtr( edgePtr );
   // Find the edge's complement
   listnodePtr1 = edgIter.NodePtr();
   assert( listnodePtr1 != 0 );

   //cout << "find complement; " << flush;
   if( edgePtr->getID()%2 == 0 ) {
     cce = edgIter.NextP();
   } else {
     cce = edgIter.PrevP();
   }
   listnodePtr2 = edgIter.NodePtr();

   // just to make sure
   assert( listnodePtr1->getDataPtrNC() == ce );
   assert( listnodePtr2->getDataPtrNC() == cce );
   assert( ce->getComplementEdge() == cce );
   assert( cce->getComplementEdge() == ce );

   //Need to make sure that "edg" member of node was not pointing
   //to one of the edges that will be removed.  Also, may be implications
   //for some types of subnodes, so take care of that also.
   // WarnSpokeLeaving is virtual:
   cce->getOriginPtrNC()->WarnSpokeLeaving( cce );

   // Find the triangle that points to the edges complement
   //cout << "find other triangle; " << flush;
   triPtr1 = TriWithEdgePtr( cce );

   // set edges' tri ptrs to zero:
   ce->setTri( 0 );
   cce->setTri( 0 );

   //if triangles exist, delete them
   //cout << "conditionally calling deletetri from extricateedge\n";
   if( triPtr0 != 0 )
       if( !DeleteTriangle( triPtr0 ) ) return 0;
   if( triPtr1 != 0 )
       if( !DeleteTriangle( triPtr1 ) ) return 0;
   //update complement's origin's spokelist
   spokIter.Reset( cce->getOriginPtrNC() );
   if( spokIter.Get( cce ) )
       if( spokIter.Next() )
	 spokIter.removePrev();

   //Since WarnSpokeLeaving can set a node to a boundary node if
   //There is no longer a legit place to flow, we need to check
   //to see if nodece and nodecce are now boundaries, and take proper action.

   if( ce->getBoundaryFlag() )
   {
      //move edges to back of list
      edgeList.moveToBack( listnodePtr1 );
      edgeList.moveToBack( listnodePtr2 );
   }
   else
   {
      //move edges to front of list
      edgeList.moveToFront( listnodePtr2 );
      edgeList.moveToFront( listnodePtr1 );
   }
   nedges-=2;
   return 1;
}

//*************************************************************************
// ClearEdge: function to remove edge and its complement from the
//   spokes of the origin nodes, call WarnSpokeLeaving for the origins,
//   and call ClearTriangle for the triangles associated with the edges.
//   Leaves the edges still pointing to the nodes and triangles and
//   to each other.
// 3/99 SL
// 4/2003 AD
//*************************************************************************
template< class tSubNode >
int tMesh< tSubNode >::
ClearEdge( tEdge* ce ) const
{
   tEdge *cce;
   tSpkIter spokIter;
   tTriangle* triPtr0, * triPtr1;
   //Need to make sure that "edg" member of node was not pointing
   //to one of the edges that will be removed.  Also, may be implications
   //for some types of subnodes, so take care of that also.
   // WarnSpokeLeaving is virtual:
   ce->getOriginPtrNC()->WarnSpokeLeaving( ce );
   // Remove the edge from it's origin's spokelist:
   spokIter.Reset( ce->getOriginPtrNC() );
   if( spokIter.Get( ce ) )
       if( spokIter.Next() )
	 spokIter.removePrev();

   // Find the triangle that points to the edge
   triPtr0 = TriWithEdgePtr( ce );

   // Find the edge's complement
   cce = ce->getComplementEdge();

   //Need to make sure that "edg" member of node was not pointing
   //to one of the edges that will be removed.  Also, may be implications
   //for some types of subnodes, so take care of that also.
   // WarnSpokeLeaving is virtual:
   cce->getOriginPtrNC()->WarnSpokeLeaving( cce );

   // Find the triangle that points to the edges complement
   triPtr1 = TriWithEdgePtr( cce );
   //if triangles exist, delete them
   if( triPtr0 != 0 )
       if( !ClearTriangle( triPtr0 ) )
       {
          cerr << "from ClearEdge(...), ClearTriangle(...) failed \n";
          return 0;
       }

   if( triPtr1 != 0 )
       if( !ClearTriangle( triPtr1 ) )
       {
          cerr << "from ClearEdge(...), ClearTriangle(...) failed \n";
          return 0;
       }

     //update complement's origin's spokelist
   spokIter.Reset( cce->getOriginPtrNC() );
   if( spokIter.Get( cce ) )
       if( spokIter.Next() )
           spokIter.removePrev();
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
LocateTriangle( double x, double y, bool useFuturePosn)
{
   if (0) //DEBUG
     cout << "\nLocateTriangle (" << x << "," << y << ")\n";
   tListIter< tTriangle > triIter( triList );  //lt
   tTriangle *lt = ( mSearchOriginTriPtr != 0 ) ? mSearchOriginTriPtr
       : triIter.FirstP();
   int online = -1;

   // it starts from the first triangle,
   // searches through the triangles until the point is on
   // the same side of all the edges of a triangle.
   // "lt" is the current triangle and "lv" is the edge number.
   int n, lv=0;
   for (n=0 ;lv!=3 && lt; n++)
   {
      const tArray< double > xy1 = useFuturePosn ?
	lt->pPtr(lv)->FuturePosn():
	lt->pPtr(lv)->get2DCoords();
      const tArray< double > xy2 = useFuturePosn ?
	lt->pPtr( (lv+1)%3 )->FuturePosn():
	lt->pPtr( (lv+1)%3 )->get2DCoords();
      const double XY[] = {x, y};
      const double c =
	predicate.orient2d(xy1.getArrayPtr(), xy2.getArrayPtr(), XY);

      if ( c < 0.0 )
      {
         lt=lt->tPtr( (lv+2)%3 );
         lv=0;
         online = -1;
      }
      else
      {
         if( c == 0.0 ) online = lv;
         ++lv;
      }

      assert( n < 3*ntri );
   }
   if( online != -1 )
       if( lt->pPtr(online)->getBoundaryFlag() != kNonBoundary &&
           lt->pPtr( (online+1)%3 )->getBoundaryFlag() != kNonBoundary ) //point on bndy
           return 0;
   return lt;
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
   if (0) //DEBUG
     cout << "LocateNewTriangle" << endl;
   return
     LocateTriangle( x, y, true );
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
TriWithEdgePtr( tEdge *edgPtr ) const
{
   assert( edgPtr != 0 );
   return edgPtr->TriWithEdgePtr();
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
DeleteTriangle( tTriangle const * triPtr )
{
   if (0) //DEBUG
     cout << "DeleteTriangle(...) " << triPtr->getID() << endl;
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

//*************************************************************************
// ClearTriangle: function to "un-point" neighbor triangles. Leaves all
//   the triangle's own pointers unaffected. Does not un-point the edges
//   pointing to the triangle because this is called in the process of
//   flipping an edge and, to reinitialize the triangles, we'll need the
//   edges to point to triangles in order to set the triangle's nbr ptrs.
// 3/99 SL
// 4/2003 AD
//*************************************************************************
template< class tSubNode >
int tMesh< tSubNode >::
ClearTriangle( tTriangle const *triPtr ) const
{
  int i, j;
  for( i=0; i<3; i++ )
    for( j=0; j<3; j++ )
      if( triPtr->tPtr(i) != 0 )
	if( triPtr->tPtr(i)->tPtr(j) == triPtr )
	  triPtr->tPtr(i)->setTPtr( j, 0 );
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
ExtricateTriangle( tTriangle const *triPtr )
{
   if (0) //DEBUG
     cout << "ExtricateTriangle" << endl;
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct;

   // Find the triangle on the list
   for( ct = triIter.FirstP(); ct != triPtr && !( triIter.AtEnd() );
        ct = triIter.NextP() );
   if( ( triIter.AtEnd() ) ) return 0;

   // Tell your neighboring triangles that you're about to disappear by
   // setting their corresponding tPtr's to zero
   int i, j;
   for( i=0; i<3; i++ )
     for( j=0; j<3; j++ )
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
   if (0) //DEBUG
     cout << "RepairMesh: " << endl;
   if( nbrList.getSize() < 3 ) return 0;
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
         int ret = AddEdgeAndMakeTriangle( nbrList, nbrIter );
	 assert( ret );
           //remove "closed off" pt
	 /* tSubNode * meshnodePtr = */
         nbrList.removeNext( nbrIter.NodePtr() );
      }
        //else cout << "Not delaun\n";

      nbrIter.Next();                    //step forward once in nbrList
   }
   assert( nbrList.getSize() == 3 );
   assert( ntri == triList.getSize() );
   assert( nedges == edgeList.getSize() );
   assert( nnodes == nodeList.getSize() );       //make sure numbers are right
   int ret =    MakeTriangle( nbrList, nbrIter );             //make final triangle
   assert( ret );

   if (0) //DEBUG
     cout << "done" << endl;
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
**   - GT 2/99 -- added calls to WelcomeCCWNeighbor and AttachNewSpoke
**     to update CCW edge connectivity
**
\**************************************************************************/
//vertices of tri in ccw order; edges are added between node1 and node2
//TODO: comment/document this fn
template< class tSubNode >
int tMesh< tSubNode >::
AddEdge( tSubNode *node1, tSubNode *node2, tSubNode const *node3 )
{
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   if (0) //DEBUG
     cout << "AddEdge"
	  << "between nodes " << node1->getID()
	  << " and " << node2->getID() << " w/ ref to node "
	  << node3->getID() << endl;
   tEdge *ce, *nle, *le;

   {
     tEdge
       tempEdge1(miNextEdgID++, node1, node2),
       tempEdge2(miNextEdgID++, node2, node1);

     // Place new edge pair on the list: active back if not a boundary
     // edge, back otherwise
     if( tempEdge1.FlowAllowed() )
       {
	 edgeList.insertAtActiveBack( tempEdge1 );  //put edge1 active in list
	 tEdge *e1 = edgeList.getLastActive()->getDataPtrNC();
	 edgeList.insertAtActiveBack( tempEdge2 );  //put edge2 active in list
	 tEdge *e2 = edgeList.getLastActive()->getDataPtrNC();
	 e1->setComplementEdge(e2);
	 e2->setComplementEdge(e1);
	 nle = e1;
	 le = e2;
       }
     else
       {
	 edgeList.insertAtBack( tempEdge1 );        //put edge1 in list
	 tEdge *e1 = edgeList.getLast()->getDataPtrNC();
	 edgeList.insertAtBack( tempEdge2 );        //put edge2 in list
	 tEdge *e2 = edgeList.getLast()->getDataPtrNC();
	 e1->setComplementEdge(e2);
	 e2->setComplementEdge(e1);
	 nle = e1;
	 le = e2;
       }
   }

   //add pointers to the new edges to nodes' spokeLists:
   tSpkIter spokIter;
   spokIter.Reset( node2 );
   if( spokIter.isEmpty()
       // only one spoke
       || ( spokIter.ReportNextP() == spokIter.CurSpoke() )
       )
       spokIter.insertAtFront( le );
   else
   {
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      //make sure we found the right spoke; if not:
      if( spokIter.AtEnd() )
          for( ce = spokIter.FirstP();
               !( spokIter.AtEnd() );
               ce = spokIter.NextP() )
              if( PointsCCW( UnitVector( ce ),
                             UnitVector( le ),
                             UnitVector( spokIter.ReportNextP() ) ) )
                  break;
      //put edge2 in SPOKELIST:
      spokIter.insertAtNext( le );
   }
   spokIter.Reset( node1 );
   if( spokIter.isEmpty()
       // only one spoke
       || ( spokIter.ReportNextP() == spokIter.CurSpoke() )
       )
     spokIter.insertAtFront( nle );
   else
   {
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      //make sure we found the right spoke; if not:
      if( spokIter.AtEnd() )
          for( ce = spokIter.FirstP();
               !( spokIter.AtEnd() );
               ce = spokIter.NextP() )
              if( PointsCCW( UnitVector( ce ),
                             UnitVector( nle ),
                             UnitVector( spokIter.ReportNextP() ) ) )
              {
                 spokIter.Next();
                 break;
              }
      //put edge1 in SPOKELIST:
      spokIter.insertAtPrev( nle );
   }

   nedges+=2;

   // Reset edge id's
   ResetEdgeID();

   if (0) //DEBUG
     cout << "AddEdge() done\n" << flush;
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
AddEdgeAndMakeTriangle( tPtrList< tSubNode > & /*nbrList*/,
                        tPtrListIter< tSubNode > &nbrIter )
{
   if (0) //DEBUG
     cout << "AddEdgeAndMakeTriangle" << endl;
   tSubNode *cn, *cnn, *cnnn;
   cn = nbrIter.DatPtr();
   cnn = nbrIter.NextP();
   cnnn = nbrIter.ReportNextP();
   nbrIter.Prev();
   return AddEdgeAndMakeTriangle( cn, cnn, cnnn );
}

template< class tSubNode >
int tMesh< tSubNode >::
AddEdgeAndMakeTriangle( tSubNode* cn, tSubNode* cnn, tSubNode* cnnn )
{
   if( !AddEdge( cnnn, cn, cnn ) ) return 0;
   if( !MakeTriangle( cn, cnn, cnnn ) ) return 0;
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
MakeTriangle( tPtrList< tSubNode > const &nbrList,
              tPtrListIter< tSubNode > &nbrIter )
{
   assert( nbrList.getSize() == 3 );
   tSubNode *cn, *cnn, *cnnn;
   cn = nbrIter.FirstP();      // cn, cnn, and cnnn are the 3 nodes in the tri
   cnn = nbrIter.NextP();
   cnnn = nbrIter.NextP();
   nbrIter.Next();
   return MakeTriangle(cn, cnn, cnnn);
}

template< class tSubNode >
int tMesh< tSubNode >::
MakeTriangle( tSubNode *cn, tSubNode *cnn, tSubNode *cnnn )
{
   assert( cn != 0 && cnn != 0 && cnnn != 0 );
   assert( cn != cnn && cn != cnnn && cnn != cnnn );

   if (0) //DEBUG
     cout << "MakeTriangle" << endl;
   const tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
       p2( cnnn->get2DCoords() );

   // error check
   if( !PointsCCW( p0, p1, p2 ) )
   {
      cerr << "in MT nodes not CCW: " << cn->getID() << ", "
           << cnn->getID() << ", " << cnnn->getID()
	   << "; nor are new coords CCW " << endl;
      return 0;
   }

   // Create the new triangle and insert a pointer to it on the list.
   // Here, the triangle constructor takes care of setting pointers to
   // the 3 vertices and 3 edges. The neighboring triangle pointers are
   // initialized to zero.
   triList.insertAtBack( tTriangle( miNextTriID++, cn, cnn, cnnn ) );//put
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct;
   ct = triIter.LastP();            //ct now points to our new triangle
   assert( cn == ct->pPtr(0) );     //make sure we're where we think we are

   // To speed up future searches in calls to LocateTriangle, assign the
   // starting triangle, mSearchOriginTriPtr, to the our new triangle.
   // The idea here is that there's a good chance that the next point
   // to be added will be close to the current location. (added 1/2000)
   mSearchOriginTriPtr = ct;

   // Now we assign the neighbor triangle pointers. The loop successively
   // gets the spokelist for (p0,p1,p2) and sets cn to the next ccw point
   // (p1,p2,p0). It then finds the edge (spoke) that joins the two points
   // (p0->p1, p1->p2, p2->p0). These are the edges that are shared with
   // neighboring triangles (t2,t0,t1) and are pointed to by the neighboring
   // triangles. This means that in order to find neighboring triangle t2,
   // we need to find the triangle that points to edge (p0->p1), and so on.
   // In general, t((j+2)%3) is the triangle that points to edge
   // p(j)->p((j+1)%3).
   tTriangle *nbrtriPtr = 0;
   int j;
   for( j=0; j<3; j++ )
   {
      // Find edge ce that connects p(j)->p(j+1)
      tEdge *ce = ct->pPtr(j)->EdgToNod( ct->pPtr( (j+1)%3 ) );

      if( nbrtriPtr != 0 && TriWithEdgePtr( ce ) == nbrtriPtr )
      {
         if( PointsCCW( p0, p1, p2 ) )
             cerr << "something FUNNY going on";
         else cerr << "tri not CCW: " << nbrtriPtr->getID() << endl;
      }

      // Find the triangle, if any, that shares (points to) this edge
      // and assign it as the neighbor triangle t((j+2)%3).
      nbrtriPtr = TriWithEdgePtr( ce );

      ct->setTPtr( (j+2)%3, nbrtriPtr );      //set tri TRI ptr (j+2)%3

      // If a neighboring triangle was found, tell it that the current
      // new triangle is its neighbor too. We need to tell it which
      // neighbor we are (0, 1, or 2), and the mapping is like this:
      // if the nbr tri calls the shared edge (0,1,2) then we are its
      // nbr (1,2,0). (ie, tri_number = (edg_number+1)%3 )
      if( nbrtriPtr != 0 )
      {
     	int i;
	for( i=0; i<3; i++ )
	  {
            assert( nbrtriPtr->ePtr(i) != 0 );
            assert( ce != 0 );
            if( nbrtriPtr->ePtr(i) == ce ) break;
	  }
	assert( nbrtriPtr->ePtr(i) == ce );
	nbrtriPtr->setTPtr( (i+1)%3, ct );  //set NBR TRI ptr to tri
      }
   }
   ++ntri;

   //reset triangle id's (why needed??) because when we make a new item of any kind we
   //give it an id; how do we know what id to use (i.e., what's large enough but not
   //too large)? we find the id of the last item in the list and add one; if the items
   //in the list have been "mixed up", then we could assign an id already in use;
   //also, if for some reason numbers are systematically skipped, the id could blow up;
   //this step
   //may not be strictly necessary for triangles (it is for nodes and edges where
   //we have active and inactive members), but I'm sure it doesn't hurt; better safe
   //than sorry...
   ResetTriangleID();
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
AddNode( tSubNode &nodeRef, bool updatemesh, double time )
{
   const tArray< double > xyz( nodeRef.get3DCoords() );

   if (0) //DEBUG
     cout << "AddNode at " << xyz[0] << ", " << xyz[1]
	  << ", " << xyz[2] << " time "<<time<<endl;

   // Assign ID to the new node and insert it at the back of either the active
   // portion of the node list (if it's not a boundary) or the boundary
   // portion (if it is)
   nodeRef.setID( miNextNodeID );
   miNextNodeID++;

   if (0) //DEBUG
     cout << "call InsertNode" << endl;
   tSubNode* newNodePtr = InsertNode(&nodeRef, time);
   if(newNodePtr == 0)
     return 0;
   if (0) //DEBUG
     cout << "call CheckTrianglesAt" << endl;
   if( xyz.getSize() == 3 )
     CheckTrianglesAt( newNodePtr );

   //reset node id's
   ResetNodeID();
   newNodePtr->InitializeNode();

   if( updatemesh ) UpdateMesh();
   return newNodePtr;  // Return ptr to new node
}

/**************************************************************************\
** CheckTrianglesAt(tNode*)
**  SL 3/99
**  AD 4/2003
\**************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
CheckTrianglesAt( tSubNode* nPtr )
{
  if (0) //DEBUG
    cout << "CheckTrianglesAt()"<< endl;

  tPtrList< tTriangle > triptrList;
  // put all triangle surrounding the node nPtr in triptrlist
  {
    tSpkIter spokIter( nPtr );
    tEdge* ce;
    tTriangle *ct;
    for( ce = spokIter.FirstP(); !spokIter.AtEnd(); ce = spokIter.NextP() )
      if( (ct = ce->TriWithEdgePtr()) != 0 )
	triptrList.insertAtBack( ct );
  }

  // re-triangulate to fix mesh (make it Delaunay)
  MakeDelaunay( triptrList );
}

/*****************************************************************************\
**
**  tMesh::MakeDelaunay
**
**  Flip edges where necessary by calling iteratively CheckForFlip.
**  triPtrList contains a pointer list to triangles.
**
**  Moved from CheckTrianglesAt and CheckLocallyDelaunay - AD 5/2003
**
\*****************************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
MakeDelaunay( tPtrList< tTriangle > &triPtrList)
{
  tTriangle *at;
  int ctr = 0;
  while( (at = triPtrList.removeFromFront()) != 0 )
    {
      ++ctr;
      if( ctr > kLargeNumber ) // Make sure to prevent endless loops
	{
	  cerr << "MakeDelaunay(): flip checking forever" << endl;
	  ReportFatalError( "Bailing out." );
	}
      for( int i=0; i<3; i++ )
	{
	  // If a neighboring triangle exists across this face, check for flip
	  if( at->tPtr(i) != 0 )
	    {
	      tTriangle *tp = at->tPtr(i);
	      // Check triangle _at_ for a flip across face opposite vertex i,
	      // and do the flip if needed
	      if( CheckForFlip( at, i, true ) )
		{
		  // A flip occurred, insert the two new triangles
		  // which have been rebuilt with the same pointers
		  tPtrListIter< tTriangle > triPtrIter( triPtrList );
		  if( triPtrIter.Get( tp ) )
		    triPtrList.moveToBack( triPtrIter.NodePtr() );
		  else
		    triPtrList.insertAtBack( tp );
		  // at has been removed from triptrList, insert it directly
		  triPtrList.insertAtBack( at );
		  break;
		}
	    }
	}
    }
}

/**************************************************************************\
**  InsertNode: (Takes over some of the functionality of AddNode of older
**    versions of CHILD) Locates triangle in which new node falls, puts it
**    there and makes three new triangles (at end of list).
**    Does NOT check the mesh for edge flips to enforce Delaunay-ness!!!
**
**  SL 1/99
**  AD 4/2003
\**************************************************************************/
template< class tSubNode >
tSubNode * tMesh< tSubNode >::
InsertNode( tSubNode* newNodePtr, double time )
{
   if (0) //DEBUG
     cout << "tMesh::InsertNode()" << endl;
   tTriangle *tri = LocateTriangle( newNodePtr->getX(), newNodePtr->getY() );
   if( tri == 0 )
   {
      cerr << "tMesh::InsertNode(...): node coords out of bounds: "
	   << newNodePtr->getX() << " " << newNodePtr->getY() << endl;
      return 0;
   }
   if( layerflag && time > 0. )
       newNodePtr->LayerInterpolation( tri, newNodePtr->getX(), newNodePtr->getY(), time );

   // Insert and retrieve a pointer to the new node
   tSubNode *cn = AddToList(*newNodePtr);

   tSubNode *newNodePtr2 = AttachNode( cn, tri);

   return newNodePtr2;
}

/**************************************************************************\
** AddToList( tNode* newNodePtr ): adds new node to appropriate place in
**   nodeList.
** SL 3/99
** AD 4/2003
\**************************************************************************/
template< class tSubNode >
tSubNode * tMesh< tSubNode >::
AddToList( tSubNode const & newNode )
{
  // insert node at the back of either the
  // active portion of the node list (if it's not a boundary) or the
  // boundary portion (if it is)
  tMeshListIter< tSubNode > nodIter( nodeList );
  tSubNode *cn;
  if( newNode.getBoundaryFlag()==kNonBoundary ) {
    nodeList.insertAtActiveBack( newNode );
    cn = nodIter.LastActiveP();
  } else if( newNode.getBoundaryFlag() == kOpenBoundary ) {
    nodeList.insertAtBoundFront( newNode );
    cn = nodIter.FirstBoundaryP();
  } else {
    nodeList.insertAtBack( newNode );
    cn = nodIter.LastP();
  }
  assert( nodeList.getSize() == nnodes + 1 );
  ++nnodes;
  return cn;
}

/**************************************************************************\
** AttachNode( node, tri ): makes necessary edges and triangles to incor-
**   porate the node in the triangulation. Assumes node is already added
**   to nodeList. Triangle, tri, is deleted. Triangulation is not checked
**   for Delaunay-ness.
**
** SL 3/99
** AD 4/2003
\**************************************************************************/
template< class tSubNode >
tSubNode* tMesh< tSubNode >::
AttachNode( tSubNode* cn, tTriangle* tri )
{
  assert( tri != 0 && cn != 0 );
  int i;
  const tArray< double > xyz( cn->get3DCoords() );

  // flush its spoke list
  cn->setEdg( 0 );

  //make ptr list of triangle's vertices:
  tPtrList< tSubNode > bndyList;
  tSubNode *tmpPtr;
  for( i=0; i<3; i++ )
    {
      tmpPtr = static_cast<tSubNode *>(tri->pPtr(i));
      bndyList.insertAtBack( tmpPtr );
    }
  bndyList.makeCircular();

  //make 3 new triangles
  //cout << "creating new triangles\n" << flush;
  tPtrListIter< tSubNode > bndyIter( bndyList );
  tSubNode *node3 = bndyIter.FirstP();     // p0 in original triangle
  tSubNode *node2 = cn;                    // new node
  tSubNode *node1 = bndyIter.NextP();      // p1 in orig triangle
  tSubNode *node4 = bndyIter.NextP();      // p2 in orig triangle
  tArray< double > p1( node1->get2DCoords() ),
    p2( node2->get2DCoords() ), p3( node3->get2DCoords() ),
    p4( node4->get2DCoords() );
  int colinearedg = -1;


  if( xyz.getSize() == 3) // why would this ever not be the case?
                          // If we need to access new coords:
    //size of xyz is basically the flag; the 4th element is never used o.w.
    {
      //cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
      //     << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;
      if( !PointsCCW( p3, p1, p2 ) ||
	  !PointsCCW( p2, p1, p4 ) ||
	  !PointsCCW( p2, p4, p3 ) )
	{
	  cout << "new tri not CCW" << endl;
	  if( Orientation( p3, p1, p2 ) == 0 ) colinearedg = 1;
	  if( Orientation( p2, p1, p4 ) == 0 ) colinearedg = 2;
	  if( Orientation( p2, p4, p3 ) == 0 ) colinearedg = 0;
	}
    }
  else
    {
      // use virtual function that will return new coords for nodes
      p1 = node1->FuturePosn();
      p2 = node2->FuturePosn();
      p3 = node3->FuturePosn();
      p4 = node4->FuturePosn();
      if( !PointsCCW( p3, p1, p2 ) ||
	  !PointsCCW( p2, p1, p4 ) ||
	  !PointsCCW( p2, p4, p3 ) )
	cout << "new tri not CCW" << endl;
    }

  if( colinearedg < 0 ) {

     // Clear triangle in which the node falls:
     ClearTriangle( tri );

     // Here's how the following works. Let the old triangle vertices be A,B,C
     // and the new node N. The task is to create 3 new triangles ABN, NBC, and
     // NCA, and 3 new edge-pairs AN, BN, and CN.
     // First, edge pair BN is added. Then AEMT is called to create triangle
     // ABN and edge pair AN. AEMT is called again to create tri NBC and edge
     // pair CN. With all the edge pairs created, it remains only to call
     // MakeTriangle to create tri NCA.
     //cout << "calling AE, AEMT, AEMT, and MT\n" << flush;

     AddEdge( node1, node2, node3 );  //add edge between node1 and node2
     AddEdgeAndMakeTriangle( node3, node1, node2 ); // ABN
     AddEdgeAndMakeTriangle( node2, node1, node4 ); // NBC
     // use old triangle for last one instead of MakeTriangle:
     tri->InitializeTriangle( node2, node4, node3 );
  } else  {
     // need to make four new triangles:
     // first, delete the edge with which the new point is colinear
     // --this will delete the other, nbr, triangle:
     node1 = static_cast<tSubNode *>(tri->pPtr( (colinearedg+2)%3 ));
     node3 = static_cast<tSubNode *>(tri->pPtr( colinearedg ));
     node4 = static_cast<tSubNode *>(tri->pPtr( (colinearedg+1)%3 ));
     tTriangle* top = tri->tPtr( (colinearedg+1)%3 );
     tSubNode* node5 = static_cast<tSubNode *>(top->pPtr( top->nVOp( tri ) ));
     assert( node1 != 0 && node2 != 0 && node3 != 0 && node4 != 0
	     && node5 != 0 );
     i = DeleteEdge( tri->ePtr( colinearedg ) );
     assert( i > 0 );
     AddEdge( node1, node2, node4 );
     AddEdgeAndMakeTriangle( node2, node1, node5 );
     AddEdgeAndMakeTriangle( node2, node5, node3 );
     AddEdgeAndMakeTriangle( node2, node3, node4 );
     MakeTriangle( node2, node4, node1 );
  }

  return node2;
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
**  05/2003 AD: make it call AddToList, AttachNodem, etc.
**
\**************************************************************************/
//TODO: ; just assign coords
// to a dummy new node and call AddNode
template< class tSubNode >
tSubNode *tMesh< tSubNode >::
AddNodeAt( tArray< double > &xyz, double time )
{
   if (0) //DEBUG
     cout << "AddNodeAt " << xyz[0] << ", " << xyz[1] << ", "
	  << xyz[2] <<" time "<<time<< endl;
   if (0) //DEBUG
     cout << "locate tri" << endl;
   tTriangle *tri;
   if( xyz.getSize() == 3 ) tri = LocateTriangle( xyz[0], xyz[1] );
   else tri = LocateNewTriangle( xyz[0], xyz[1] );
   if( tri == 0 )
     return 0;

   tSubNode tempNode;
   tempNode.set3DCoords( xyz[0], xyz[1], xyz[2]  );
   if( layerflag && time > 0.0) tempNode.LayerInterpolation( tri, xyz[0], xyz[1], time );
   if( xyz.getSize() != 3 ) tempNode.setNew2DCoords( xyz[0], xyz[1] );
   tempNode.setBoundaryFlag( kNonBoundary );

   // Assign ID to the new node and insert it at the back of the active
   // portion of the node list.
   tempNode.setID( miNextNodeID );
   miNextNodeID++;
   tSubNode *cn = AddToList( tempNode );

   tSubNode *newNodePtr2 = AttachNode( cn, tri);
   if(newNodePtr2 == 0)
     return 0;
   //put 3 resulting triangles in ptr list
   if( xyz.getSize() == 3 )
   {
     CheckTrianglesAt( newNodePtr2 );
   }
   //reset node id's
   ResetNodeID();
   //nmg uncommented line below and added initialize line
   newNodePtr2->InitializeNode();

   UpdateMesh();

   if (0)//DEBUG
     cout << "AddNodeAt finished, " << nnodes << endl;
   return newNodePtr2;
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
getEdgeComplement( tEdge *edge ) const
{
   return edge->getComplementEdge();
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
   if (0) //DEBUG
     cout << "UpdateMesh()" << endl;

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
      assert( curedg != 0 ); // failure = complementary edges not consecutive
      curedg->setLength( len );
   } while( (curedg=elist.NextP()) );

   setVoronoiVertices();
   CalcVoronoiEdgeLengths();
   CalcVAreas();
   CheckMeshConsistency( false );  // debug only -- remove for release

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
bool tMesh< tSubNode >::
CheckForFlip( tTriangle * tri, int nv, bool flip, bool useFuturePosn )
{
   if( tri == 0 )  // TODO: is this just a bug check?
   {
      cout << "CheckForFlip: tri == 0" << endl;
      return false;
   }
   assert( nv < 3 );
   if (0) //DEBUG
     cout << "THIS IS CheckForFlip(...) " << tri->getID() << endl;
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
   tArray< double >
     p0( node0->get2DCoords() ),
     p1( node1->get2DCoords() ),
     p2( node2->get2DCoords() ),
     ptest( node3->get2DCoords() );

   // If "flip" flag isn't set and the node is a moving node, use "new"
   // coordinates rather than current coordinates
   if( !flip && useFuturePosn)
   {
     // use virtual function that will return new coords for nodes
     p0 = node0->FuturePosn();
     p1 = node1->FuturePosn();
     p2 = node2->FuturePosn();
     ptest = node3->FuturePosn();
   }

   // If p0-p1-p2 passes the test, no flip is necessary
   if( TriPasses( ptest, p0, p1, p2 ) ) return false;

   // Otherwise, a flip is needed, provided that the new triangles are
   // counter-clockwise (would this really ever happen??) and that the
   // node isn't a moving node (ie "flip" is true)
   if( flip )                     //and make sure there isn't already an edge?
   {
      if( !PointsCCW( p0, p1, ptest ) || !PointsCCW( p0, ptest, p2 ) )
          return false;
      //cout << "calling Flip edge from cff" << endl;
      FlipEdge( tri, triop, nv, nvop );
   }
     //cout << "finished" << endl;
   return true;
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
**                d
**               / \
**       tri->  /   \
**             /     \
**            a-------c
**             \     /
**              \   / <-triop
**               \ /
**                b
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
** Edges and triangles are re-used (SL)
** 5/2003 AD
\*******************************************************************/
template< class tSubNode >
void tMesh< tSubNode >::
FlipEdge( tTriangle * tri, tTriangle * triop ,int nv, int nvop )
{
   if (0) //DEBUG
     cout << "FlipEdge(...)..." << endl;
   tSubNode* na = static_cast<tSubNode *>(tri->pPtr(nv));
   tSubNode* nb = static_cast<tSubNode *>(tri->pPtr((nv+1)%3));
   tSubNode* nc = static_cast<tSubNode *>(triop->pPtr( nvop ));
   tSubNode* nd = static_cast<tSubNode *>(tri->pPtr((nv+2)%3));
   tEdge* edg = tri->ePtr( (nv+2)%3 );
   tEdge* edgop = triop->ePtr( (nvop+2)%3 );
   const bool move =
     BOOL(
	  tEdge::isFlowAllowed(na, nc) == !tEdge::isFlowAllowed(nb, nd)
	  );

   tListNode< tEdge >* enodePtr1=0;
   tListNode< tEdge >* enodePtr2=0;
   if( move )
   {
      if( edg->getID()%2 == 0 )
          enodePtr1 = edgeList.getListNode( edg );
      else
          enodePtr1 = edgeList.getListNode( edgop );
      enodePtr2 = enodePtr1->getNextNC();
      assert( enodePtr1 != 0 && enodePtr2 != 0 );
   }

   // does edg and edgop, tri and triop:
   ClearEdge( edg );
   // give edges' initialization routine nodes of tri in cw order:
   edg->InitializeEdge( nc, na, nd );
   edgop->InitializeEdge( na, nc, nb );

   if( move )
   {
      if( edg->FlowAllowed() )
      {
         edgeList.moveToActiveBack( enodePtr1 );
	 edgeList.setNActiveNodes(edgeList.getActiveSize()+1);
         edgeList.moveToActiveBack( enodePtr2 );
	 edgeList.setNActiveNodes(edgeList.getActiveSize()+1);
      }
      else
      {
         edgeList.moveToBack( enodePtr1 );
	 edgeList.setNActiveNodes(edgeList.getActiveSize()-1);
         edgeList.moveToBack( enodePtr2 );
	 edgeList.setNActiveNodes(edgeList.getActiveSize()-1);
      }
   }
   // give triangles' initialization routine nodes of tri in ccw order:
   tri->InitializeTriangle( nd, na, nc );
   triop->InitializeTriangle( nb, nc, na );
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
  if (0) //DEBUG
    cout << "CheckLocallyDelaunay()" << endl;
  tPtrList< tTriangle > triPtrList;

  // Search through tri list to find triangles with at least one
  // moving vertex, and put these on triPtrList
  //put each triangle into the stack
  {
    tListIter< tTriangle > triIter( triList );
    tSubNode *nodPtr;
    for( tTriangle *at = triIter.FirstP(); !( triIter.AtEnd() );
	 at = triIter.NextP() )
      {
	bool change = false;
	for( int i = 0; i < 3; i++ )
	  {
	    nodPtr = static_cast< tSubNode * >(at->pPtr(i));
	    if( nodPtr->isMobile() ) change = true;
	  }
	if( change ) triPtrList.insertAtBack( at );
      }
  }

  // re-triangulate to fix mesh (make it Delaunay)
  MakeDelaunay( triPtrList);
}

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
   if (0) //DEBUG
     cout << "CheckTriEdgeIntersect()..." << endl;
   //DumpNodes();
   int i, j, nv, nvopp;
   bool flipped = true;
   bool crossed;
   tSubNode *subnodePtr;
   tEdge *cedg, *ce;
   tTriangle *ct, *ctop, *rmtri;
   tListIter< tTriangle > triIter( triList );
   tMeshListIter< tSubNode > nodIter( nodeList );
   tMeshList< tSubNode > tmpNodeList;
   tMeshListIter< tSubNode > tmpIter( tmpNodeList );
   tSubNode *cn;
   tPtrList< tTriangle > triptrList;
   tPtrListNode< tTriangle > *tpListNode;
   tPtrListIter< tTriangle > tpIter( triptrList );

   //check for triangles with edges which intersect (an)other edge(s)
   //newedg = new tEdge;
   while( flipped )
   {
      flipped = false;

      // Make a list of triangles containing at least one moving vertex
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      {
         for( i=0; i<3; i++ )
         {
            cn = static_cast<tSubNode *>(ct->pPtr(i));
            if( cn->isMobile() ) break;
         }
         if( i!=3 ) triptrList.insertAtBack( ct );
      }
        //for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      for( ct = tpIter.FirstP(); !(triptrList.isEmpty());
           ct = triptrList.removeFromFront(), ct = tpIter.FirstP() )
      {
           //cout<<"PA: check triangle "<<ct->id<<", w edges "
           //<<ct->e[0]->id<<", "<<ct->e[1]->id<<", "<<ct->e[2]->id<<endl;
         if( !NewTriCCW( ct ) )
         {
            flipped = true;
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
               crossed = false;
               for( i=0; i<3; i++ )
               {
                  cn = static_cast<tSubNode *>(ct->pPtr(i));
                  if( cn->isMobile() )
                  {
                     cedg = ct->ePtr( (i+2)%3 );
		     tSpkIter spokIter( cn );
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
                              crossed = true;
                              ctop = ct->tPtr(i);
			      const tArray< double > xy = cn->getNew2DCoords();
                                //check to make sure the opposite tri is still CCW;
                                //if so, check whether the point has moved into it;
                                //otherwise delete node and re-add it
                              if( NewTriCCW( ctop ) && InNewTri( xy, ctop ) )
                              {
                                   //if node has simply moved into 'opposite' triangle;
                                   //remove opposite tri from ptr list, flipedge,
                                   //add two new tri's to ptr list.
				 for( tpIter.First();
                                      tpIter.ReportNextP() != ctop && !(tpIter.AtEnd());
                                      tpIter.Next() );
                                 if( !(tpIter.AtEnd()) ) //ctop is in tri ptrlist
                                 {
                                    tpListNode = tpIter.NodePtr();
                                    triptrList.removeNext( tpListNode );
                                 }
                                 nv = ct->nVOp( ctop );
                                 nvopp = ctop->nVOp( ct );
				 if (0) //DEBUG
				   cout << "call FlipEdge from CTEI for edge between nodes "
					<< ct->pPtr( (nv+1)%3 )->getID() << " and "
					<< ct->pPtr( (nv+2)%3 )->getID() << endl;
                                 FlipEdge( ct, ctop, nv, nvopp );
                                 triptrList.insertAtBack( ct );
                                 triptrList.insertAtBack( ctop );
                              }
                              else
                                    //things have gotten complicated and it's probably
                                    //easier to just delete the node and add it again
                                    //at the new location
                              {
                                 if( LocateTriangle( xy[0], xy[1] ) != 0 )
                                 {
                                      //find spoke tri's in tri ptr list and remove them
                                    for( ce = spokIter.FirstP(); !(spokIter.AtEnd());
                                         ce = spokIter.NextP() )
                                    {
                                       rmtri = TriWithEdgePtr( ce );
                                       for( tpIter.First();
					    tpIter.ReportNextP() != rmtri &&
					      !(tpIter.AtEnd());
					    tpIter.Next() );
                                       if( !(tpIter.AtEnd()) ) //rmtri is in tri ptrlist
                                       {
                                          tpListNode = tpIter.NodePtr();
                                          triptrList.removeNext( tpListNode );
                                       }
                                    }
				    //delete the node;
				    if (0) {//DEBUG
				      const tArray< double > xyz = cn->getNew3DCoords();
				      cout << "delete node at " << xyz[0] << ", " << xyz[1]
					   << ", " << xyz[2] << endl;
				    }
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

   // Update coordinates of moving nodes. (UpdateCoords is virtual)
   for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
     cn->UpdateCoords();//Nic, here is where x&y change
   for( cn = tmpIter.FirstP(); !(tmpIter.AtEnd()); cn = tmpIter.NextP() )
   {
     cn->UpdateCoords();//Nic, here is where x&y change
     //cout << "add node at " << cn->getX() << ", " << cn->getY() << ", "
     //     << cn->getZ() << endl;
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
MoveNodes( double time, bool interpFlag )
{
   if (1) //DEBUG
     cout << "MoveNodes()... time " << time << endl;

   //Before any edges and triangles are changed, layer interpolation
   //must be performed.
   if( layerflag && time > 0.0 ) {
     tSubNode * cn;
     tMeshListIter< tSubNode > nodIter( nodeList );
     for(cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP()){
       const tArray<double> newxy = cn->getNew2DCoords();
       if( (cn->getX()!=newxy[0]) || (cn->getY()!=newxy[1]) ){
	 //Nic - there may be some issues if boundary nodes make up
	 //the triangle.
	 //cout<<"a point will be moved in movenodes"<<endl;
	 tTriangle *tri;
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
   if (0) //DEBUG
     cout << "MoveNodes() finished" << endl;
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
      AddNode( tmpnode, false, time );  // Add the node
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
DumpSpokes( tSubNode *cn ) const
{
   tEdge *ce;
   tSpkIter spokIter( cn );
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

/*****************************************************************************\
**
**      ResetNodeID(): reset node ID in list order
**      Created: AD 5/2003 (refactored from code lying in various locations)
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
ResetNodeID()
{
  tSubNode *cn;
  tMeshListIter< tSubNode > nodIter( nodeList );
  for( cn = nodIter.FirstP(), miNextNodeID=0; !( nodIter.AtEnd() );
       cn = nodIter.NextP(), miNextNodeID++ )
    {
      cn->setID( miNextNodeID );
    }
}

/*****************************************************************************\
**
**      ResetEdgeID(): reset edge ID in list order
**      Created: AD 5/2003 (refactored from code lying in various locations)
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
ResetEdgeID()
{
  tEdge *ce;
  tMeshListIter< tEdge > edgIter( edgeList );
  for( ce = edgIter.FirstP(), miNextEdgID = 0; !( edgIter.AtEnd() );
       ce = edgIter.NextP(), miNextEdgID++ )
    {
      ce->setID( miNextEdgID );
    }
}

/*****************************************************************************\
**
**      ResetTriangleID(): reset triangle ID in list order
**      Created: AD 5/2003 (refactored from code lying in various locations)
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
ResetTriangleID()
{
  tTriangle *ct;
  tListIter< tTriangle > triIter( triList );
  for( ct = triIter.FirstP(), miNextTriID=0; !( triIter.AtEnd() );
       ct = triIter.NextP(), miNextTriID++ )
    {
      ct->setID( miNextTriID );
    }
}

/*****************************************************************************\
**
**      SetmiNextNodeID(): set miNextNodeID
**      Created: AD 5/2003
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
SetmiNextNodeID(int n_)
{
  miNextNodeID = n_;
}

/*****************************************************************************\
**
**      SetmiNextEdgID(): set miNextEdgID
**      Created: AD 5/2003
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
SetmiNextEdgID(int n_)
{
  miNextEdgID = n_;
}

/*****************************************************************************\
**
**      SetmiNextTriID(): set miNextTriID
**      Created: AD 5/2003
**
\*****************************************************************************/
template<class tSubNode>
void tMesh<tSubNode>::
SetmiNextTriID(int n_)
{
  miNextTriID = n_;
}

#include "tMesh2.cpp"

