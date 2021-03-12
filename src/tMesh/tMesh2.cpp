/***************************************************************************/
/**
**  @file
**  @brief interface to Tipper triangulator
*/
/***************************************************************************/

#include <stdlib.h>
#include "TipperTriangulator.h"

template< class tSubNode >
void tMesh< tSubNode >::
RemovePointDuplicates( tSubNode &tempnode, tList<double> &xList, 
		       tList<double> &yList, tList<double> &zList )
{ // use global functions used with Tipper triangulator
  // to eliminate duplicate points:
  double x(0.0);
  double y(0.0);
  double z(0.0);
  int npoints = xList.getSize();
  point *p = new point[npoints];
  vector<double> zVec( npoints );
  for( int i=0; i<npoints; ++i )
    {
      xList.removeFromFront( x );
      yList.removeFromFront( y );
      zList.removeFromFront( z );
      p[i] = point( x, y, i );
      zVec[i] = z;
      //      if( zVec[i] == 0.0 ) std::cout << "zero elev\n";
    }
  //sort the points and move duplicated points to the end:
  int npoints_unique;
  tt_sort_only( npoints, p, npoints_unique );
  // add only unique points to nodeList:
  miNextNodeID = 0;
  for( int i=0; i<npoints_unique; ++i )
    {
      tempnode.set3DCoords( p[i].x(), p[i].y(), zVec[p[i].id()] );
      tempnode.setID( miNextNodeID++ );
//       if( tempnode.getZ() == 0.0 )
// 	std::cout << "zero elev after sort\n";
      nodeList.insertAtActiveBack( tempnode );
    }
}      


/**************************************************************************\
**
**   tMesh::MakeMeshFromScratchTipper( infile )
**
**   as tMesh::MakeMeshFromScratch but uses Tipper's triangulation
**   algorithm.
**
**   Created: 07/2002, Arnaud Desitter
**   Calls:
**   Parameters:
**   Modified:
**
\**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromScratchTipper( const tInputFile &infile, tRand &rand )
{
   //std::cout << "In MGFS, calling node constr w/ infile\n";

   // Parameters defined in Input File
   ParamMMFS_t Param(infile);

   // Generate points
   {
     tPtrList< tSubNode > bndList;
     MakePointBoundary(Param, infile, bndList, rand);
     MakePointInterior(Param, infile, false, rand);
     nnodes = nodeList.getSize();
   }

   // call triangulator based on Tipper's method
   BuildDelaunayMeshTipper();

   std::cout<<"MakeMeshFromScratchTipper done.\n";
}

/**************************************************************************
**
**   tMesh::MakeMeshFromPointsTipper( infile )
**
**   Similar to tMesh::MakeMeshFromPoints but uses Tipper's triangulation
**   algorithm.
**
**   Created: 07/2002, Arnaud Desitter, Greg Tucker, Oxford
**   Modified: 08/2002, MIT
**
**************************************************************************/

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromPointsTipper( const tInputFile &infile ){
  {
    int numpts;                      // no. of points in mesh
    char pointFilenm[80];            // name of file containing (x,y,z,b) data
    std::ifstream pointfile;         // the file (stream) itself

    //Read Points
    infile.ReadItem( pointFilenm, sizeof(pointFilenm), "POINTFILENAME" );
    pointfile.open( pointFilenm );
    if( !pointfile.good() ){
      std::cerr << "\nPoint file name: '" << pointFilenm << std::endl;
      ReportFatalError( "I can't find a file by this name." );
    }

    std::cout<<"\nReading in '"<<pointFilenm<<"' points file..."<<std::endl;
    pointfile >> numpts;
    if( !pointfile.good() ) {
      std::cerr << "\nPoint file name: '" << pointFilenm << "'\n";
      ReportFatalError( "Error reading points file: check the file format." );
    }
    // temporary node used to create node list (creation is costly)
    const tSubNode aNode( infile );
    //Read point file, make Nodelist
    for( int i=0; i<numpts; i++ ){
      double x, y, z;
      int bnd;
      if( pointfile.eof() ) {
	std::cout << "\nReached end-of-file while reading points.\n" ;
	ReportFatalError("Invalid point file.");
      }
      pointfile >> x >> y >> z >> bnd;
      if( pointfile.fail() ) {
	std::cerr << "\nPoint file name: '" << pointFilenm
		  << "' - point " << i << std::endl;
	ReportFatalError( "I can't read the point above." );
      }

      tSubNode tempnode( aNode );
      tempnode.set3DCoords( x, y, z);
      tempnode.setBoundaryFlag( IntToBound(bnd) );
      miNextNodeID = i;
      tempnode.setID( miNextNodeID );
      if( bnd<0 || bnd>3 ){
	ReportFatalError("Invalid boundary code.");
      }

      switch(IntToBound(bnd)){
      case kNonBoundary:
	nodeList.insertAtActiveBack( tempnode );
	break;
      case kOpenBoundary:
	nodeList.insertAtBoundFront( tempnode );
	break;
      case kClosedBoundary:
	nodeList.insertAtBack( tempnode );
	break;
      }
    }
    pointfile.close();
  }
  nnodes = nodeList.getSize();

  // call triangulator based on Tipper's method
  BuildDelaunayMeshTipper();

  std::cout<<"MakeMeshFromPointsTipper done.\n";
}

/**************************************************************************\
**
**   tMesh::BuildDelaunayMeshTipper()
**
**   once nodeList is set up, build a Delaunay triangulation using Tipper's
**   algorithm.
**
**   Created: 07/2002, Arnaud Desitter
**   Calls:
**   Parameters:
**   Modified:
**
\**************************************************************************/
// edge numbering translation
static inline int e_t2c(int ei, bool o){ // Tipper to child
  return o? 2*ei : 2*ei+1;
}
static inline int e_t2c(const oriented_edge &oe){
  return e_t2c(oe.e(), oe.o());
}

template< class tSubNode >
void tMesh< tSubNode >::
BuildDelaunayMeshTipper( kUpdateMesh_t updateFlag /*= kUpdateMesh*/ )
{
   point *p = new point[nnodes];   // for Tipper triangulator
   {
     nodeListIter_t nodIter(nodeList);
     tSubNode* cn;
     int inode = 0;
     for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP()){
       p[inode] = point(cn->getX(), cn->getY(), cn->getID());
       ++inode;
     }

     if (0) { // DEBUG
       for (int i=0; i<nnodes; ++i){
	 std::cout << i << " x=" << p[i].x() << " y=" << p[i].y() << std::endl;
       }
     }
   }

   // call mesh generator based on Tipper's method
   std::cout << "Computing triangulation..." << std::flush;
   int nedgesl;
   int nelem;
   edge* edges(0);
   elem* elems(0);
   int nnodes_unique;
   tt_sort_triangulate(nnodes, p, &nnodes_unique, &nedgesl, &edges, &nelem, &elems);
   if (nnodes != nnodes_unique)
     {
      std::cerr << "\nDuplicated points: '" << std::endl;
      ReportFatalError( "Please fix your points file." );
     }

   std::cout << "done.\n";

   if (0) { // DEBUG
     std::cout << "After sort_triangulate:\n";
     for( int iedge=0; iedge < nedgesl; ++iedge) {
       std::cout << "edge " << 2*iedge
		 << " from=" << p[edges[iedge].from].id()
		 << " (" << p[edges[iedge].from].x() << ","
		 << p[edges[iedge].from].y() << ")"
		 << " to=" << p[edges[iedge].to].id()
		 << " (" << p[edges[iedge].to].x() << ","
		 << p[edges[iedge].to].y() << ")"
		 << std::endl;
     }
   }

   // set sizes
   this->nedges = 2*nedgesl;
   this->ntri = nelem;
   {
     const tIdArrayNode_t NodeTable(nodeList); // for fast lookup per ID

     // Create and initialize the edge list by creating two temporary edges
     // (which are complementary, ie share the same endpoints) and then
     // iteratively assigning values to the pair and inserting them onto the
     // back of the edgeList
     std::cout << "Creating edge list..." << std::flush;
     {
       for( int iedge = 0; iedge < nedgesl; ++iedge ) {
	 tSubNode *nodPtr1 = NodeTable[p[edges[iedge].from].id()];
	 tSubNode *nodPtr2 = NodeTable[p[edges[iedge].to].id()];
	 // Assign values: ID, origin and destination pointers,
	 // "flowallowed" status
	 tEdge
	   tempedge1( e_t2c(iedge,true) , nodPtr1, nodPtr2),
	   tempedge2( e_t2c(iedge,false), nodPtr2, nodPtr1);

	 // insert edge pair onto the list --- active
	 // part of list if flow is allowed, inactive if not

	 if (0) //DEBUG
	   std::cout << "Setting edges " << tempedge1.getID() << " and "
		     << tempedge2.getID() <<
	     (tempedge1.FlowAllowed() ? " as OPEN" : " as no-flux")
		     << std::endl;
	 if ( tempedge1.FlowAllowed() )
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
	     tEdge *e1 = edgeList.getLastNC()->getDataPtrNC();
	     edgeList.insertAtBack( tempedge2 );
	     tEdge *e2 = edgeList.getLastNC()->getDataPtrNC();
	     e1->setComplementEdge(e2);
	     e2->setComplementEdge(e1);
	   }
       }
     }
     const tIdArrayEdge_t EdgeTable(edgeList); // for fast lookup per ID
     std::cout << "done.\n";

     if (0) { //DEBUG
       std::cout << "JUST ADDED EDGES:\n";
       edgeListIter_t ei( edgeList );
       for( tEdge *ce=ei.FirstP(); !(ei.AtEnd()); ce=ei.NextP() ){
	 ce->TellCoords();
	 std::cout << tEdge::EdgeBoundName(ce->FlowAllowed()) << std::endl;
       }
     }

     // set up the lists of edges (spokes) connected to each node
     std::cout << "Setting up edg pointer..." << std::flush;
     {
       // connectivity point - sorted point
       tArray< int > p2sp(nnodes);
       for(int inodes=0;inodes!=nnodes;++inodes){
	 p2sp[p[inodes].id()] = inodes;
       }

       oriented_edge *oedge;
       tt_build_spoke(nnodes, nedgesl, edges, &oedge);

       nodeListIter_t nodIter(nodeList);
       tSubNode * cn;
       for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP()){
	 tEdge *edgPtr = EdgeTable[ e_t2c(oedge[p2sp[cn->getID()]]) ];
	 cn->setEdg( edgPtr );
       }
       delete [] oedge;
     }
     std::cout << "done.\n";

     // Assign ccwedg connectivity (that is, tell each edge about its neighbor
     // immediately counterclockwise). Likewise for cwedg connectivity.
     std::cout << "Setting up CCW and CW edges..." << std::flush;
     {
       for( int iedge=0; iedge<nedgesl; ++iedge){
	 {
	   const oriented_edge e1(iedge,true);
	   tEdge *curedg = EdgeTable[ e_t2c(e1) ];
	   const oriented_edge ccw_from = e1.ccw_edge_around_from(edges);
	   tEdge *ccwedg = EdgeTable[ e_t2c(ccw_from) ];
	   curedg->setCCWEdg( ccwedg );
	   ccwedg->setCWEdg( curedg );
	 }
	 {
	   const oriented_edge e2(iedge,false);
	   tEdge *curedg = EdgeTable[ e_t2c(e2) ];
	   const oriented_edge ccw_to = e2.ccw_edge_around_from(edges);
	   tEdge *ccwedg = EdgeTable[ e_t2c(ccw_to) ];
	   curedg->setCCWEdg( ccwedg );
	   ccwedg->setCWEdg( curedg );
	 }
       }
     }
     std::cout << "done.\n";

     std::cout << "Setting up triangle connectivity..." << std::flush;
     {
       int ielem;
       for ( ielem=0; ielem<nelem; ++ielem ) {
	 if (0) // DEBUG
	   std::cout << "TRI " << ielem << std::endl << std::flush;
	 if (0) { // DEBUG
	   std::cout << "p0=" << p[elems[ielem].p1].id() << " "
		<< "(" << p[elems[ielem].p1].x()
		<< "," << p[elems[ielem].p1].y() << "), "
		<< "p1=" << p[elems[ielem].p2].id() << " "
		<< "(" << p[elems[ielem].p2].x()
		<< "," << p[elems[ielem].p2].y() << "), "
		<< "p2=" << p[elems[ielem].p3].id() << " "
		<< "(" << p[elems[ielem].p3].x()
		<< "," << p[elems[ielem].p3].y() << "), "
		<< "e0=" << elems[ielem].e1 << " "
		<< "e1=" << elems[ielem].e2 << " "
		<< "e2=" << elems[ielem].e3 << std::endl;
	 }
	 tTriangle newtri( ielem,
			   NodeTable[p[elems[ielem].p1].id()],
			   NodeTable[p[elems[ielem].p2].id()],
			   NodeTable[p[elems[ielem].p3].id()],
			   EdgeTable[e_t2c(elems[ielem].e1,
					   elems[ielem].eo1)],
			   EdgeTable[e_t2c(elems[ielem].e2,
					   elems[ielem].eo2)],
			   EdgeTable[e_t2c(elems[ielem].e3,
					   elems[ielem].eo3)]
			   );
	 triList.insertAtBack( newtri );
       }
     }
   }
   // deallocation of some Tipper triangulator data structures
   delete [] edges; edges = 0;
   delete [] p;
   {
     const tIdArrayTri_t TriTable(triList); // for fast lookup per ID
     int ielem;
     for( ielem=0; ielem<nelem; ++ielem ) {
       tTriangle *ct = TriTable[ ielem ];
       ct->setTPtr( 0,
		    elems[ielem].t1>=0 ? TriTable[ elems[ielem].t1 ] : 0
		    );
       ct->setTPtr( 1,
		    elems[ielem].t2>=0 ? TriTable[ elems[ielem].t2 ] : 0
		    );
       ct->setTPtr( 2,
		    elems[ielem].t3>=0 ? TriTable[ elems[ielem].t3 ] : 0
		    );
     }
   }
   std::cout << "done.\n";

   // deallocation of remaining Tipper triangulator data structures
   delete [] elems;

   // set Maximum IDs
   SetmiNextEdgID( edgeList.getSize() );
   SetmiNextTriID( triList.getSize() );

   // assertions
   assert( edgeList.getSize() == 2*nedgesl );
   assert( triList.getSize() == nelem );
   if( updateFlag == kUpdateMesh )
     {
       CheckMeshConsistency();                     //remove CMC call from UM
       UpdateMesh(); //calls CheckMeshConsistency()  TODO: once bug-free,
     }
}


