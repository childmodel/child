#include <stdlib.h>
#include "TipperTriangulator.h"

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

// edge numbering translation
static inline int e_t2c(int ei, bool o){ // Tipper to child
  return o? 2*ei : 2*ei+1;
}
static inline int e_t2c(const oriented_edge &oe){
  return e_t2c(oe.e(), oe.o());
}

template< class tSubNode >
void tMesh< tSubNode >::
MakeMeshFromScratchTipper( tInputFile &infile )
{
   //cout << "In MGFS, calling node constr w/ infile\n";

   seed = infile.ReadItem( seed, "SEED" );

   // Parameters defined in Input File
   ParamMMFS_t Param(infile);

   // Generate points
   {
     tPtrList< tSubNode > bndList;
     MakePointBoundary(Param, infile, bndList);
     MakePointInterior(Param, infile, false);
     nnodes = nodeList.getSize();
   }
   point *p = new point[nnodes];   // for Tipper triangulator 
   {
     tMeshListIter< tSubNode > nodIter(nodeList);
     tSubNode* cn;
     int inode = 0;
     for( cn=nodIter.FirstP(); !(nodIter.AtEnd()); cn=nodIter.NextP()){
       p[inode].x = cn->getX();
       p[inode].y = cn->getY();
       p[inode].id = cn->getID();
       ++inode;
     }

     if (0) { // DEBUG
       for (int i=0; i<nnodes; ++i){
	 cout << i << " x=" << p[i].x << " y=" << p[i].y << endl;
       }
     }
   }
   const tIdArray< tSubNode > NodeTable(nodeList); // for fast lookup per ID

   // call mesh generator based on Tipper's method
   cout << "Computing triangulation..." << flush;
   int nedgesl;
   int nelem;
   edge* edges(0);
   elem* elems(0);
   tt_sort_triangulate(nnodes, p, &nedgesl, &edges, &nelem, &elems);
   cout << "done.\n";

   if (0) { // DEBUG
     cout << "After sort_triangulate:\n";
     for( int iedge=0; iedge < nedgesl; ++iedge) {
       cout << "edge " << 2*iedge
	    << " from=" << p[edges[iedge].from].id 
	    << " (" << p[edges[iedge].from].x << ","
	    << p[edges[iedge].from].y << ")"
	    << " to=" << p[edges[iedge].to].id
	    << " (" << p[edges[iedge].to].x << ","
	    << p[edges[iedge].to].y << ")"
	    << endl;
     }
   }

   // set sizes
   nedges = 2*nedgesl;
   ntri = nelem;

   // Create and initialize the edge list by creating two temporary edges
   // (which are complementary, ie share the same endpoints) and then
   // iteratively assigning values to the pair and inserting them onto the
   // back of the edgeList
   cout << "Creating edge list..." << flush;
   {
     for( int iedge = 0; iedge < nedgesl; ++iedge ) {
       tEdge tempedge1, tempedge2;
       int obnd, dbnd;
     
       // Assign values: ID, origin and destination pointers
       tempedge1.setID( e_t2c(iedge,true) );
       tempedge2.setID( e_t2c(iedge,false) );
       {
	 tSubNode *nodPtr1 = NodeTable[p[edges[iedge].from].id];
	 tempedge1.setOriginPtr( nodPtr1 );
	 tempedge2.setDestinationPtr( nodPtr1 );
	 obnd = (*nodPtr1).getBoundaryFlag();
       }
       {
	 tSubNode *nodPtr2 = NodeTable[p[edges[iedge].to].id];
	 tempedge1.setDestinationPtr( nodPtr2 );
	 tempedge2.setOriginPtr( nodPtr2 );
	 dbnd = (*nodPtr2).getBoundaryFlag();
       }
       
       // set the "flowallowed" status (FALSE if either endpoint is a
       // closed boundary, or both are open boundaries) 
       // and insert edge pair onto the list --- active
       // part of list if flow is allowed, inactive if not
       if( obnd == kClosedBoundary || dbnd == kClosedBoundary
	   || (obnd==kOpenBoundary && dbnd==kOpenBoundary) )
	 {
	   if (0) //DEBUG
	     cout << "setting edges " << tempedge1.getID() << " and "
		  << tempedge2.getID() << " as no-flux" << endl;
	   tempedge1.setFlowAllowed( 0 );
	   tempedge2.setFlowAllowed( 0 );
	   edgeList.insertAtBack( tempedge1 );
	   edgeList.insertAtBack( tempedge2 );
	 }
       else
	 {
	   if (0) //DEBUG
	     cout << "setting edges " << tempedge1.getID() << " and "
		  << tempedge2.getID() << " as OPEN" << endl;
	   tempedge1.setFlowAllowed( 1 );
	   tempedge2.setFlowAllowed( 1 );
	   edgeList.insertAtActiveBack( tempedge1 );
	   edgeList.insertAtActiveBack( tempedge2 );
	 }
     }
   }
   const tIdArray< tEdge > EdgeTable(edgeList); // for fast lookup per ID
   cout << "done.\n";

   if (0) { //DEBUG
     cout << "JUST ADDED EDGES:\n";
     tMeshListIter< tEdge > ei( edgeList );
     tEdge * ce;
     for( ce=ei.FirstP(); !(ei.AtEnd()); ce=ei.NextP() )
       {
	 ce->TellCoords();
	 cout << ce->FlowAllowed() << endl;
       }
   }

   // set up the lists of edges (spokes) connected to each node
   cout << "setting up spoke lists..." << flush;
   {
     // connectivity point - sorted point
     tArray< int > p2sp(nnodes);
     for(int inodes=0;inodes!=nnodes;++inodes){
       p2sp[p[inodes].id] = inodes;
     }

     tMeshListIter< tSubNode > nodIter(nodeList);
     oriented_edge *oedge;
     tt_build_spoke(nnodes, nedgesl, edges, &oedge);

     tSubNode * curnode;
     assert( nodIter.First() );
     do
       {
	 // first spoke
	 curnode = nodIter.DatPtr();
	 {
	   const int e1 = e_t2c(oedge[p2sp[curnode->getID()]]);
	   tEdge *edgPtr = EdgeTable[e1];
	   curnode->insertBackSpokeList( edgPtr );
	   curnode->setEdg( edgPtr );
	 }
	 // build rest of spoke list
	 const oriented_edge& oe_ref = oedge[p2sp[curnode->getID()]];
	 oriented_edge ccw_from = oe_ref.ccw_edge_around_from(edges);
	 while( ccw_from.e() != oe_ref.e()) {
	   assert(ccw_from.e() < nedgesl);
	   const int ne = e_t2c(ccw_from);
	   tEdge *edgPtr = EdgeTable[ne];
	   curnode->insertBackSpokeList( edgPtr );
	   ccw_from = ccw_from.ccw_edge_around_from(edges);
	 }
       }
     while( nodIter.Next() );
     delete [] oedge;
   }
   cout << "done.\n";

   // Assign ccwedg connectivity (that is, tell each edge about its neighbor
   // immediately counterclockwise)
   cout << "Setting up CCW edges..." << flush;
   {
     int iedge;
     tEdge * curedg, * ccwedg;
     tMeshListIter< tEdge > edgIter( edgeList );
     for( iedge=0, curedg=edgIter.FirstP(); iedge<nedgesl; ++iedge)
       {
	 {
	   const oriented_edge e1(iedge,true);
	   const oriented_edge ccw_from = e1.ccw_edge_around_from(edges);
	   const int ccwedgid = e_t2c(ccw_from);
	   ccwedg = EdgeTable[ccwedgid];
	   curedg->setCCWEdg( ccwedg );
	 }
	 curedg = edgIter.NextP();
	 {
	   const oriented_edge e2(iedge,false);
	   const oriented_edge ccw_to = e2.ccw_edge_around_from(edges);
	   const int ccwedgid = e_t2c(ccw_to);
	   ccwedg = EdgeTable[ccwedgid];
	   curedg->setCCWEdg( ccwedg );
	 }
	 curedg = edgIter.NextP(); 
       }
   }
   cout << "done.\n";
   
   cout << "setting up triangle connectivity..." << flush;
   {
     int ielem;
     for ( ielem=0; ielem<nelem; ++ielem ) {
       if (0) // DEBUG
	 cout << "TRI " << ielem << endl << flush;
       tTriangle newtri;
       newtri.setID( ielem );
       if (0) { // DEBUG
	 cout << "p0=" << p[elems[ielem].p1].id << " "
	      << "(" << p[elems[ielem].p1].x
	      << "," << p[elems[ielem].p1].y << "), " 
	      << "p1=" << p[elems[ielem].p2].id << " "
	      << "(" << p[elems[ielem].p2].x
	      << "," << p[elems[ielem].p2].y << "), " 
	      << "p2=" << p[elems[ielem].p3].id << " "
	      << "(" << p[elems[ielem].p3].x
	      << "," << p[elems[ielem].p3].y << "), " 
	      << "e0=" << elems[ielem].e1 << " "
	      << "e1=" << elems[ielem].e2 << " "
	      << "e2=" << elems[ielem].e3 << endl;
       }
       {
	 newtri.setPPtr( 0, NodeTable[p[elems[ielem].p1].id] );
	 newtri.setPPtr( 1, NodeTable[p[elems[ielem].p2].id] );
	 newtri.setPPtr( 2, NodeTable[p[elems[ielem].p3].id] );
       }
       {
	 newtri.setEPtr( 0, EdgeTable[e_t2c(elems[ielem].e1, 
					    elems[ielem].eo1)] );
	 newtri.setEPtr( 1, EdgeTable[e_t2c(elems[ielem].e2,
					    elems[ielem].eo2)] );
	 newtri.setEPtr( 2, EdgeTable[e_t2c(elems[ielem].e3,
					    elems[ielem].eo3)] );
       }
       triList.insertAtBack( newtri );
     }
     const tIdArray< tTriangle > TriTable(triList); // for fast lookup per ID
     
     tTriangle * ct, * nbrtri;
     tListIter< tTriangle > triIter( triList );
     for( ielem=0, ct=triIter.FirstP(); ielem<nelem; ct=triIter.NextP(), 
	    ++ielem ) {
       nbrtri = ( elems[ielem].t1>=0 ) ? TriTable[ elems[ielem].t1 ] : 0;
       ct->setTPtr( 0, nbrtri );
       nbrtri = ( elems[ielem].t2>=0 ) ? TriTable[ elems[ielem].t2 ] : 0;
       ct->setTPtr( 1, nbrtri );
       nbrtri = ( elems[ielem].t3>=0 ) ? TriTable[ elems[ielem].t3 ] : 0;
       ct->setTPtr( 2, nbrtri );
     }
   }   
   cout<<"done.\n";

   // deallocation of Tipper triangulator data structures
   delete [] edges;
   delete [] elems;
   delete [] p;
   
   // assertions
   assert( edgeList.getSize() == 2*nedgesl );
   assert( triList.getSize() == nelem );

   //
   UpdateMesh(); //calls CheckMeshConsistency()  TODO: once bug-free,
   CheckMeshConsistency();                     //remove CMC call from UM
   cout<<"MakeMeshFromScratchTipper done.\n";
}


