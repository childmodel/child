/***************************************************************************\
**
**  tGrid.cpp: Functions for class tGrid
**
**  $Id: tMesh.cpp,v 1.49 1998-10-01 18:46:16 nmgaspar Exp $
\***************************************************************************/

#include "tGrid.h"


     //global functions:


//global function; determines whether test point violates "Delaunay-ness"
//of other three; i.e., does the triangle 'pass the test' against the other
int TriPasses( tArray< double > &ptest,
               tArray< double > &p0,
               tArray< double > &p1,
               tArray< double > &p2 )
{
   assert( (&ptest != 0) && (&p0 != 0) && (&p1 != 0) && (&p1 != 0) );
     //cout << "TriPasses? ";
   double dx0, dx1, dy0, dy1;
   double crossp, dotp, angle0_2_1, angle0_test_1;
   dx0 = p0[0] - p2[0];
   dx1 = p1[0] - p2[0];
   dy0 = p0[1] - p2[1];
   dy1 = p1[1] - p2[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   dotp = dx0 * dx1 + dy0 * dy1;
   angle0_2_1 = atan2( crossp, dotp );
   dx0 = p0[0] - ptest[0];
   dx1 = p1[0] - ptest[0];
   dy0 = p0[1] - ptest[1];
   dy1 = p1[1] - ptest[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   dotp = dx0 * dx1 + dy0 * dy1;
   angle0_test_1 = atan2( crossp, dotp );
   //determine:
   if( angle0_2_1 < angle0_test_1 )
   {
        //cout << "Yes" << endl;
      return 0;
   }
   else
   {
        //cout << "No" << endl;
      return 1;
   }
}

//global function; determines whether points are counter-clockwise:
int PointsCCW( tArray< double > &p0,
               tArray< double > &p1,
               tArray< double > &p2 )
{
   assert( &p0 != 0 && &p1 != 0 && &p1 != 0 );
     //cout << "PointsCCW? ";
   int i;
   double dx0, dx1, dy0, dy1;
   double crossp;
   dx0 = p1[0] - p0[0];
   dx1 = p2[0] - p0[0];
   dy0 = p1[1] - p0[1];
   dy1 = p2[1] - p0[1];
   crossp = dx0 * dy1 - dx1 * dy0;
   if( crossp > 0 )
   {
        //cout << "Yes" << endl;
      return 1;
   }
   else
   {
        //cout << "No" << endl;
      return 0;
   }
}

//global function; determines whether triangle's "new" coords are CCW
int NewTriCCW( tTriangle *ct )
{
   assert( ct != 0 );
   tLNode *cn;
   cn = (tLNode *) ct->pPtr(0);
   tArray< double > p0( cn->get2DCoords() );
   if( cn->Meanders() ) p0 = cn->getNew2DCoords();
   cn = (tLNode *) ct->pPtr(1);
   tArray< double > p1( cn->get2DCoords() );
   if( cn->Meanders() ) p1 = cn->getNew2DCoords();
   cn = (tLNode *) ct->pPtr(2);
   tArray< double > p2( cn->get2DCoords() );
   if( cn->Meanders() ) p2 = cn->getNew2DCoords();
   if( PointsCCW( p0, p1, p2 ) ) return 1;
   else return 0;
}

//global function; determines whether coords are in "new" triangle
int InNewTri( tArray< double > &xy, tTriangle *ct )
{
   int j;
   tLNode *vtx;
   tArray< double > xy1, xy2;
   for( j=0; j<3; j++ )
   {
      vtx = (tLNode *) ct->pPtr(j);
      if( vtx->Meanders() ) xy1 = vtx->getNew2DCoords();
      else xy1 = vtx->get2DCoords();
      vtx = (tLNode *) ct->pPtr( (j+1)%3 );
      if( vtx->Meanders() ) xy2 = vtx->getNew2DCoords();
      else xy2 = vtx->get2DCoords();
      if ( ( (xy1[1] - xy[1]) * (xy2[0] - xy[0]) ) >
           ( (xy1[0] - xy[0]) * (xy2[1] - xy[1])) )
          break;
   }
   if( j == 3) return 1;
   else return 0;
}


//global function; determines whether nbr node currently pointed to
//by iterator and the next two in the nbr list form a Delaunay triangle:
template< class tSubNode >
int Next3Delaunay( tPtrList< tSubNode > &nbrList,
                   tPtrListIter< tSubNode > &nbrIter )
{
   static ncalls = 0;
   ncalls++;
   
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   tSubNode *cn, *nbrnd;
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
   tPtrListNode< tSubNode > *tempptrnode = nbrIter.NodePtr();
   tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   tArray< double > p1( nbrIterCopy.NextP()->get2DCoords() );
   tArray< double > p2( nbrIterCopy.NextP()->get2DCoords() );
     //cout << "N3D: point B\n";
   if( !PointsCCW( p0, p1, p2 ) ) return 0;
   tArray< double > ptest;
   cn = nbrIterCopy.NextP();
   while( cn != nbrnd )
   {
      ptest = cn->get2DCoords();
      if( !TriPasses( ptest, p0, p1, p2 ) )
      {
           //cout << "Next3Delaunay? No" << endl;
         return 0;
      }
        //else cout << "Next3Del? this tri passed..\n";
      
      cn = nbrIterCopy.NextP();
   }
   //cout << "Next3Delaunay? Yes" << endl;
   return 1;
}

//global function; determines whether nbr node currently pointed to
//by iterator and the next two in the nbr list form a Delaunay triangle:
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
   tPtrListNode< tSubNode > *tempptrnode = nbrIter.NodePtr();
   tArray< double > p0( nbrIterCopy.DatPtr()->get2DCoords() );
   assert( nbrIterCopy.Next() );
   tArray< double > p1( nbrIterCopy.DatPtr()->get2DCoords() );
   tArray< double > p2( testNode.get2DCoords() );
   if( !PointsCCW( p0, p1, p2 ) ) return 0;
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

/*****************************************************************************\
**
**      Intersect: uses newx, newy
**      Data members updated: Grid
**      Called by: 
**      Calls:  
**
\*****************************************************************************/
int Intersect( tEdge * ae, tEdge * be )
{
   //cout << "Intersect(...)..." << endl;
   double dxa, dxb, dya, dyb,
       a, b, c, f, g, h, x, y,
       rangemina, rangeminb, rangemaxa, rangemaxb, rangemin, rangemax;
   double smallnum = 0.0000000001;
   tLNode * lnode;
   if( !ae || !be )
   {
      cout<<"Intersect: Warning: invalid edge(s)"<<endl<<flush;
      return( 0 );
   }
   if( !ae->getOriginPtr() || !ae->getDestinationPtr() ||
       !be->getOriginPtr() || !be->getOriginPtr() )
   {
      cout<<"Intersect: Warning: invalid org or dest"<<endl<<flush;
      return( 0 );
   }
   lnode = (tLNode *) ae->getOriginPtrNC();
   tArray< double > aoc( lnode->get2DCoords() );
   if( lnode->Meanders() ) aoc = lnode->getNew2DCoords();
   lnode = (tLNode *) ae->getDestinationPtrNC();
   tArray< double > adc( lnode->get2DCoords() );
   if( lnode->Meanders() ) adc = lnode->getNew2DCoords();
   lnode = (tLNode *) be->getOriginPtrNC();
   tArray< double > boc( lnode->get2DCoords() );
   if( lnode->Meanders() ) boc = lnode->getNew2DCoords();
   lnode = (tLNode *) be->getDestinationPtrNC();
   tArray< double > bdc( lnode->get2DCoords() );
   if( lnode->Meanders() ) bdc = lnode->getNew2DCoords();
   
//   cout << "Intersect: edge endpoints, edge " << ae->getID()
//        << " org " << aoc[0] << " " << aoc[1] << ", dest " << adc[0]
//        << " " << adc[1] << "; edge " << be->getID() << " org " << boc[0]
//        << " " << boc[1] << ", dest " << bdc[0] << " " << bdc[1] << endl << flush;
   
   dxa = adc[0] - aoc[0] + smallnum;
   dxb = bdc[0] - boc[0] + smallnum;
   dya = adc[1] - aoc[1] + smallnum;
   dyb = bdc[1] - boc[1] + smallnum;
   a = dya;
   b = -dxa;
   c = dxa * aoc[1] - dya * aoc[0];
   f = dyb;
   g = -dxb;
   h = dxb * boc[1] - dyb * boc[0];
   if( fabs(dxa) > 0 && fabs(dxb) > 0 )
   {
      if( fabs(f - g * a / b) > smallnum )
      {
         x = (g * c / b - h) / (f - g * a / b);
         y = (-c - a * x) / b;
         if( adc[0] < aoc[0] )
         {
            rangemina = adc[0];
            rangemaxa = aoc[0];
         }
         else
         {
            rangemina = aoc[0];
            rangemaxa = adc[0];
         }
         if( bdc[0] < boc[0] )
         {
            rangeminb = bdc[0];
            rangemaxb = boc[0];
         }
         else
         {
            rangeminb = boc[0];
            rangemaxb = bdc[0];
         }
         if( rangemina < rangeminb ) rangemin = rangeminb + smallnum;
         else rangemin = rangemina + smallnum;
         if( rangemaxa > rangemaxb ) rangemax = rangemaxb - smallnum;
         else rangemax = rangemaxa - smallnum;
         if( x > rangemin && x < rangemax )
             return( 1 );
         else
             return( 0 );
      }
      else return (0);
   }
   else
   {
      if( fabs(dya) > 0 && fabs(dyb) > 0 )
      {
         if( fabs(g - f * b / a) > smallnum )
         {
            y = (f * c / a - h) / (g - f * b / a);
            x = (-c - b * y) / a;
            if( adc[1] < aoc[1] )
            {
               rangemina = adc[1];
               rangemaxa = aoc[1];
            }
            else
            {
               rangemina = aoc[1];
               rangemaxa = adc[1];
            }
            if( bdc[1] < boc[1] )
            {
               rangeminb = bdc[1];
               rangemaxb = boc[1];
            }
            else
            {
               rangeminb = boc[1];
               rangemaxb = bdc[1];
            }
            if( rangemina < rangeminb ) rangemin = rangeminb + smallnum;
            else rangemin = rangemina + smallnum;
            if( rangemaxa > rangemaxb ) rangemax = rangemaxb - smallnum;
            else rangemax = rangemaxa - smallnum;
            if( y > rangemin && y < rangemax ) return( 1 );
            else return( 0 );
         }
         else return 0;
      }
      else //one horiz. line and one vert. line:
      {
         if( fabs(dya) == 0 )
         {
            x = boc[0];
            y = aoc[1];
            if( adc[0] < aoc[0] )
            {
               rangemina = adc[0] + smallnum;
               rangemaxa = aoc[0] - smallnum;
            }
            else
            {
               rangemina = aoc[0] + smallnum;
               rangemaxa = adc[0] - smallnum;
            }
            if( bdc[1] < boc[1] )
            {
               rangeminb = bdc[1] + smallnum;
               rangemaxb = boc[1] - smallnum;
            }
            else
            {
               rangeminb = boc[1] + smallnum;
               rangemaxb = bdc[1] - smallnum;
            }
            if( x > rangemina && x < rangemaxa && y > rangeminb && y < rangemaxb )
                return 1;
            else return 0;
         }
         else
         {
            x = aoc[0];
            y = boc[1];
            if( adc[1] < aoc[1] )
            {
               rangemina = adc[1] + smallnum;
               rangemaxa = aoc[1] - smallnum;
            }
            else
            {
               rangemina = aoc[1] + smallnum;
               rangemaxa = adc[1] - smallnum;
            }
            if( bdc[0] < boc[0] )
            {
               rangeminb = bdc[0] + smallnum;
               rangemaxb = boc[0] - smallnum;
            }
            else
            {
               rangeminb = boc[0] + smallnum;
               rangemaxb = bdc[0] - smallnum;
            }
            if( y > rangemina && y < rangemaxa && x > rangeminb && x < rangemaxb )
                return 1;
            else return 0;
         }
      }
   }
}

/**************************************************************************\
**
**         Utilities for class tGrid
**
\**************************************************************************/
//default constructor
template< class tSubNode >
tGrid< tSubNode >::
tGrid() {nnodes = nedges = ntri = seed = 0;cout<<"tGrid()"<<endl;}     //tGrid

/**************************************************************************\
**
**   tGrid( infile ): Reads from infile whether it is to reconstruct a grid
**                    from input, construct a grid from a list of (x,y,z,b)
**                    points (b=boundary code), or construct a grid from
**                    scratch, given parameters in infile.
**
**   Created: 2/11/98, SL
**   Modified: 4/98 GT added MakeGridFromPoints function
**   Calls: tInputFile::ReadItem,
**        MakeGridFromInputData( infile ), MakeGridFromScratch( infile ),
**        of MakeGridFromPoints( infile )
**   Parameters: infile -- main input file containing option for input
**                         reading and any other needed parameters (but
**                         not mesh point coords and connectivity data;
**                         if needed, these are in separate files
**   Notes: needs to find 0, 1 or 2 under the heading of "OPTREADINPUT"
**               in infile.
**   Added 7/98 - will read in layering information from a Child output file
**         if OPTREADLAYER is set to 1.  
**
\**************************************************************************/
template< class tSubNode >
tGrid< tSubNode >::
tGrid( tInputFile &infile )
{
   int read = infile.ReadItem( read, "OPTREADINPUT" );
   if( read<0 || read>2 )
   {
      cerr << "Valid options for reading mesh input are:\n"
           << "  0 -- create rectangular offset mesh\n"
           << "  1 -- read mesh from input data files\n"
           << "  2 -- create mesh from a list of (x,y,z,b) points\n";
      ReportFatalError( "Invalid mesh input option requested." );
   }
   
   if( read==1 ){
      MakeGridFromInputData( infile ); //create grid by reading data files
      int lay = infile.ReadItem( lay, "OPTREADLAYER" );
      if( lay == 1 )
          MakeLayersFromInputData( infile );
   }
   else if( read==2 )
       MakeGridFromPoints( infile );  //create new mesh from list of points
   else
       MakeGridFromScratch( infile ); //create new grid with parameters

   layerflag=infile.ReadItem( layerflag, "OPTINTERPLAYER");
   
}

//destructor
template< class tSubNode >
tGrid< tSubNode >::
~tGrid() {cout << "    ~tGrid()" << endl;}                    //tGrid


/************************************************************************\
  MakeLayersFromInputData( infile ):

  This function reads in layer information from a CHILD output file and
  copies it to the nodes which have already been created.

\************************************************************************/

template< class tSubNode >
void tGrid< tSubNode >::
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
         cout << "from file, time = " << time << endl;
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
   int numg = infile.ReadItem( numg, "NUMGRNSIZE" );
   layhelp.setDgradesize(numg);

   int g;
   tLNode * cn;
   int nActNodes = getNodeList()->getActiveSize();
   //int NNodes = getNodeList()->get
   tGridListIter<tLNode> ni( getNodeList() );

   for( cn = ni.FirstP(); ni.IsActive(); cn = ni.NextP() )
   {
      layerinfile >> numl;
      for(i = 1; i<=numl; i++){
         layerinfile >> ditem;
         layhelp.setCtime(ditem);
         layerinfile >> ditem;
         layhelp.setRtime(ditem);
         layerinfile >> item;
         layhelp.setFlag(item);
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
   layhelp.setFlag(0);
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
**   MakeGridFromInputData( infile ):
**          formerly tGrid( tListInputData &, tlnodflag ). Constructs
**          tListInputData obj., makes grid from data in that object.
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
**
\**************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
MakeGridFromInputData( tInputFile &infile )
{
   int i;
   tListInputData< tSubNode > input( infile );
   seed = 0;
   // set the number of nodes, edges, and triangles in the grid mesh
   //assert( lnodflag );
   nnodes = input.x.getSize();
   nedges = input.orgid.getSize();
   ntri = input.p0.getSize();
   cout << "nnodes, nedges, ntri: " << nnodes << " " << nedges << " " << ntri << endl << flush;
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
   tGridListIter< tSubNode > nodIter( nodeList );
   tEdge tempedge1, tempedge2;
   int obnd, dbnd;
   for( i = 0; i< nedges-1; i+=2 )
   {
      // Assign values: ID, origin and destination pointers
      tempedge1.setID( i );
      tempedge2.setID( i + 1 );
      //cout << input.orgid[i] << " " << input.destid[i] << endl;
      //cout << nodIter.Get( input.orgid[i] ) << " ";
      //cout << nodIter.Get( input.destid[i] ) << endl;
      assert( nodIter.Get( input.orgid[i] ) );
          //{
         tempedge1.setOriginPtr( &(nodIter.DatRef()) );
         tempedge2.setDestinationPtr( &(nodIter.DatRef()) );
         obnd = nodIter.DatRef().getBoundaryFlag();
         //cout << nodIter.DatRef().getID() << "->";
         //}
         assert( nodIter.Get( input.destid[i] ) );
          //{
         tempedge1.setDestinationPtr( &(nodIter.DatRef()) );
         tempedge2.setOriginPtr( &(nodIter.DatRef()) );
         dbnd = nodIter.DatRef().getBoundaryFlag();
         //cout << nodIter.DatRef().getID() << endl;
         //}
      // set the "flowallowed" status (FALSE if either endpoint is a
         // closed boundary) and insert edge pair onto the list --- active
         // part of list if flow is allowed, inactive if not
         //cout << "BND: " << obnd << " " << dbnd << " " << kClosedBoundary
         //     << endl;
      if( obnd == kClosedBoundary || dbnd == kClosedBoundary )
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

   // set up the lists of edges (spokes) connected to each node
   // (GT added code to also assign the 1st edge to "edg" as an alternative
   // to spokelist implementation)
   cout << "setting up spoke lists..." << flush;
   int e1;
   int ne;
   tGridListIter< tEdge >
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
   tGridListIter<tEdge> ccwIter( edgeList ); // 2nd iter for performance
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
   tTriangle temptri;
   for ( i=0; i<ntri; i++ )
   {
      //cout << "TRI " << i << endl << flush;
      temptri.setID( i );
      if( nodIter.Get( input.p0[i] ) )
          temptri.setPPtr( 0, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p1[i] ) )
          temptri.setPPtr( 1, &(nodIter.DatRef()) );
      if( nodIter.Get( input.p2[i] ) )
          temptri.setPPtr( 2, &(nodIter.DatRef()) );
      if( edgIter.Get( input.e0[i] ) )
          temptri.setEPtr( 0, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e1[i] ) )
          temptri.setEPtr( 1, &(edgIter.DatRef()) );
      if( edgIter.Get( input.e2[i] ) )
          temptri.setEPtr( 2, &(edgIter.DatRef()) );
      triList.insertAtBack( temptri );
   }
   
   tListIter< tTriangle >
       trIter1( triList );
   tListIter< tTriangle >
       trIter2( triList );
   i = 0;
   do                                                //for ( i=0; i<ntri; i++ )
   {
      //cout << "tGrid: Tri loop, i: " << i << endl;
      if( trIter2.Get( input.t0[i] ) )
          trIter1.DatRef().setTPtr( 0, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 0, 0 );
      if( trIter2.Get( input.t1[i] ) )
          trIter1.DatRef().setTPtr( 1, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 1, 0 );
      if( trIter2.Get( input.t2[i] ) )
          trIter1.DatRef().setTPtr( 2, &(trIter2.DatRef()) );
      else trIter1.DatRef().setTPtr( 2, 0 );
      i++;
   }
   while( trIter1.Next() );
   cout<<"done.\n";

   UpdateMesh();
   CheckMeshConsistency();
   cout << "end of tGrid( input )" << endl << flush;
}


/**************************************************************************\
**
**   MakeGridFromScratch( infile ):formerly tGrid( infile ). Makes
**       grid from scratch; reads parameters
**        from input file to get grid size, spacing, method of point
**        placement.
**       Could probably be done more gracefully, but here's how it does it:
**        1) makes boundary nodes and edges between them;
**        2) makes triangulation with only boundary nodes;
**        3) adds the rest of the nodes one at a time, i.e., triangulation
**           redone for each node added.
**
**       A nicer way (who wants to code it? oo-oo! stephen does, stephen does!):
**        1) make boundary nodes and edges between them;
**        2) add all other nodes at once;
**        3) do triangulation once by making a temporary list of boundary
**           edges and, for each edge, finding the node which would make
**           the biggest angle with the endpoints of the edge; remove edge(s)
**           from temp. list, make new edge(s) and triangle; add new edge(s)
**           to temp. list. (from Du)
**                    
**
**   Created: 2/11/98, SL
**   Calls: tInputFile::ReadItem, MakeCCWEdges(),
**          UpdateMesh(), CheckMeshConsistency()
**   Parameters: xGrid, yGrid, boundType, mElev, ptPlace, delGrid, numPts,
**               upperZ, xout, yout
**   Modified: 3/13/98--now makes rows offset so pattern is "hexagonal"
**      rather than square
**               
**
\**************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
MakeGridFromScratch( tInputFile &infile )
{
   int i, j, id, n, nx, ny;
   int numPts;
   double dist;
   double delGrid, slope;
   double upperZ;
   tArray< double > xyz(3);
   cout << "In MGFS, calling node constr w/ infile\n";
   tSubNode tempnode( infile ),  // temporary node used to create node list
       *cn, *node0, *node1, *node2, *node3;
   tEdge *ce;
   tGridListIter< tEdge > edgIter( edgeList );
   tGridListIter< tSubNode > nodIter( nodeList );
   tPtrList< tSubNode > bndList;
   seed = infile.ReadItem( seed, "SEED" );
     //reads in size of grid (meters)
   double xGrid = infile.ReadItem( xGrid, "X_GRID_SIZE" );
   double yGrid = infile.ReadItem( yGrid, "Y_GRID_SIZE" );
     //read type of open boundary:
     //  0 = corner outlet (lower left)
     //  1 = one open side (lower)
     //  2 = two opposite sides (upper and lower)
     //  3 = all sides
     //  4 = specify outlet coordinates
   int boundType = infile.ReadItem( boundType, "TYP_BOUND" );
     //read mean elevation (also determines elev variation)
   double mElev = infile.ReadItem( mElev, "MEAN_ELEV" );
     //reads method of point placement:
     //  0 = uniform grid;
     //  1 = perturbed grid;
     //  2 = random scatter;
   int ptPlace = infile.ReadItem( ptPlace, "OPT_PT_PLACE" );
     //read grid spacing or number of points (excluding four boundary points)
   if( ptPlace == kUniformGrid || ptPlace == kPerturbedGrid )
   {
      delGrid = infile.ReadItem( delGrid, "GRID_SPACING" );
   }
   else
   {
      numPts = infile.ReadItem( numPts, "NUM_PTS" );
      delGrid = sqrt( xGrid * yGrid / numPts );
   }

   //MAKE BOUNDARY
   if( boundType == kCornerOutlet )
   {
      id = 0;
      tempnode.setBoundaryFlag( kOpenBoundary );
      tempnode.set3DCoords( 0, 0, 0 );
      tempnode.setID( id );
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      nodeList.insertAtBack( tempnode );
      bndList.insertAtBack( nodIter.LastP() );
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=1, id++; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( boundType == kOpenSide )
   {
      cout << "OPEN SIDE boundary\n";
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=1, id=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( i - 1 );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=n; i>=0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   if( boundType == kOppositeSidesOpen )
   {
      upperZ = infile.ReadItem( upperZ, "UPPER_BOUND_Z" );
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=1, id=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = yGrid / delGrid;
      for( i=0; i<=n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      tempnode.setBoundaryFlag( kOpenBoundary );
      n = xGrid / delGrid;
      for( i=n-1; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, upperZ );
         tempnode.setID( id );
         nodeList.insertAtBoundFront( tempnode );
         bndList.insertAtBack( nodIter.GetP( id ) );
      }
      tempnode.setBoundaryFlag( kClosedBoundary );
      n = yGrid / delGrid;
      for( i=n; i>=0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( boundType == kAllSidesOpen )
   {
      id = 0;
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kOpenBoundary );
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
      n = yGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
      }
   }
   else if( boundType == kSpecifyOutlet )
   {
      double xout = infile.ReadItem( xout, "OUTLET_X_COORD" );
      double yout = infile.ReadItem( yout, "OUTLET_Y_COORD" );
      n = xGrid / delGrid;
      tempnode.setBoundaryFlag( kClosedBoundary );
      for( i=0, id=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.01 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, 0, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( yout == 0 && xout > dist && xout < dist + delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
      n = yGrid / delGrid;
      for( i=0; i<n; i++, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( xGrid, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( xout == xGrid && yout > dist && yout < dist + delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
      n = xGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( dist, yGrid, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( yout == yGrid && xout < dist && xout > dist - delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
      n = yGrid / delGrid;
      for( i=n; i>0; i--, id++ )
      {
         dist = i * delGrid + 0.0001 * delGrid * ( ran3( &seed ) - 0.5 );
         tempnode.set3DCoords( 0, dist, 0 );
         tempnode.setID( id );
         nodeList.insertAtBack( tempnode );
         bndList.insertAtBack( nodIter.LastP() );
         if( xout == 0 && yout < dist && yout > dist - delGrid )
         {
            tempnode.set3DCoords( xout, yout, 0 );
            tempnode.setBoundaryFlag( kOpenBoundary );
            id++;
            tempnode.setID( id );
            nodeList.insertAtBoundFront( tempnode );
            bndList.insertAtBack( nodIter.GetP( id ) );
            tempnode.setBoundaryFlag( kClosedBoundary );
         }
      }
   }
   bndList.makeCircular();
   cout << "made points; now adding edges\n";
   
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
   cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }
   cout << "calling repair mesh for initial boundary\n";
   int meshok = RepairMesh( bndList );
   assert( meshok );
   cout << "filling in points\n";
   
     //FILL IN POINTS
   //gt modified to use AddNode instead of AddNodeAt
   tempnode.setBoundaryFlag( kNonBoundary );
   if( ptPlace == kUniformGrid || ptPlace == kPerturbedGrid )
   {
      nx = xGrid / delGrid;
      ny = yGrid / delGrid;
      for( i=1; i<nx; i++ )
      {
         for( j=1; j<ny; j++, id++ )
         {
            //rows are offset such that there should be an
            //edge leading to a corner outlet
            xyz[0] = i * delGrid - 0.25 * delGrid * (j%2)
                + 0.25 * delGrid * ((j+1)%2);
            xyz[1] = j * delGrid;
            if( ptPlace == kPerturbedGrid )
            {
               xyz[0] += 0.01 * delGrid * ( ran3( &seed ) - 0.5 );
               xyz[1] += 0.01 * delGrid * ( ran3( &seed ) - 0.5 );
            }
            xyz[2] = mElev + mElev * ( ran3( &seed ) - 0.5 );
            if( boundType == kOppositeSidesOpen )
            {
               slope = upperZ / yGrid;
               xyz[2] += slope * xyz[1] - mElev;
            }
            tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
            tempnode.setID( id );
            AddNode( tempnode );
            //XAddNodeAt( xyz );
         }
      }
   }
   else if( ptPlace == kRandomGrid )
   {
      for( i=0; i<numPts; i++ )
      {
         xyz[0] = ran3(&seed) * xGrid;
         xyz[1] = ran3(&seed) * yGrid;
         xyz[2] = mElev + mElev * ( ran3( &seed ) - 0.5 );
         if( xyz[0] != 0 && xyz[0] != xGrid && xyz[1] != 0 && xyz[1] != yGrid )
         {
            tempnode.set3DCoords( xyz[0], xyz[1], xyz[2] );
            tempnode.setID( id );
            AddNode( tempnode );
            id++;
         }
            //XAddNodeAt( xyz );
      }
   }
   MakeCCWEdges();
   UpdateMesh(); //calls CheckMeshConsistency()  TODO: once bug-free,
   CheckMeshConsistency();                     //remove CMC call from UM
}


/**************************************************************************\
**
**   MakeGridFromPoints
**
**   Constructs a mesh from a given set of (x,y,z,b) points, where
**   b=boundary code. Unlike MakeGridFromInputData no connectivity
**   information is needed, just the coordinates and boundary codes.
**
**   The format of the file containing points is ...
**
**   Variations: to reduce memory overhead, the routine could be modified
**   to read on (x,y,z,b) set at a time rather than reading and
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
void tGrid< tSubNode >::
MakeGridFromPoints( tInputFile &infile )
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
   dx = maxx - minx;
   dy = maxy - miny;

   // Create the 3 nodes that form the supertriangle and place them on the
   // node list in counter-clockwise order. (Note that the base and height
   // of the supertriangle are 5 times the
   // width and height, respectively, of the rectangle that encloses the
   // points.) Assigning the IDs allows us to retrieve and delete these
   // nodes when we're done creating the mesh.
   tempnode.set3DCoords( minx-2*dx, miny-2*dy, 0.0 );
   tempnode.setBoundaryFlag( kClosedBoundary );
   tempnode.setID( -1 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( maxx+2*dx, miny-2*dy, 0.0 );
   tempnode.setID( -2 );
   nodeList.insertAtBack( tempnode );
   tempnode.set3DCoords( minx+0.5*dx, maxy+2*dy, 0.0 );
   tempnode.setID( -3 );
   nodeList.insertAtBack( tempnode );

   // set # of nodes, edges, and triangles
   nnodes = 3;
   nedges = ntri = 0;

   // Create the edges that connect the supertriangle vertices and place
   // them on the edge list.
   // (To do this we need to retrieve pointers from the nodeList)
   tGridListIter<tSubNode> nodIter( nodeList );
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

   //cout << "1 NN: " << nnodes << " (" << nodeList.getActiveSize() << ")  NE: " << nedges << " NT: " << ntri << endl;

   // Now add the points one by one to construct the mesh.
   for( i=0; i<numpts; i++ )
   {
      //cout << "IN MGFP, ADDING NODE " << i << endl;
      tempnode.setID( i );
      tempnode.set3DCoords( x[i], y[i], z[i] );
      tempnode.setBoundaryFlag( bnd[i] );
      AddNode( tempnode );
   }

   //cout << "\n2 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: " << nedges << " NT: " << ntri << endl;

   // We no longer need the supertriangle, so remove it by deleting its
   // vertices.
   DeleteNode( stp1, kNoRepair );
   DeleteNode( stp2, kNoRepair );
   DeleteNode( stp3, kNoRepair );
   
   //cout << "3 NN: " << nnodes << " (" << nodeList.getActiveSize() << ") NE: " << nedges << " NT: " << ntri << endl;
   
   // Update Voronoi areas, edge lengths, etc., and test the consistency
   // of the new mesh.
   UpdateMesh();
   CheckMeshConsistency( );

   // Clean up (probably not strictly necessary bcs destructors do the work)
   supertriptlist.Flush();
   
}


/*****************************************************************************\
**
**  CheckMeshConsistency
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
void tGrid< tSubNode >::
CheckMeshConsistency( int boundaryCheckFlag ) /* default: TRUE */
{
   tGridListIter<tSubNode> nodIter( nodeList );
   tGridListIter<tEdge> edgIter( edgeList );
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
   cout << "MESH PASSED\n";
   return;
   
  error:
   ReportFatalError( "Error in mesh consistency." );
   
}
#undef kMaxSpokes

template< class tSubNode >
void tGrid< tSubNode >::
Print()                                                  //tGrid
{
   triList.print();
   nodeList.print();
   edgeList.print();
}

template< class tSubNode >
void tGrid< tSubNode >::
MakeCCWEdges()
{
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn;
   for( cn = nodIter.FirstP(); !( nodIter.AtEnd() ); cn = nodIter.NextP() )
   {
      cn->makeCCWEdges();
   }
}


/*****************************************************************************\
**
**  tGrid::setVoronoiVertices
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
**     - copied function to tGrid member from tStreamNet member, gt 3/98.
**       Other fns now use "right-hand" definition; this fn may have to
**       be changed.
**
\*****************************************************************************/
template <class tSubNode>
void tGrid<tSubNode>::setVoronoiVertices()
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
**  tGrid::CalcVoronoiEdgeLengths
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
void tGrid<tSubNode>::CalcVoronoiEdgeLengths()
{
	tEdge *ce;
	double vedglen;
	tGridListIter<tEdge> edgIter( edgeList );

	for( ce=edgIter.FirstP(); edgIter.IsActive(); ce=edgIter.NextP() )
	{
		vedglen = ce->CalcVEdgLen();  // Compute Voronoi edge length
		ce = edgIter.NextP();         // Advance to complement edge and
		ce->setVEdgLen( vedglen );    //   and assign the same edge length.
	}
}


// NB: TODO: move to tGrid
// TODO: change L-hand to R-hand orientation of Voronoi vertices and
// create faster Voronoi edge length computation scheme
template <class tSubNode>
 void tGrid<tSubNode>::CalcVAreas()
{
   //cout << "CalcVAreas()..." << endl << flush;
   double area;
   tLNode * curnode; //X *cn;
   tEdge *ce;
   tArray< double > xy;
   tGridListIter<tLNode> nodIter( nodeList );
   //XnI( nodeList );
   //XtPtrListIter< tEdge > sI;
   
   for( curnode = nodIter.FirstP(); nodIter.IsActive();
        curnode = nodIter.NextP() )
   {
         //area = VoronoiArea( curnode );
      curnode->ComputeVoronoiArea();
      
/*      for( cn = nI.FirstP(); nI.IsActive(); cn = nI.NextP() )
      {
         sI.Reset( cn->getSpokeListNC() );
         cout << "node " << cn->getID() << ": ";
         for( ce = sI.FirstP(); !(sI.AtEnd()); ce = sI.NextP() )
         {
            xy = ce->getRVtx();
            cout << xy[0] << " " << xy[1] << "; ";
         }
         cout << endl;
      }*/
         /*assert( area > 0 );
         curnode->setVArea( area );
         curnode->setVArea_Rcp( 1.0 / area );*/
   }
   //cout << "CalcVAreas() finished" << endl;
}


/* TODO: consolidate w/ the other deletenode (one calls the other) */
template< class tSubNode >
int tGrid< tSubNode >::
DeleteNode( tListNode< tSubNode > *nodPtr, int repairFlag )
{
   //cout << "DeleteNode: " << nodPtr->getDataPtr()->getID() << endl;
   int i;
   tPtrList< tSubNode > nbrList;
     //tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode nodeVal;
   if( !( ExtricateNode( nodPtr->getDataPtrNC(), nbrList ) ) ) return 0;
   nbrList.makeCircular();
   if( nodPtr->getDataRefNC().getBoundaryFlag() )
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
   //    << nodeVal.getX() << ", " << nodeVal.getY() << "; " << endl;
   nnodes = nodeList.getSize();
   nedges = edgeList.getSize();
   ntri = triList.getSize();
   tPtrListIter< tSubNode > nbrIter( nbrList );
     /*cout << "leaving hole defined by " << endl << "   Node  x  y " << endl;
   for( i=0, nbrIter.First(); nbrIter.NextIsNotFirst(); i++ )
   {
      if( i>0 ) nbrIter.Next();
      cout << "   " << nbrIter.DatPtr()->getID() << "     "
           << nbrIter.DatPtr()->getX() << "  "
           << nbrIter.DatPtr()->getY() << endl;
   }*/
   if( !RepairMesh( nbrList ) ) return 0;
   //reset node id's
   tGridListIter< tSubNode > nodIter( nodeList );
   assert( nodIter.First() );
   i = 0;
   do
   {
      nodIter.DatRef().setID( i );
      i++;
   }
   while( nodIter.Next() );
   
   UpdateMesh();
     //cout << "Mesh repaired" << endl;
   return 1;
}


/**************************************************************************\
**
**  tGrid::DeleteNode( tSubNode *, int =1 )
**
**  Deletes a node from the mesh. This is done by first calling
**  ExtricateNode to detach the node by removing its edges and their
**  associated triangles, then removing the node from the nodeList.
**  Normally, RepairMesh is then called to retriangulate the "hole" left
**  behind in the mesh. (However, if the node was on the hull of the
**  mesh there's no "hole" to fix --- the caller is assumed to be smart
**  enough to recognize this case and let us know about it by setting
**  repairFlag to kNoRepair. This is the case, for example, when deleting
**  the nodes that form a "supertriangle" as in MakeGridFromPoints).
**
**  Once the mesh is repaired, the nodes are renumbered and as a safety
**  measure for debugging/testing purposes, UpdateMesh is called.
**
**  Data mbrs modified:  nnodes, nedges, and ntri are updated;
**                       the node is deleted from nodeList; other edges &
**                       triangles are removed and/or modified by
**                       ExtricateNode and RepairMesh (qv)
**  Calls:  tGrid::ExtricateNode, tGrid::RepairMesh, plus utility member
**               functions of tNode, tGridList, etc. (and temporarily,
**               UpdateMesh)
**  Returns:  error code: 0 if either ExtricateNode or RepairMesh fails,
**            1 otherwise.
**  Assumes:  
**  Notes:
**    - repairFlag defaults to TRUE if not specified
**    - if node is on the hull and repairFlag is TRUE, the result can be
**        an infinite loop in RepairMesh as it tries to mend a hole that
**        doesn't exist. This condition isn't tested for, so be careful.
**  Modifications: added repairFlag 4/98 GT
**
\**************************************************************************/
template< class tSubNode >
int tGrid< tSubNode >::
DeleteNode( tSubNode *node, int repairFlag )
{
   int i;
   tPtrList< tSubNode > nbrList;
   tListNode< tSubNode > *nodPtr;
   tGridListIter< tSubNode > nodIter( nodeList );
   nodIter.Get( node->getID() );
   tSubNode nodeVal;
   
   //cout << "DeleteNode: " << node->getID() << " at " << node->getX() << " "
   //    << node->getY() << " " << node->getZ() << endl;
   //assert( repairFlag || node->getBoundaryFlag()==kClosedBoundary );
   
   nodPtr = nodIter.NodePtr();
   if( !( ExtricateNode( node, nbrList ) ) ) return 0;
   nbrList.makeCircular();
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
   //   << nodeVal.getX() << ", " << nodeVal.getY() << "; " << endl;
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
       if( !RepairMesh( nbrList ) ) return 0;

   //reset node id's
   assert( nodIter.First() );
   i = 0;
   do
   {
      nodIter.DatRef().setID( i );
      i++;
   }
   while( nodIter.Next() );
     //cout << "Mesh repaired" << endl;
   /*tSubNode *cn;
   for( cn = nodIter.FirstP(); nodIter.IsActive(); cn = nodIter.NextP() )
   {
      cout << "node " << cn->getID() << endl;
   }*/

   UpdateMesh();
   return 1;
}


/**************************************************************************\
**
**  tGrid::ExtricateNode
**
**  Detaches a node from the mesh by deleting all of its edges (which in
**  turn removes the affected triangles). Returns a list of the node's
**  former neighbors by modifying the nbrList input parameter. Also
**  returns a code that indicates failure if the node still has a non-empty
**  spoke list after edge deletion.
**
**  Data mbrs modified:  nnodes; edges and triangles are removed from
**                       edgeList and triList
**  Calls:  tGrid::DeleteEdge and utility member functions of tNode,
**               tPtrList, tPtrListIter
**  Output:  list of node's (former) neighbors, in nbrList
**  Returns:  1 if all edges successfully deleted, 0 if not
**  Assumes:  
**  Notes:
**  Modifications: if node is a closed boundary, any of its neighbors that
**             are non-boundaries are switched to closed boundaries, so
**             that nodes along the edge of the domain (including nodes of
**             a "supertriangle" used in MakeGridFromPoints) may be removed
**             without causing errors, GT 4/98
**
\**************************************************************************/
template< class tSubNode >
int tGrid< tSubNode >::
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
      nbrPtr = ( tSubNode * ) ce->getDestinationPtrNC();
      nbrList.insertAtBack( nbrPtr );
      if( node->getBoundaryFlag()                      // If node is a bdy make
          && nbrPtr->getBoundaryFlag()==kNonBoundary )// sure nbrs are also
      {                                                // boundaries.
         nbrPtr->ConvertToClosedBoundary();
         nodeList.moveToBack( nbrPtr );
      }
      DeleteEdge( ce );
   }  
   nnodes--;
   if( node->getSpokeList().isEmpty() ) return 1;
   return 0;
}
   
template< class tSubNode >
int tGrid< tSubNode >::
DeleteEdge( tEdge * edgePtr )
{
   //cout << "DeleteEdge(...) " << edgePtr->getID() << endl;
   tEdge edgeVal1, edgeVal2;
   if( !ExtricateEdge( edgePtr ) ) return 0;
   //Note, extricate edge does not actually remove the edge from
   //the edge list, it only moves the two edges (one 'line' that
   //has two directions) to the end or front of the list.
   //These two edges are then actually removed here.
   if( edgePtr->getBoundaryFlag() )
   {
      assert( edgeList.removeFromBack( edgeVal1 ) );
      assert( edgeList.removeFromBack( edgeVal2 ) );
   }
   else
   {
      assert( edgeList.removeFromFront( edgeVal1 ) );
      assert( edgeList.removeFromFront( edgeVal2 ) );
   }
//    cout << "  edges " << edgeVal1.getID() << " and "
//         <<  edgeVal2.getID() << " between nodes "
//         << edgeVal1.getOriginPtr()->getID() << " and "
//         << edgeVal2.getOriginPtr()->getID() << " removed" << endl;
   
   if( &edgeVal1 == 0 || &edgeVal2 == 0 ) return 0;
   return 1;
}


template< class tSubNode >
int tGrid< tSubNode >::
ExtricateEdge( tEdge * edgePtr )
{
   //cout << "ExtricateEdge: " << edgePtr->getID() << endl;
   assert( edgePtr != 0 );
     //temporary objects:
   tEdge *tempedgePtr, *ce, *cce, *spk;
   tGridListIter< tEdge > edgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   tPtrList< tEdge > *spkLPtr;
   tListNode< tEdge > *listnodePtr;
   tTriangle triVal1, triVal2;
   tArray< tTriangle * > triPtrArr(2);
     //cout << "find edge in list; " << flush;
   ce = edgIter.GetP( edgePtr->getID() );  //NB: why necessary? isn't ce the
                                       // same as edgePtr??  Yes, why???

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
    tSubNode * nodece = (tSubNode *) ce->getOriginPtrNC();
    nodece->WarnSpokeLeaving( ce );
    tSubNode * nodecce = (tSubNode *) cce->getOriginPtrNC();
    nodecce->WarnSpokeLeaving( cce );
   
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
**  tGrid::LocateTriangle
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
tTriangle * tGrid< tSubNode >::
LocateTriangle( double x, double y )
{
   //cout << "\nLocateTriangle (" << x << "," << y << ")\n";
   int n, lv=0;
   tListIter< tTriangle > triIter( triList );  //lt
   tTriangle *lt = &(triIter.DatRef());
   double a, b, c;
   int online;
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

      /*cout << "find tri for point w/ x, y, " << x << ", " << y
           << "; no. tri's " << ntri << "; now at tri " << lt->getID() << endl;
      lt->TellAll();
      cout << flush;*/

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
                                    //if lt == 0, i.e., point is out of bounds,
                                    //and we don't want that;
                                    //calling code is built to deal with lt == 0.
   if( online != -1 && lt->pPtr(online)->getBoundaryFlag() != kNonBoundary &&
       lt->pPtr( (online+1)%3 )->getBoundaryFlag() != kNonBoundary ) //point on bndy
       return 0;
   //else cout << "location out of bounds\n";
   return(lt);
}

template< class tSubNode >
tTriangle * tGrid< tSubNode >::
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
      p1 = (tSubNode *) lt->pPtr(lv);
      if( p1->Meanders() ) xy1 = p1->getNew2DCoords();
      else xy1 = p1->get2DCoords();
      p2 = (tSubNode *) lt->pPtr( (lv+1)%3 );
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

template< class tSubNode >
tTriangle *tGrid< tSubNode >::
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

template< class tSubNode >
int tGrid< tSubNode >::
DeleteTriangle( tTriangle * triPtr )
{
   //cout << "DeleteTriangle(...) " << triPtr->getID() << endl;
   //triPtr->TellAll();
   tTriangle triVal;
   if( !ExtricateTriangle( triPtr ) ) return 0;
   if( !( triList.removeFromFront( triVal ) ) ) return 0;
   if( &triVal == 0 ) return 0;
   return 1;
}
      
template< class tSubNode >
int tGrid< tSubNode >::
ExtricateTriangle( tTriangle *triPtr )
{
   //cout << "ExtricateTriangle" << endl;
   tListIter< tTriangle > triIter( triList );
   tTriangle *ct;
   for( ct = triIter.FirstP(); ct != triPtr && !( triIter.AtEnd() );
        ct = triIter.NextP() );
   if( ( triIter.AtEnd() ) ) return 0;
   int i, j;
   for( i=0; i<3; i++ ) for( j=0; j<3; j++ )
       if( triPtr->tPtr(i) != 0 )
           if( triPtr->tPtr(i)->tPtr(j) == triPtr ) triPtr->tPtr(i)->setTPtr( j, 0 );
   triList.moveToFront( triIter.NodePtr() );
   ntri--;
   return 1;
}

/*template< class tSubNode >
int tGrid< tSubNode >::
MakeMesh( tPtrList< tSubNode > &bndList )
{
   tPtrListIter< tSubNode > bndIter( bndList );
   tPtrList< tSubNode > nptrList;
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn, *cnn;
   for( cn = nodIter.FirstP(); !( nodIter.AtEnd() ); cn = nodIter.NextP() )
   {
      nptrList.insertAtBack( cn );
   }
   for( cn = bndIter.FirstP(); bndIter.NodePtr() != 0; cn = bndIter.NextP() )
   {
      cnn = bndIter.ReportNextP();
   }
   
}*/

template< class tSubNode >
int tGrid< tSubNode >::
RepairMesh( tPtrList< tSubNode > &nbrList )
{
   assert( &nbrList != 0 );
   //cout << "RepairMesh: " << endl;
   if( nbrList.getSize() < 3 ) return 0;
   int flowflag, i, j;
   tSubNode * gridnodePtr;
   nbrList.makeCircular();
   tPtrListIter< tSubNode > nbrIter( nbrList );
   
   while( nbrList.getSize() > 3 )
   {
      //cout << "in loop, nbr size = " << nbrList.getSize() << endl;
      
      flowflag = 1;
      if( Next3Delaunay( nbrList, nbrIter ) ) //checks for ccw and Del.
      {
           //cout << "found 3 Delaun!\n";
         
         AddEdgeAndMakeTriangle( nbrList, nbrIter );
           //remove "closed off" pt
         nbrList.removeNext( gridnodePtr, nbrIter.NodePtr() );
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

//vertices of tri in ccw order; edges are added between node1 and node2
//TODO: comment/document this fn
template< class tSubNode >
int tGrid< tSubNode >::
AddEdge( tSubNode *node1, tSubNode *node2, tSubNode *node3 ) 
{
   assert( node1 != 0 && node2 != 0 && node3 != 0 );
   /*cout << "AddEdge" << 
     "between nodes " << node1->getID()
        << " and " << node2->getID() << " w/ ref to node " << node3->getID() << endl;*/
   int flowflag = 1;
   int i, j, newid;
   tEdge tempEdge1, tempEdge2;
   tEdge *ce, *le;
   tGridListIter< tEdge > edgIter( edgeList );
   tGridListIter< tSubNode > nodIter( nodeList );
   tPtrListIter< tEdge > spokIter;
        //DumpEdges:
     /*cout << "edges:" << endl;
   for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
   {
      cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
           << " to " << ce->getDestinationPtrNC()->getID() << endl;
   }*/

     //deal with new edges and new triangle:
   tempEdge1.setOriginPtr( node1 );   //set edge1 ORG
   tempEdge2.setDestinationPtr( node1 );//set edge2 DEST
   if( node1->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   tempEdge2.setOriginPtr( node2 );   //set edge2 ORG
   tempEdge1.setDestinationPtr( node2 );//set edge1 DEST
   if( node2->getBoundaryFlag() == kClosedBoundary ) flowflag = 0;
   if( node1->getBoundaryFlag()==kOpenBoundary     // Also no-flow if both
       && node2->getBoundaryFlag()==kOpenBoundary ) //  nodes are open bnds
       flowflag = 0;
   
   if( !( edgeList.isEmpty() ) )
       newid = edgIter.LastP()->getID() + 1;
   else newid = 0;
   tempEdge1.setID( newid );                     //set edge1 ID
   newid++;
   tempEdge2.setID( newid );                     //set edge2 ID
   tempEdge1.setFlowAllowed( flowflag );         //set edge1 FLOWALLOWED
   tempEdge2.setFlowAllowed( flowflag );         //set edge2 FLOWALLOWED
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
     //add pointers to the new edges to nodes' spokeLists:
   spokIter.Reset( node2->getSpokeListNC() );
     //cout << "node " << node2->getID() << ": ";
   if( node2->getSpokeListNC().isEmpty() )
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " in otherwise empty list" << endl;
      node2->insertFrontSpokeList( le);
      node2->getSpokeListNC().makeCircular();
   }
   else if( spokIter.ReportNextP() == spokIter.DatPtr() )
   {
      node2->insertFrontSpokeList( le);
        //node2->getSpokeListNC().makeCircular();
   }
   else
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " w/ reference to node " << node3->getID() << endl;
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
      node2->getSpokeListNC().insertAtNext( le,
                                            spokIter.NodePtr() ); //put edge2 in SPOKELIST
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
   }
   else if( spokIter.ReportNextP() == spokIter.DatPtr() )
   {
      node1->insertFrontSpokeList( le );
        //node2->getSpokeListNC().makeCircular();
   }
   else
   {
        //cout << "place spoke " << edgIter.DatRef().getID()
        //   << " w/ reference to node " << node3->getID() << endl;
      for( ce = spokIter.FirstP();
           ce->getDestinationPtr() != node3 && !( spokIter.AtEnd() );
           ce = spokIter.NextP() );
      assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
      node1->getSpokeListNC().insertAtPrev( le,
                                            spokIter.NodePtr() );//put edge1 in SPOKELIST
   }
   
   nedges+=2;
   //Xcout << "2 edges added, edge list NOW:" << endl;
   //reset edge id's
   for( ce = edgIter.FirstP(), i = 0; !( edgIter.AtEnd() ); ce = edgIter.NextP(), i++ )
   {
      ce->setID( i );
      /*cout << "    Edg " << i << " (" << ce->getOriginPtr()->getID() << "->"
           << ce->getDestinationPtr()->getID() << ")\n";*/
   }
   return 1;
}

   
template< class tSubNode >
int tGrid< tSubNode >::
AddEdgeAndMakeTriangle( tPtrList< tSubNode > &nbrList,
                        tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   //cout << "AddEdgeAndMakeTriangle" << endl;
     //cout << "aemt 0\n" << flush;
   int flowflag = 1;
     //cout << "aemt 0.1\n" << flush;
   int i, j, newid;
   tSubNode *cn, *cnn, *cnnn;
   cn = nbrIter.DatPtr();
   cnn = nbrIter.NextP();
   cnnn = nbrIter.ReportNextP();
   nbrIter.Prev();
   tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
       p2( cnnn->get2DCoords() );
   if( !PointsCCW( p0, p1, p2 ) )
   {
       cout << "in AEAMT nodes not CCW: " << cn->getID() << ", "
            << cnn->getID() << ", " << cnnn->getID() << endl;
       cn->TellAll();
       cnn->TellAll();
       cnnn->TellAll();
   }
   
     //cout << "aemt 1/2\n" << flush;
   tEdge tempEdge1, tempEdge2, *ce, *le;
   tGridListIter< tEdge > edgIter( edgeList );
   tTriangle tempTri;
   tTriangle *nbrtriPtr, *ct;
   tListIter< tTriangle > triIter( triList );
   tPtrListIter< tEdge > spokIter;

   /*Xcout << "In AEMT, the three pts are:\n";
   cn->TellAll();
   cnn->TellAll();
   cnnn->TellAll();*/
   
        //DumpEdges:
//    cout << "edges:" << endl;
//    for( ce = edgIter.FirstP(); !( edgIter.AtEnd() ); ce = edgIter.NextP() )
//    {
//       cout << ce->getID() << " from " << ce->getOriginPtrNC()->getID()
//            << " to " << ce->getDestinationPtrNC()->getID() << endl;
//    }

   // Deal with new edges and new triangle:
   // Here, tempEdge1 and tempEdge2 are the complementary edges that connect
   // points cn and cnnn, which correspond to points 0 and 2 in the new
   // triangle. We set up the 3 vertex ptrs in the new triangle and set the
   // endpoints of the new edges.
   tempTri.setPPtr(0, cn );                      //set tri VERTEX ptr 0
   tempTri.setPPtr(1, cnn );                     //set tri VERTEX ptr 1
   tempTri.setPPtr(2, cnnn );                    //set tri VERTEX ptr 2
   tempEdge1.setOriginPtr( cn );                 //set edge1 ORG
   tempEdge1.setDestinationPtr( cnnn );          //set edge1 DEST
   tempEdge2.setOriginPtr( cnnn );               //set edge2 ORG
   tempEdge2.setDestinationPtr( cn );            //set edge2 DEST

   // get new unique ID's for the new edge pair
   le = edgIter.LastP();
   newid = le->getID() + 1;
   tempEdge2.setID( newid );                     //set edge1 ID
   newid++;
   tempEdge1.setID( newid );                     //set edge2 ID

   // set the boundary status of the edge pair: flow is allowed unless
   // one of the endpoints is a closed boundary, or both are open boundaries.
   if( cn->getBoundaryFlag() == kClosedBoundary ||
       cnnn->getBoundaryFlag() == kClosedBoundary ||
       (cn->getBoundaryFlag() == kOpenBoundary &&
        cnnn->getBoundaryFlag() == kOpenBoundary ) )
       flowflag = 0;
   tempEdge1.setFlowAllowed( flowflag );         //set edge1 FLOWALLOWED
   tempEdge2.setFlowAllowed( flowflag );         //set edge2 FLOWALLOWED
     //cout << "aemt 2\n" << flush;

   // The following block modified by gt 3/98 to change triangle edge-ptr
   // definition from "left-hand" to "right-hand". 
   
   // Add the new edges either to the back of the edge list (if they're
   // boundary edges) or to the back of the active portion of the edge list
   // (if not, ie if flow is allowed)
   if( flowflag == 1 )
   {
      edgeList.insertAtActiveBack( tempEdge2 );  //put edge2 active in list
      edgeList.insertAtActiveBack( tempEdge1 );  //put edge1 active in list
      le = edgIter.LastActiveP();                //set edgIter to lastactive
   }
   else
   {
      edgeList.insertAtBack( tempEdge2 );        //put edge2 in list
      edgeList.insertAtBack( tempEdge1 );        //put edge1 in list
      le = edgIter.LastP();                      //set edgIter to last
   }

   // le now points to edge1, P0->P2. Assign this edge as ePtr(0): the
   // clw-oriented edge that begins at vertex zero.
   tempTri.setEPtr(0, le );                      //set tri EDGE ptr 0
   
   /*Xcout << "In AEMT, adding edges: last = " << le->getID()
        << " (" << le->getOriginPtr()->getID() << "->"
        << le->getDestinationPtr()->getID() << "),\n";
   le->TellCoords();
   le = edgIter.PrevP();
   cout << "next to last = " << le->getID()
        << " (" << le->getOriginPtr()->getID() << "->"
        << le->getDestinationPtr()->getID() << ")\n";
   le->TellCoords();
   le = edgIter.NextP();*/
   
   // Now we add pointers to the new edges in the spokeLists of nodes 0 (cn)
   // and 2 (cnnn), and set the triangle's edge pointers ePtr(1) and ePtr(2)
   // to the clockwise-oriented edges P1->P0 and P2->P1.
   //   First add edge1 (P0->P2) to the spoke list for node cn (P0). Start by
   // finding the edge that runs cn->cnn (P0->P1), then insert tempEdge1 AFTER
   // (ccw from) it on the list
   spokIter.Reset( cn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   
   assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
   cn->getSpokeListNC().insertAtNext( le, spokIter.NodePtr() ); //put edge1 in SPOKELIST

   //   Next, find the spoke of cnnn (P2) whose destination is cnn (P1) (this
   // is the edge P2->P1 in the new triangle), and insert edge2
   // BEFORE this (ie clockwise of it) in the spoke list. Because le now points
   // to edge1, and edge2 comes just before edge1, we start by "backing up"
   // so that le points to edge2
     //cout << "aemt 3\n" << flush;
   le = edgIter.PrevP();  // point le to edge2 (P2->P0)
   spokIter.Reset( cnnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
   cnnn->getSpokeListNC().insertAtPrev( le, spokIter.NodePtr() );//put edge2 in SPOKELIST

   // ce now points to the edge P2->P1 (cnnn->cnn), which is also the
   // clw-oriented edge #2 for the triangle, so assign it as such
   tempTri.setEPtr(2, ce );                      //set tri EDGE ptr 2

   // Now we add the last remaining clw-oriented edge, cnn->cn (P1->P0),
   // to the triangle
   spokIter.Reset( cnn->getSpokeListNC() );
     //cout << "aemt 4\n" << flush;
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   tempTri.setEPtr(1, ce );                      //set tri EDGE ptr 1

   /* original code, unmodified follows
   if( flowflag == 1 )
   {
      edgeList.insertAtActiveBack( tempEdge1 );  //put edge1 active in list
      edgeList.insertAtActiveBack( tempEdge2 );  //put edge2 active in list
      le = edgIter.LastActiveP();                //set edgIter to lastactive
   }
   else
   {
      edgeList.insertAtBack( tempEdge1 );        //put edge1 in list
      edgeList.insertAtBack( tempEdge2 );        //put edge2 in list
      le = edgIter.LastP();                      //set edgIter to last
   }
   
   // Assign tempEdge2 as the ...
   tempTri.setEPtr(2, le );                      //set tri EDGE ptr 2
   
   //add pointers to the new edges to nodes' spokeLists:
   // First, find the spoke of cnnn whose destination is cnn (this is the
   // edge from pt 2 to pt 1 of the new triangle), and insert tempEdge2
   // BEFORE this (ie clockwise of it) in the spoke list
     //cout << "aemt 3\n" << flush;
   spokIter.Reset( cnnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
   cnnn->getSpokeListNC().insertAtPrev( le, spokIter.NodePtr() );//put edge2 in SPOKELIST
   le = edgIter.PrevP();

   // Now add tempEdge1 to the spoke list for node cn. Start by finding the
   // edge that runs cn to cnn, then insert tempEdge1 AFTER (counter-clw)
   // it on the list
   spokIter.Reset( cn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );  //make sure we found the right spoke
   cn->getSpokeListNC().insertAtNext( le, spokIter.NodePtr() ); //put edge1 in SPOKELIST

   // ce now contains the edge cn->cnn, which is also the ccw-oriented edge
   // #0 for the triangle, so assign it as such
   tempTri.setEPtr(0, ce );                      //set tri EDGE ptr 0

   // Now we add the last remaining ccw-oriented edge, cnn->cnnn,
   // to the triangle
   spokIter.Reset( cnn->getSpokeListNC() );
     //cout << "aemt 4\n" << flush;
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   tempTri.setEPtr(1, ce );                      //set tri EDGE ptr 1
*/
   
   // Give the new triangle a unique ID, initialize its tPtr's to NULL,
   // and insert it at the back of the triangle list
   tempTri.setID( ntri );                        //set tri ID
   for( i=0; i<3; i++ ) tempTri.setTPtr( i, 0 ); //initialize tri TRI ptrs to nul
   triList.insertAtBack( tempTri );              //put tri in list

   // The next task is to assign pointers to neighboring triangles. We start
   // by retrieving a ptr to the new triangle from the back of the list.
   ct = triIter.LastP();                         //set triIter to last

   /*cout <<"In AEMT, adding triangle:\n";
   ct->TellAll();*/

   // TODO: this op could be encapsulated in a tTriangle mbr function that
   // seeks a neighbor.
   // Find neighbor triangle 0. This is the triangle that includes edge
   // P1->P2 (cnn->cnnn) as one of its 3 clw-oriented edges. If such a
   // triangle exists, point to it as nbr tri #0, and also ask it to point
   // back to the new triangle as one of its neighbors.
   spokIter.Reset( cnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );
     //cout << "aemt 5\n" << flush;
   nbrtriPtr = TriWithEdgePtr( ce );
   ct->setTPtr( 0, nbrtriPtr );                   //set tri TRI ptr 0
   if( nbrtriPtr != 0 )
   {
      // Tell the nbr tri to point back to us. If the edge we share is i,
      // then we are nbri (i+1)%3
      for( i=0; nbrtriPtr->ePtr(i) != ce && i<3; i++ );
      assert( i<3 );
      nbrtriPtr->setTPtr( (i+1)%3, ct );         //set NBR TRI ptr to tri
   }

   // Neighbor triangle 2 is the one that points to edge P0-P1 (cn->cnn)
   spokIter.Reset( cn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );
   nbrtriPtr = TriWithEdgePtr( ce );
   ct->setTPtr( 2, nbrtriPtr );                   //set tri TRI ptr 2
     //cout << "aemt 6\n" << flush;
   if( nbrtriPtr != 0 )
   {
      for( i=0; nbrtriPtr->ePtr(i) != ce && i<3; i++ );
      assert( i<3 );
      nbrtriPtr->setTPtr( (i+1)%3, ct );         //set NBR TRI ptr to tri
   }

   /* original unmodified second block follows:
   // The next task is to assign pointers to neighboring triangles...
   // (I don't understand this: why wd any other triangle share P2->P0 etc
   // if it's lefthand-oriented?)
   ct = triIter.LastP();                         //set triIter to last
   spokIter.Reset( cnnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cnn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );
     //cout << "aemt 5\n" << flush;
   nbrtriPtr = TriWithEdgePtr( ce );
   ct->setTPtr(0, nbrtriPtr );                   //set tri TRI ptr 0
   if( nbrtriPtr != 0 )
   {
      for( i=0; nbrtriPtr->ePtr(i) != ce && i<3; i++ );
      assert( i<3 );
      nbrtriPtr->setTPtr( (i+2)%3, ct );         //set NBR TRI ptr to tri
   }
   spokIter.Reset( cnn->getSpokeListNC() );
   for( ce = spokIter.FirstP();
        ce->getDestinationPtr() != cn && !( spokIter.AtEnd() );
        ce = spokIter.NextP() );
   assert( !( spokIter.AtEnd() ) );
   nbrtriPtr = TriWithEdgePtr( ce );
   ct->setTPtr(2, nbrtriPtr );                   //set tri TRI ptr 2
     //cout << "aemt 6\n" << flush;
   if( nbrtriPtr != 0 )
   {
      for( i=0; nbrtriPtr->ePtr(i) != ce && i<3; i++ );
      assert( i<3 );
      nbrtriPtr->setTPtr( (i+2)%3, ct );         //set NBR TRI ptr to tri
   }*/
   
   ntri++;                          //increment tGrid::ntri by one
   nedges+=2;                       //increment tGrid::nedges by two
   
     //cout << "aemt 7\n" << flush;
   //reset edge id's (NB: why needed?)
   for( ce = edgIter.FirstP(), i=0; !( edgIter.AtEnd() );
        ce = edgIter.NextP(), i++ )
       ce->setID( i );
     //cout << "aemt 8\n" << flush;
   
   //reset triangle id's (NB: why needed?)
   for( ct = triIter.FirstP(), i=0; !( triIter.AtEnd() );
        ct = triIter.NextP(), i++ )
   {
      ct->setID( i );
      assert( ( ( ct->tPtr(0) != ct->tPtr(1) && ct->tPtr(0) != ct->tPtr(2) ) ||
                ct->tPtr(0) == 0 ) &&
              ( ( ct->tPtr(1) != ct->tPtr(0) && ct->tPtr(1) != ct->tPtr(2) ) ||
                ct->tPtr(1) == 0 ) &&
              ( ( ct->tPtr(2) != ct->tPtr(0) && ct->tPtr(2) != ct->tPtr(1) ) ||
                ct->tPtr(2) == 0 ) );
   }

     //cout << "  finished" << endl << flush;
   return 1;
}

// TODO: consolidate w/ AEMT to have AEMT call this as a helper fn
template< class tSubNode >
int tGrid< tSubNode >::
MakeTriangle( tPtrList< tSubNode > &nbrList,
              tPtrListIter< tSubNode > &nbrIter )
{
   assert( (&nbrList != 0) && (&nbrIter != 0) );
   assert( nbrList.getSize() == 3 );
   //cout << "MakeTriangle" << endl;
   int i, j;
   int newid;                          // ID of new triangle
   tTriangle tempTri;
   tTriangle *nbrtriPtr;
   tSubNode *cn, *cnn, *cnnn;
   tEdge *ce, *dce;
   tTriangle *ct;
   tListIter< tTriangle > triIter( triList );
   tGridListIter< tEdge > edgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   assert( nbrList.getSize() == 3 );
   //Xcn = nbrIter.DatPtr(); Below is bug fix:
   cn = nbrIter.FirstP();      // cn, cnn, and cnnn are the 3 nodes in the tri
   cnn = nbrIter.NextP();
   cnnn = nbrIter.NextP();
   nbrIter.Next();
   tArray< double > p0( cn->get2DCoords() ), p1( cnn->get2DCoords() ),
       p2( cnnn->get2DCoords() );
   if( !PointsCCW( p0, p1, p2 ) )
   {
      cout << "in MT nodes not CCW: " << cn->getID() << ", "
           << cnn->getID() << ", " << cnnn->getID();
      if( cn->Meanders() ) p0 = cn->getNew2DCoords();
      else p0 = cn->get2DCoords();
      if( cnn->Meanders() ) p1 = cnn->getNew2DCoords();
      else p1 = cnn->get2DCoords();
      if( cnnn->Meanders() ) p2 = cnnn->getNew2DCoords();
      else p2 = cnnn->get2DCoords();
      if( !PointsCCW( p0, p1, p2 ) )
          cout << "; nor are new coords CCW ";
      cout << endl;
   }

   /*cout << "In MT, the 3 nodes are: " << cn->getID() << " " << cnn->getID()
        << " " << cnnn->getID() << endl;*/
   
     //for debugging:
     //DumpEdges();
     //DumpSpokes:
   
   // set the ID for the new triangle based on the ID of the last triangle
   // on the list plus one, or if there are no triangles on the list yet
   // (which happens when we're creating an initial "supertriangle" as in
   // MakeGridFromPoints), set the ID to zero.
   ct = triIter.LastP();
   if( ct ) newid = ct->getID() + 1;
   else newid = 0;
   tempTri.setID( newid );

   // set edge and vertex ptrs & add to triList: We go through each point,
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

      /*cout << "SEEKing edg from " << ce->getOriginPtr()->getID()
           << " to " << cnn->getID() << " and found it in edg " << ce->getID()
           << endl;
      ce->TellCoords();*/
      
      // 4 debug
      /*if( ( spokIter.AtEnd() ) )
      {
         cout << "dest node " << cn->getID() << " not found for node "
              << ce->getOriginPtrNC()->getID() << endl;
         DumpNodes();
         ReportFatalError( "failed: !( spokIter.AtEnd() )" );
      }*/
      
      assert( !( spokIter.AtEnd() ) );

      // Assign edge p(i)->p(i+2) as the triangle's clockwise edge #i
      // (eg, ePtr(0) is the edge that connects p0->p2, etc)
      tempTri.setEPtr( i, ce );      //set tri EDGE ptr i
      tempTri.setTPtr( i, 0 );       //initialize tri TRI ptrs to nul
   }
   triList.insertAtBack( tempTri );       //put tri in list
   ct = triIter.LastP();                  //set triIter to last
   assert( cn == ct->pPtr(0) );           //make sure we're where we
                                          //think we are

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
   dce = 0;
   nbrtriPtr = 0;
   cn = nbrIter.FirstP();
   //cout << "starting w/ node " << cn->getID();
   for( j=0; j<3; j++ )
   {
      // get spokelist for p(j) and advance to p(j+1)
      spokIter.Reset( cn->getSpokeListNC() );
      cn = nbrIter.NextP();               //step forward once in nbrList
      if( j>0 ) dce = ce;

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
             cout << "something FUNNY going on";
         else cout << "tri not CCW: " << nbrtriPtr->getID() << endl;
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
   for( ct = triIter.FirstP(), i=0; !( triIter.AtEnd() );
        ct = triIter.NextP(), i++ )
   {
      ct->setID( i );
   }
   
   return 1;
}


/**************************************************************************\
**
**   tGrid::AddNode ( tSubNode nodeRef& )
**
**   Adds a new node with the properties of nodRef to the mesh.
**
**   Calls: tGrid::LocateTriangle, tGrid::DeleteTriangle, tGrid::AddEdge,
**            tGrid::AddEdgeAndMakeTriangle, tGrid::MakeTriangle,
**            tGrid::CheckForFlip; various member functions of tNode,
**            tGridList, tGridListIter, tPtrList, etc. Also tLNode
**            functions (TODO: this needs to be removed somehow),
**            and temporarily, tGrid::UpdateMesh
**   Parameters: nodeRef -- reference to node to be added (really,
**                          duplicated)
**   Returns:  (always TRUE: TODO make void return type)
**   Assumes:
**   Modifications:
**        - 4/98: node is no longer assumed to be a non-boundary (GT)
**        - 7/98: changed return type from int (0 or 1) to ptr to
**                the new node (GT)
**
\**************************************************************************/
#define kLargeNumber 1000000000
template< class tSubNode >
tSubNode * tGrid< tSubNode >::
AddNode( tSubNode &nodeRef )
{
   int i, j, k, ctr;
   tTriangle *tri;
   tSubNode *cn;
   tEdge * tedg1;
   tEdge * tedg3;
   tArray< double > xyz( nodeRef.get3DCoords() );
   tGridListIter< tSubNode > nodIter( nodeList );
   assert( &nodeRef != 0 );

   //cout << "AddNode at " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;

   //cout << "locate tri" << endl << flush;
   tri = LocateTriangle( xyz[0], xyz[1] );
   assert( tri != 0 );  //if( tri == 0 ) return 0;
   
   // Assign ID to the new node and insert it at the back of either the active
   // portion of the node list (if it's not a boundary) or the boundary
   // portion (if it is)
   int newid = nodIter.LastP()->getID() + 1;
   nodeRef.setID( newid );
   if( nodeRef.getBoundaryFlag()==kNonBoundary )
       nodeList.insertAtActiveBack( nodeRef );
   else
       nodeList.insertAtBack( nodeRef );
   assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
   
   // Retrieve a pointer to the new node and flush its spoke list
   if( nodeRef.getBoundaryFlag()==kNonBoundary )
       cn = nodIter.LastActiveP();
   else
       cn = nodIter.LastP();
   assert( cn!=0 );
   cn->getSpokeListNC().Flush();
   
   //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = (tSubNode *) tri->pPtr(i);
      bndyList.insertAtBack( tmpPtr );
   }
   bndyList.makeCircular();


   // Delete the triangle in which the node falls
   //Xcout << "deleting tri in which new node falls\n";
   i = DeleteTriangle( tri );
   assert( i != 0 );  //if ( !DeleteTriangle( tri ) ) return 0;


   //make 3 new triangles
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
         {
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
   for( cn = nodIter.FirstP(), i=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), i++ )
   {
      cn->setID( i );
   }
   node2->makeCCWEdges();

   UpdateMesh();
   //cout << "AddNode finished" << endl; 

   return node2;  // Return ptr to new node
}


/**************************************************************************\
**
**   tGrid::AddNode ( tSubNode nodeRef&, int dum )
**
**   Adds a new node with the properties of nodRef to the mesh.
**   Same as AddNode above, it just takes a dummy integer to indicate
**   that layer interpolation should be done.
**   This is not a good solution, this function should only be
**   temporary until a better way is thought of to indicate that
**   layer interpolation is necessary.
**
**   Calls: tGrid::LocateTriangle, tGrid::DeleteTriangle, tGrid::AddEdge,
**            tGrid::AddEdgeAndMakeTriangle, tGrid::MakeTriangle,
**            tGrid::CheckForFlip; various member functions of tNode,
**            tGridList, tGridListIter, tPtrList, etc. Also tLNode
**            functions (TODO: this needs to be removed somehow),
**            and temporarily, tGrid::UpdateMesh
**   Parameters: nodeRef -- reference to node to be added (really,
**                          duplicated)
**   Returns:  (always TRUE: TODO make void return type)
**   Assumes:
**   Modifications:
**        - 4/98: node is no longer assumed to be a non-boundary (GT)
**        - 7/98: changed return type from int (0 or 1) to ptr to
**                the new node (GT)
**
\**************************************************************************/
#define kLargeNumber 1000000000
template< class tSubNode >
tSubNode * tGrid< tSubNode >::
AddNode( tSubNode &nodeRef, int dum )
{
   int i, j, k, ctr;
   tTriangle *tri;
   tSubNode *cn;
   tArray< double > xyz( nodeRef.get3DCoords() );
   tGridListIter< tSubNode > nodIter( nodeList );
   assert( &nodeRef != 0 );

   //cout << "AddNode at " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;

   //cout << "locate tri" << endl << flush;
   tri = LocateTriangle( xyz[0], xyz[1] );
   assert( tri != 0 );  //if( tri == 0 ) return 0;

   // Assign ID to the new node and insert it at the back of either the active
   // portion of the node list (if it's not a boundary) or the boundary
   // portion (if it is)
   int newid = nodIter.LastP()->getID() + 1;
   nodeRef.setID( newid );
   if( nodeRef.getBoundaryFlag()==kNonBoundary )
       nodeList.insertAtActiveBack( nodeRef );
   else
       nodeList.insertAtBack( nodeRef );
   assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
   
   // Retrieve a pointer to the new node and flush its spoke list
   if( nodeRef.getBoundaryFlag()==kNonBoundary )
       cn = nodIter.LastActiveP();
   else
       cn = nodIter.LastP();
   assert( cn!=0 );
   cn->getSpokeListNC().Flush();
   
     //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = (tSubNode *) tri->pPtr(i);
      bndyList.insertAtBack( tmpPtr );
   }
   bndyList.makeCircular();
   
   // Delete the triangle in which the node falls
   //Xcout << "deleting tri in which new node falls\n";
   i = DeleteTriangle( tri );
   assert( i != 0 );  //if ( !DeleteTriangle( tri ) ) return 0;

   //make 3 new triangles
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
   
   //put 3 resulting triangles in ptr list
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
         {
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
   for( cn = nodIter.FirstP(), i=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), i++ )
   {
      cn->setID( i );
   }
   //node2->makeCCWEdges();
   UpdateMesh();
   //cout << "AddNode finished" << endl;
   return node2;  // Return ptr to new node
}


/*   AddNodeAt: add a node with referenced coordinates to mesh   */
//TODO: this fn duplicates functionality of AddNode; just assign coords
// to a dummy new node and call AddNode
template< class tSubNode >
tSubNode *tGrid< tSubNode >::
AddNodeAt( tArray< double > &xyz )
{
   assert( &xyz != 0 );
   cout << "AddNodeAt " << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << endl;
   tTriangle *tri;
   //cout << "locate tri" << endl << flush;
   if( xyz.getSize() == 3 ) tri = LocateTriangle( xyz[0], xyz[1] );
   else tri = LocateNewTriangle( xyz[0], xyz[1] );
   if( tri == 0 ) return 0;
   int i, j, k, ctr;
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode tempNode, *cn;
   tempNode.set3DCoords( xyz[0], xyz[1], xyz[2]  );
   if( layerflag )
       tempNode.LayerInterpolation( tri, xyz[0], xyz[1] );
   if( xyz.getSize() != 3 ) tempNode.setNew2DCoords( xyz[0], xyz[1] );
   tempNode.setBoundaryFlag( 0 );

   // Assign ID to the new node and insert it at the back of the active
   // portion of the node list (NOTE: node is assumed NOT to be a boundary)
   int newid = nodIter.LastP()->getID() + 1;
   tempNode.setID( newid );
   nodeList.insertAtActiveBack( tempNode );
   assert( nodeList.getSize() == nnodes + 1 );
   nnodes++;
     //make ptr list of triangle's vertices:
   tPtrList< tSubNode > bndyList;
   tSubNode *tmpPtr;
   for( i=0; i<3; i++ )
   {
      tmpPtr = (tSubNode *) tri->pPtr(i);
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
      /*cout << "   in triangle w/ vtcs. at " << p3[0] << " " << p3[1] << "; "
           << p1[0] << " " << p1[1] << "; " << p4[0] << " " << p4[1] << endl;*/
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
   for( cn = nodIter.FirstP(), i=0; !( nodIter.AtEnd() ); cn = nodIter.NextP(), i++ )
   {
      cn->setID( i );
   }
   //node2->makeCCWEdges();
   UpdateMesh();
   //cout << "AddNodeAt finished, " << nnodes << endl;
   return node2;
}
#undef kLargeNumber

template <class tSubNode>
tGridList<tEdge> * tGrid<tSubNode>::
getEdgeList() {return &edgeList;}

template <class tSubNode>
tGridList<tSubNode> * tGrid<tSubNode>::
getNodeList() {return &nodeList;}

template <class tSubNode>
tList< tTriangle > * tGrid<tSubNode>::
getTriList() {return &triList;}

template< class tSubNode >
tEdge *tGrid< tSubNode >::
getEdgeComplement( tEdge *edge )
{
   tGridListIter< tEdge > edgIter( edgeList );
   int edgid = edge->getID();
   assert( edgIter.Get( edgid ) );
   edgIter.Get( edgid );
   if( edgid%2 == 0 ) return edgIter.GetP( edgid + 1 );
   if( edgid%2 == 1 ) return edgIter.GetP( edgid - 1 );
}

template <class tSubNode>
void tGrid<tSubNode>::
UpdateMesh()
{
   //cout << "UpdateMesh()" << endl;
   
   //tListIter<tTriangle> tlist( triList );
   tGridListIter<tEdge> elist( edgeList );
   //tGridListIter<tSubNode> nlist( nodeList );
   tEdge * curedg = 0;
   double len;
   tTriangle * curtri;
   
   // Edge lengths
   curedg = elist.FirstP();
   do
   {
      len = curedg->CalcLength();
      //Xcout << "Edge " << curedg->getID() << " length: " << curedg->getLength() << endl;
      curedg = elist.NextP();
      assert( curedg > 0 ); // failure = complementary edges not consecutive
      curedg->setLength( len );
   } while( curedg=elist.NextP() );

   MakeCCWEdges();

   setVoronoiVertices();
   CalcVoronoiEdgeLengths();
   CalcVAreas();
   CheckMeshConsistency( 0 );

// Triangle areas
/*   for( tlist.First(); !tlist.AtEnd(); tlist.Next() )
   {
      curtri = tlist.DatPtr();
      curtri->length_sides();
      curtri->CalcArea();
      curtri = curtri->next;
   }
   */
   
   // Voronoi vertices
   //getVoronoiVertices();

   

   // Voronoi Areas


}


/*****************************************************************************\
**
**      CheckForFlip: checks whether edge between two triangles should be
**                    flipped; may either check, flip, and report, or just
**                    check and report.
**                    Checks whether the present angle or the possible angle
**                    is greater. Greater angle wins. Also uses flip variable
**                    to determine whether to use newx, newy, or x, y.
**
**      Data members updated: Grid
**      Called by: 
**      Calls:  
**        
**      Created: 8/28/97 SL
**      Modified: 12/16/97 SL                                               
**                                                              
**
\*****************************************************************************/
template< class tSubNode >
int tGrid< tSubNode >::
CheckForFlip( tTriangle * tri, int nv, int flip )
{
   if( tri == 0 )
   {
      cout << "CheckForFlip: tri == 0" << endl;
      return 0;
   }
   assert( nv < 3 );
   //cout << "THIS IS CheckForFlip(...) " << tri->getID() << endl;
   tSubNode *node0, *node1, *node2, *node3;
   node0 = ( tSubNode * ) tri->pPtr(nv);
   //cout<<"node0 id "<<node0->getID()<<endl;
   node1 = ( tSubNode * ) tri->pPtr((nv+1)%3);
   //cout<<"node1 id "<<node1->getID()<<endl;
   node2 = ( tSubNode * ) tri->pPtr((nv+2)%3);
   //cout<<"node2 id "<<node2->getID()<<endl;
   tTriangle *triop = tri->tPtr(nv);
   int nvop = triop->nVOp( tri );
   node3 = ( tSubNode * ) triop->pPtr( nvop );
   tArray< double > ptest( node3->get2DCoords() ), p0( node0->get2DCoords() ),
       p1( node1->get2DCoords() ), p2( node2->get2DCoords() );
   if( !flip )
   {
      if( node0->Meanders() ) p0 = node0->getNew2DCoords();
      if( node1->Meanders() ) p1 = node1->getNew2DCoords();
      if( node2->Meanders() ) p2 = node2->getNew2DCoords();
      if( node3->Meanders() ) ptest = node3->getNew2DCoords();
   }
   if( TriPasses( ptest, p0, p1, p2 ) ) return 0;
   if( flip )                     //and make sure there isn't already an edge?
   {
      if( !PointsCCW( p0, p1, ptest ) || !PointsCCW( p0, ptest, p2 ) )return 0;
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
   Note on notation in flip edge

                d
               /|\
       tri->  / | \ <-triop
             /  |  \
            a   |   c
             \  |  /
              \ | /
               \|/
                b
        Edge bd will be removed
        and an edge ac will be made.
        nbrList contains the points a, b, c, d     
\******************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
FlipEdge( tTriangle * tri, tTriangle * triop ,int nv, int nvop )
{
   //cout << "FlipEdge(...)..." << endl;
   tSubNode *cn;
   tPtrList< tSubNode > nbrList;
   //DumpTriangles();
   nbrList.insertAtBack( (tSubNode *) tri->pPtr(nv) );
   nbrList.insertAtBack( (tSubNode *) tri->pPtr((nv+1)%3) );
   nbrList.insertAtBack( (tSubNode *) triop->pPtr( nvop ) );
   nbrList.insertAtBack( (tSubNode *) tri->pPtr((nv+2)%3) );
   nbrList.makeCircular();
   //cout << "calling deleteedge from flipedge\n";
   //XDeleteEdge( tri->ePtr( (nv+1)%3 ) );
   DeleteEdge( tri->ePtr( (nv+2)%3 ) );  // Changed for right-hand data struc
   tPtrListIter< tSubNode > nbrIter( nbrList );
   AddEdgeAndMakeTriangle( nbrList, nbrIter );
   nbrIter.First();
   nbrList.removeNext( cn, nbrIter.NodePtr() );
   MakeTriangle( nbrList, nbrIter );
   //cout << "finished" << endl;
}


/*****************************************************************************\
**
**      CheckLocallyDelaunay : updates the triangulation after moving
**             some points.
**             only uses x and y values, which have already been updated in
**             MoveNodes (frmr PreApply).
**             PREAPPLY SHOULD BE CALLED BEFORE THIS FUNCTION IS CALLED
**      Data members updated: Grid
**      Called by: MoveNodes
**      Calls:  
**        
\*****************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
CheckLocallyDelaunay()
{
   //cout << "CheckLocallyDelaunay()" << endl;
   tTriangle *at, * trop[3];
   tPtrList< tTriangle > triPtrList;
   tPtrListIter< tTriangle > triPtrIter( triPtrList );
   tListIter< tTriangle > triIter( triList );
   int i, change, flipped;
   int id0, id1, id2;
   tArray< int > npop(3);
   tSubNode *nodPtr;
   
   int flip = 1;
   //first find a way to go through all the points an the line
   //put each triangle into the stack
     //flipped = TRUE;
     /*do
   {*/
     //flipped = FALSE;
   for( at = triIter.FirstP(); !( triIter.AtEnd() ); at = triIter.NextP() )
   {
      change = FALSE;
      for( i = 0; i < 3; i++ )
      {
         nodPtr = ( tSubNode * ) at->pPtr(i);
         if( nodPtr->Meanders() ) change = TRUE;
      }
      if( change ) triPtrList.insertAtBack( at );
   }

   //check list for flips; if flip, put new triangles at end of list
   tPtrListIter< tTriangle > duptriPtrIter( triPtrList );
   tTriangle *tn, *tp;
   while( !( triPtrList.isEmpty() ) )
   {
      at = triPtrIter.FirstP();
      for( i=0; i<3; i++ )
      {
         if( at->tPtr(i) != 0 )
         {
            tp = at->tPtr(i);
            for( tn = duptriPtrIter.FirstP();
                 duptriPtrIter.ReportNextP() != tp &&
                     !( duptriPtrIter.AtEnd() );
                 tn = duptriPtrIter.NextP() );
            tn = 0;
            if( !( duptriPtrIter.AtEnd() ) )
            {
               tn = duptriPtrIter.ReportNextP();
            }
            if( at->tPtr(0) != 0 ) id0 = at->tPtr(0)->getID();
            else id0 = -1;
            if( at->tPtr(1) != 0 ) id1 = at->tPtr(1)->getID();
            else id1 = -1;
            if( at->tPtr(2) != 0 ) id2 = at->tPtr(2)->getID();
            else id2 = -1;
            /*cout << "check tri " << at->getID() << " with nbrs "
                 << id0 << ", " << id1
                 << ", and " << id2;*/
               
            if( tp->tPtr(0) != 0 ) id0 = tp->tPtr(0)->getID();
            else id0 = -1;
            if( tp->tPtr(1) != 0 ) id1 = tp->tPtr(1)->getID();
            else id1 = -1;
            if( tp->tPtr(2) != 0 ) id2 = tp->tPtr(2)->getID();
            else id2 = -1;
            /*cout << " against tri " << tp->getID() << " with nbrs "
                 << id0 << ", " << id1
                 << ", and " << id2 << endl;*/
            //cout << "call cff from cld\n";
            if( CheckForFlip( at, i, flip ) )
            {
               //cout << "flipped tri's, got tri ";
               if( tn != 0 )
                   triPtrList.removeNext( tn, duptriPtrIter.NodePtr() );
               tn = triIter.LastP();
               if( tn->tPtr(0) != 0 ) id0 = tn->tPtr(0)->getID();
               else id0 = -1;
               if( tn->tPtr(1) != 0 ) id1 = tn->tPtr(1)->getID();
               else id1 = -1;
               if( tn->tPtr(2) != 0 ) id2 = tn->tPtr(2)->getID();
               else id2 = -1;
               /*cout << tn->getID() << " with nbrs "
                    << id0 << ", " << id1
                    << ", and " << id2;*/
               triPtrList.insertAtBack( tn );
               tn = triIter.PrevP();
               if( tn->tPtr(0) != 0 ) id0 = tn->tPtr(0)->getID();
               else id0 = -1;
               if( tn->tPtr(1) != 0 ) id1 = tn->tPtr(1)->getID();
               else id1 = -1;
               if( tn->tPtr(2) != 0 ) id2 = tn->tPtr(2)->getID();
               else id2 = -1;
               /*cout << " and tri " << tn->getID() << " with nbrs "
                    << id0 << ", " << id1
                    << ", and " << id2 << endl;*/
               triPtrList.insertAtBack( tn );
               break;
            }
         }
      }
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
**      IntersectsAnyEdge: returns the first edge in the list which intersects
**                         "edge" or NULL if "edge" intersects no other edges
**      Data members updated: Grid
**      Called by: 
**      Calls:  
**
\*****************************************************************************/
template< class tSubNode >
tEdge *tGrid< tSubNode >::
IntersectsAnyEdge( tEdge * edge )
{
   //cout << "IntersectsAnyEdge( tEdge * edge )..." << endl;
   int i;
   tEdge * ce;
   tGridListIter< tEdge > edgIter( edgeList );
   if( !edge )
   {
      cout<<"IntersectsAnyEdge: Warning: invalid edge"<<endl<<flush;
      return( NULL );
   }
     //cout<<"IAE: nedges "<<nedges<<endl<<flush;
     //cout << "call Intersect for edges " << edge->getID()
     //   << " from nodes " << edge->getOriginPtr()->getID()
     //   << " to " << edge->getDestinationPtr()->getID() << "; " << endl;
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
     /*if( i < nedges - 1 )
       cout<<"IntersectsAnyEdge: Warning: whole list not checked"<<endl<<flush;*/
   return( NULL );
}

/*****************************************************************************\
**
**      CheckTriEdgeIntersect():
**        We want to know if the moving point has passed beyond the polygon
**        defined by its spoke edges; if it has, then we will have edges
**        intersecting one another. In the case where the point has simply
**        passed into one of the 'opposite' triangles, then we can just do a
**        flip operation. In the other case, the remedial action is much more
**        complicated, so we just delete the point and add it again.
**
**  Modifications:
**   - minor change from i = AddNode to cn = AddNode to handle changed
**     return type (GT 7/98)
**
\*****************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
CheckTriEdgeIntersect()
{
   //cout << "CheckTriEdgeIntersect()..." << flush << endl;
     //DumpNodes();
   int i, j, nv, nvopp, id0, id1, id2;
   int flipped = TRUE;
   int crossed;
   tSubNode *subnodePtr, tempNode, newNode;  
   tEdge * newedg, * cedg, * ccedg, *fedg, *ce, *cex;
   tTriangle * ct, * ctop, *rmtri, *tri;
   tListIter< tTriangle > triIter( triList );
   tGridListIter< tEdge > edgIter( edgeList );
   tGridListIter< tSubNode > nodIter( nodeList );
   tGridListIter< tEdge > xedgIter( edgeList );
   tPtrListIter< tEdge > spokIter;
   tGridList< tSubNode > tmpNodeList;
   tGridListIter< tSubNode > tmpIter( tmpNodeList );
   tArray< double > p0, p1, p2, xy, xyz, xy1, xy2;
   tSubNode *cn, *vtx;
   tPtrList< tTriangle > triptrList;
   tPtrListNode< tTriangle > *tpListNode;
   tPtrListIter< tTriangle > tpIter( triptrList );
     //check for triangles with edges which intersect (an)other edge(s)
     //newedg = new tEdge;
   while( flipped )
   {
      flipped = FALSE;
      for( ct = triIter.FirstP(); !( triIter.AtEnd() ); ct = triIter.NextP() )
      {
         for( i=0; i<3; i++ )
         {
            cn = (tSubNode *) ct->pPtr(i);
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
                  subnodePtr = (tSubNode *) ct->pPtr(i);
                  subnodePtr->RevertToOldCoords();
               }
            }
            else
            {   
               crossed = FALSE;
               for( i=0; i<3; i++ )
               {
                  cn = (tSubNode *) ct->pPtr(i);
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
                              subnodePtr = (tSubNode *) ct->pPtr(i);
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
                                       for( tri = tpIter.FirstP();
                                            tpIter.ReportNextP() != rmtri &&
                                                !(tpIter.AtEnd());
                                            tri = tpIter.NextP() );
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
                                    subnodePtr = (tSubNode *) ct->pPtr(i);
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
**      MoveNodes (formerly PreApply) :
**             The function that deleted some points so that 
**             they won't be a problem in ApplyChanges
**               uses newx newy x y.
**             Separate call to ApplyChanges now unnecessary!!!
**      Data members updated: Grid
**      Called by: 
**      Calls:  ApplyChanges
**        
**      Created: SL
**                                                              
**
\*****************************************************************************/
template< class tSubNode >
void tGrid< tSubNode >::
MoveNodes()
{

   cout << "MoveNodes()..." << flush << endl;
   tSubNode * cn;  
   tGridListIter< tSubNode > nodIter( nodeList );
   //Before any edges and triangles are changed, layer interpolation
   //must be performed.
   if( layerflag ){   
      tTriangle *tri;
      tArray<double> newxy(2);
      for(cn=nodIter.FirstP(); nodIter.IsActive(); cn=nodIter.NextP()){
         newxy=cn->getNew2DCoords();
         if( (cn->getX()!=newxy[0]) && (cn->getY()!=newxy[1]) ){
            //Nic - there may be some issues if boundary nodes make up
            //the triangle.
            cout<<"a point will be moved in movenodes"<<endl;
            cout<<cn->getX()<<" "<<newxy[0]<<" "<<cn->getY()<<" "<<newxy[1]<<endl;
            tri = LocateTriangle( newxy[0], newxy[1] );
            cn->LayerInterpolation( tri, newxy[0], newxy[1] );
         }
      }
   }
   
   //check for triangles with edges which intersect (an)other edge(s)
   CheckTriEdgeIntersect(); //calls tLNode::UpdateCoords() for each node
   //resolve any remaining problems after points moved
   CheckLocallyDelaunay();
   UpdateMesh();
   CheckMeshConsistency();
   //cout << "MoveNodes() finished" << endl;
}


#ifndef NDEBUG
/*****************************************************************************\
**
**      DumpEdges(), DumpSpokes(), DumpTriangles(), DumpNodes(): debugging
**         routines which simply write out information pertaining to the grid;
**      DumpNodes() calls DumpSpokes for each node;
**      DumpSpokes() takes a pointer to a node as an argument.
**
**      Created: SL 1/98
**
\*****************************************************************************/
template<class tSubNode>
void tGrid<tSubNode>::
DumpEdges()
{
   tGridListIter< tEdge > edgIter( edgeList );
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
void tGrid<tSubNode>::
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
void tGrid<tSubNode>::
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
void tGrid<tSubNode>::
DumpNodes()
{
   tGridListIter< tSubNode > nodIter( nodeList );
   tSubNode *cn;
   cout << "nodes: " << endl;
   for( cn = nodIter.FirstP(); !(nodIter.AtEnd()); cn = nodIter.NextP() )
   {
      cout << " at " << cn->getX() << ", " << cn->getY() << ", " << cn->getZ() << "; ";
      DumpSpokes( cn );
   }
}
#endif
