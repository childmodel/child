/**************************************************************************\
**
**  tListInputData.cpp
**
**  Functions for classes tListInputData and tListIFStreams.
**
**  Modifications:
**   - changed .tri file format from points-edges-triangles to
**     points-triangles-edges, compatible with earlier format (gt 1/98)
**
**  $Id: tListInputData.cpp,v 1.3 1998-02-01 00:56:42 stlancas Exp $
\**************************************************************************/

#include "tListInputData.h"

/**************************************************************************\
**
**         Utilities for classes tListIFStreams &  tListInputData
**
\**************************************************************************/
/*tListIFStreams::tListIFStreams()                               //tListIFStreams
{
   cout << "enter no. nodes" << endl;
   cin >> nnodes;
   cout << "enter no. edges" << endl;
   cin >> nedges;
   cout << "enter no. triangles" << endl;
   cin >> ntri;
   cout << "tInputChannel()" << endl;
}*/
tListIFStreams::tListIFStreams( tInputFile &infile )
{
   char thestring[80], inname[80];
   infile.ReadItem( thestring, "INPUTDATAFILE" );
   strcpy( inname, thestring );
   strcat( inname, ".nodes" );
   nodeinfile.open(inname); /* Node input file pointer */
   assert( nodeinfile.good() );

   strcpy( inname, thestring );
   strcat( inname, ".edges" );
   edgeinfile.open(inname); /* Edge input file pointer */
   assert( edgeinfile.good() );

   strcpy( inname, thestring );
   strcat( inname, ".tri" );
   triinfile.open( inname ); //Triangle input file pointer
   assert( triinfile.good() );
   
   nodeinfile >> nnodes;
   edgeinfile >> nedges;
   triinfile >> ntri;
}

tListIFStreams::tListIFStreams( const char * argv = 0 )        //tListIFStreams
{
   if( argv == 0 )
   {
      cout << "enter no. nodes" << endl;
      cin >> nnodes;
      cout << "enter no. edges" << endl;
      cin >> nedges;
      cout << "enter no. triangles" << endl;
      cin >> ntri;
      //nodeinfile = 0;
   }
   else
   {
      char inname[80];
      strcpy( inname, argv );
      strcat( inname, ".nodes" );
      nodeinfile.open(inname); /* Node input file pointer */
      assert( nodeinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, ".edges" );
      edgeinfile.open(inname); /* Edge input file pointer */
      assert( edgeinfile.good() );
   
      strcpy( inname, argv );
      strcat( inname, ".tri" );
      triinfile.open( inname ); //Triangle input file pointer
      assert( triinfile.good() );
   
      nodeinfile >> nnodes;
      edgeinfile >> nedges;
      triinfile >> ntri;
   }
   cout << "tListIFStreams( argv )" << endl;
}

tListIFStreams::~tListIFStreams()                              //tListIFStreams
{   
   nodeinfile.close();
   edgeinfile.close();
   triinfile.close();
   cout << "    ~tListIFStreams()" << endl;
}

int tListIFStreams::getNNodes() const {return nnodes;}        //tListIFStreams

int tListIFStreams::getNEdges() const {return nedges;}        //tListIFStreams

int tListIFStreams::getNTri() const {return ntri;}            //tListIFStreams

ifstream &tListIFStreams::getNodeInFile()                     //tListIFStreams
{return nodeinfile;}

ifstream &tListIFStreams::getEdgeInFile()                     //tListIFStreams
{return edgeinfile;}

ifstream &tListIFStreams::getTriInFile()                      //tListIFStreams
{return triinfile;}

//default constructor
//template< class tSubNode >
//tListInputData< tSubNode >::tListInputData()                //tListInputData
//constructor
template< class tSubNode >
tListInputData< tSubNode >::
tListInputData( tListIFStreams &inchan )                             //tListInputData
        : x( inchan.getNNodes() ), y( inchan.getNNodes() ),
          z( inchan.getNNodes() ), edgid( inchan.getNNodes() ),
          boundflag( inchan.getNNodes() ),
          orgid( inchan.getNEdges() ), destid( inchan.getNEdges() ),
          nextid( inchan.getNEdges() ),
          p0( inchan.getNTri() ), p1( inchan.getNTri() ), p2( inchan.getNTri() ),
          e0( inchan.getNTri() ), e1( inchan.getNTri() ), e2( inchan.getNTri() ),
          t0( inchan.getNTri() ), t1( inchan.getNTri() ), t2( inchan.getNTri() )
{
   if( !inchan.getNodeInFile().good() ) GetKeyEntry( inchan );
   else GetFileEntry( inchan );
   cout << "tListInputData( inchan )" << endl;
}

template< class tSubNode >
void tListInputData< tSubNode >::
GetKeyEntry( tListIFStreams &inchan )                   //tListInputData
{
   int i;
   for( i=0; i < inchan.getNNodes(); i++ )
   {
      cout << "x y z edgid:" << endl;
      cin >> x[i] >> y[i] >> z[i] >> edgid[i];
   }
   for( i=0; i < inchan.getNEdges(); i++ )
   {
      cout << "orgid destid nextid" << endl;
      cin >> orgid[i] >> destid[i] >> nextid[i];
   }
   for( i=0; i< inchan.getNTri(); i++ )
   {
      cout << "nodeids (3), edgids (3), triangleids (3)" << endl;
      cin >> p0[i] >> p1[i] >> p2[i]
          >> e0[i] >> e1[i] >> e2[i]
          >> t0[i] >> t1[i] >> t2[i];
   }
}

template< class tSubNode >
void tListInputData< tSubNode >::
GetFileEntry( tListIFStreams &inchan )                  //tListInputData
{
   int i;
   for( i=0; i< inchan.getNNodes(); i++ )
       inchan.getNodeInFile() >> x[i] >> y[i] >> z[i]
                              >> edgid[i] >> boundflag[i];
   
   for( i=0; i< inchan.getNEdges(); i++ )
       inchan.getEdgeInFile() >> orgid[i] >> destid[i] >> nextid[i];
   
   for( i=0; i< inchan.getNTri(); i++ )
       inchan.getTriInFile() >> p0[i] >> p1[i] >> p2[i]
                             >> t0[i] >> t1[i] >> t2[i]
                             >> e0[i] >> e1[i] >> e2[i];
   
}

template< class tSubNode >
tListInputData< tSubNode >::
~tListInputData()                                                   //tInput
{
   cout << "    ~tListInputData()" << endl;
}

